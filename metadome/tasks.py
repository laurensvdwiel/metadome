from metadome.domain.repositories import GeneRepository, RepositoryException
from metadome.domain.models.entities.gene_region import GeneRegion, GenomeBuild
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant
from metadome.domain.models.entities.meta_domain import MetaDomain,\
    UnsupportedMetaDomainIdentifier
from metadome.domain.services.computation.gene_region_computations import compute_tolerance_landscape
from metadome.domain.services.annotation.gene_region_annotators import GRCh37_annotateTranscriptWithClinvarData,\
    GRCh38_annotateTranscriptWithClinvarData
from metadome.domain.services.annotation.annotation import annotateSNVs,\
    convertNucleotide
from metadome.controllers.job import store_error, store_visualization
from metadome.domain.error import RecoverableError
import numpy as np

from celery import current_app as celery_app
from celery.signals import after_task_publish
from flask import current_app as flask_app

import traceback
import json
import os
import logging
from metadome.domain.models.gene import Strand

_log = logging.getLogger(__name__)



@after_task_publish.connect
def update_sent_state(sender=None, body=None, **kwargs):
    # the task may not exist if sent using `send_task` which
    # sends tasks by name, so fall back to the default result backend
    # if that is the case.
    task = celery_app.tasks.get(sender)
    backend = task.backend if task else celery_app.backend
    backend.store_result(kwargs['headers']['id'], None, "SENT")

    _log.debug("set task {} status to SENT".format(body['id']))

    # After this, celery tasks that do exist have status 'SENT'.


def get_task(job_name):
    """
    Get the task for the given input_type and output_type combination.
    If the combination is not allowed, a ValueError is raised.
    """
    _log.debug(get_celery_worker_status())

    task = None
    if job_name == 'create':
        task = create_prebuild_visualization
    elif job_name == 'retrieve':
        task = retrieve_prebuild_visualization
    else:
        raise ValueError("Unexpected input_type '{}'".format(job_name))

    _log.debug("Got task '{}'".format(task.__name__))
    return task

def get_celery_worker_status():
    """Used for retrieving celery status. Usefule for debug purposes"""
    ERROR_KEY = "ERROR"
    try:
        from celery.task.control import inspect
        insp = inspect()
        d = insp.stats()
        if not d:
            d = { ERROR_KEY: 'No running Celery workers were found.' }
    except IOError as e:
        from errno import errorcode
        msg = "Error connecting to the backend: " + str(e)
        if len(e.args) > 0 and errorcode.get(e.args[0]) == 'ECONNREFUSED':
            msg += ' Check that the RabbitMQ server is running.'
        d = { ERROR_KEY: msg }
    except ImportError as e:
        d = { ERROR_KEY: str(e)}
    return d

@celery_app.task(bind=True,
                 autoretry_for=(RecoverableError,),
                 retry_kwargs={'max_retries': 50})
def initialize_metadomain(self, domain_id):
    return MetaDomain.initializeFromDomainID(domain_id)

@celery_app.task(bind=True,
                 autoretry_for=(RecoverableError,),
                 retry_kwargs={'max_retries': 50})
def retrieve_prebuild_visualization(self, transcript_id):
    _log.info("Getting prebuild visualization for '{}'".format(transcript_id))

    # Determine path to the file and check that it exists.
    visualization_path = flask_app.config['PRE_BUILD_VISUALIZATION_DIR'] + transcript_id + '/' + flask_app.config['PRE_BUILD_VISUALIZATION_FILE_NAME']
    if not os.path.exists(visualization_path):
        raise RuntimeError("File not found: '{}'".format(visualization_path))

    # Unzip the file and return the contents
    _log.info("Reading '{}'".format(visualization_path))
    with open(visualization_path) as f:
        visualization_content = json.load(f)
    return visualization_content

@celery_app.task(bind=True,
                 autoretry_for=(RecoverableError,),
                 retry_kwargs={'max_retries': 50})
def create_prebuild_visualization(self, transcript_id, genome_build):
    try:
        _log.info("Attempting to create visualization for '{}'".format(transcript_id))

        result = analyse_transcript(transcript_id, genome_build)
        if 'error' in result.keys():
            _log.info("Something went wrong while trying to create visualization for '{}'".format(transcript_id))
            raise RuntimeError(result['error'])

        store_visualization(transcript_id, genome_build, result)
        return result
    except:
        store_error(transcript_id, genome_build, traceback.format_exc())
        raise

def retrieve_metadomain_annotation(transcript_id, protein_position, domain_positions):
    # first correct the protein_position
    protein_position -= 1
        
    domain_results = {}

    for domain_id in domain_positions.keys():
        # add new key to domain results
        domain_results[domain_id] = {}

        # create the values that are to be returned
        normal_variants = []
        pathogenic_variants = []
        alignment_depth = 0

        # retrieve the metadomain
        meta_domain = MetaDomain.initializeFromDomainID(domain_id)

        # retrieve the codon
        current_codon = meta_domain.get_codon_for_transcript_and_position(transcript_id, protein_position)

        for consensus_position in domain_positions[domain_id]:
            # first correct the consensus_position
            consensus_position -= 1

            # Retrieve the meta codons for this position
            meta_codons = meta_domain.get_codons_aligned_to_consensus_position(consensus_position)
            alignment_depth += len(meta_codons)

            # Retrieve the meta SNVs for this position
            meta_snvs = meta_domain.get_annotated_SNVs_for_consensus_position(consensus_position)

            # Retrieve the matching gene names for the transcripts
            transcript_ids = [meta_snvs[meta_snv_repr][0]['gencode_transcription_id'] for meta_snv_repr in meta_snvs.keys()]
            transcripts_to_gene = GeneRepository.retrieve_gene_names_for_multiple_transcript_ids(transcript_ids)

            # iterate over meta_codons and add to metadom_entry
            for meta_snv_repr in meta_snvs.keys():
                if not current_codon.unique_str_representation() in meta_snv_repr:
                    # unique variant at homologous position, can just take the first from the list
                    meta_snv = meta_snvs[meta_snv_repr][0]

                    # initiate the SNV variant
                    snv_variant = SingleNucleotideVariant.initializeFromDict(meta_snv)

                    # Convert to the original nucleotide
                    if snv_variant.strand == Strand.minus:
                        snv_variant.ref_nucleotide = convertNucleotide(snv_variant.ref_nucleotide)
                        snv_variant.alt_nucleotide = convertNucleotide(snv_variant.alt_nucleotide)

                    # start the variant entry and add the codon based information
                    variant_entry = snv_variant.toCodonJson()

                    # Add the gene name
                    variant_entry['gene_name'] = transcripts_to_gene[snv_variant.gencode_transcription_id]

                    # Add the variant specific information
                    if meta_snv['variant_source'] == 'gnomAD':
                        # convert the variant to the expected format
                        gnomad_json = snv_variant.toGnommADJson(allele_number=meta_snv['allele_number'], allele_count=meta_snv['allele_count'])
                        for key in gnomad_json.keys():
                            variant_entry[key] = gnomad_json[key]

                        # append to the list of variants
                        normal_variants.append(variant_entry)
                    elif meta_snv['variant_source'] == 'ClinVar':
                        # convert the variant to the expected format
                        clinvar_json = snv_variant.initializeFromDict(meta_snv).toClinVarJson(ClinVar_id=meta_snv['clinvar_ID'])
                        for key in clinvar_json.keys():
                            variant_entry[key] = clinvar_json[key]

                        # append to the list of variants
                        pathogenic_variants.append(variant_entry)

        domain_results[domain_id]["pathogenic_variants"] = pathogenic_variants
        domain_results[domain_id]["normal_variants"] = normal_variants
        domain_results[domain_id]["alignment_depth"] = alignment_depth

    return domain_results

def analyse_transcript(transcript_id, genome_build):
    # Retrieve the gene from the database
    try:
        gene = GeneRepository.retrieve_gene(transcript_id, genome_build)
    except RepositoryException as e:
        return {'error': 'No gene region could be build for transcript {}, genome_build {}, reason: {}'.format(transcript_id, genome_build, e)}

    # build the gene region
    gene_region = GeneRegion(gene)
    
    # Retrieve the refseq ids    
    refseq_ids = retrieve_refseq_identifiers_for_transcript(transcript_id)

    if not gene_region is None:
        # generate the positional annotation for this gene by first computing the tolerance landscape
        region_positional_annotation = compute_tolerance_landscape(gene_region, flask_app.config['SLIDING_WINDOW_SIZE'], flask_app.config['ALLELE_FREQUENCY_CUTOFF'])

        # Annotate Pfam domains
        Pfam_domains = []
        meta_domains = {}
        for domain in gene_region.interpro_domains:
            if domain.ext_db_id.startswith('PF'):
                # we have a Pfam domain
                pfam_domain = {}
                pfam_domain["ID"] = domain.ext_db_id
                pfam_domain["Name"] = domain.region_name
                pfam_domain["start"] = domain.uniprot_start
                pfam_domain["stop"] = domain.uniprot_stop

                try:
                    if not pfam_domain['ID'] in meta_domains.keys():
                        # construct a meta-domain if possible
                        temp_meta_domain = MetaDomain.initializeFromDomainID(domain.ext_db_id, gene_region.genome_build)

                        # Ensure there are enough instances to actually perform the metadomain trick
                        if temp_meta_domain.n_instances < 2:
                            pfam_domain["metadomain"] = False
                            meta_domains[pfam_domain['ID']] = None
                        else:
                            pfam_domain["metadomain"] = True
                            meta_domains[pfam_domain['ID']] = temp_meta_domain
                            pfam_domain['meta_domain_alignment_depth'] = temp_meta_domain.get_max_alignment_depth()
                    else:
                        pfam_domain["metadomain"] = not(meta_domains[pfam_domain['ID']] is None)
                except UnsupportedMetaDomainIdentifier as e:
                    _log.error(str(e))
                    # meta domain is not possible
                    meta_domains[pfam_domain['ID']] = None

                # Add the domain to the domain list
                Pfam_domains.append(pfam_domain)

        # Based on genome build, select the correct gnomAD annotation function
        if gene_region.genome_build == GenomeBuild.GRCh37:
            clinvar_annotation_function = GRCh37_annotateTranscriptWithClinvarData
        elif gene_region.genome_build == GenomeBuild.GRCh38:
            clinvar_annotation_function = GRCh38_annotateTranscriptWithClinvarData
        else:
            raise Exception("Could not determine ClinVar annotation function based on genome build '" + str(
                gene_region.genome_build.value) + "'")

        # Annotate the clinvar variants for the current gene
        ClinVar_annotation = annotateSNVs(clinvar_annotation_function,
                                         mappings_per_chr_pos=gene_region.retrieve_mappings_per_chromosome(),
                                         strand=gene_region.strand,
                                         chromosome=gene_region.chr,
                                         regions=gene_region.regions)

        # retrieve the mappings per chromosome position
        _mappings_per_chromosome = gene_region.retrieve_mappings_per_chromosome()

        for chrom_pos in ClinVar_annotation.keys():
            for variant in ClinVar_annotation[chrom_pos]:
                protein_pos = _mappings_per_chromosome[chrom_pos]['amino_acid_position']

                if not 'ClinVar' in region_positional_annotation[protein_pos].keys():
                    region_positional_annotation[protein_pos]['ClinVar'] = []

                codon = gene_region.retrieve_codon_for_protein_position(protein_pos)

                # create new entry for this variant
                variant_entry = SingleNucleotideVariant.initializeFromVariant(_codon=codon, _chr_position=chrom_pos, _alt_nucleotide=variant['ALT'], _variant_source='ClinVar').toClinVarJson(ClinVar_id=variant['ID'])
                region_positional_annotation[protein_pos]['ClinVar'].append(variant_entry)

        # annotate the positions further
        for d in region_positional_annotation:
            # retrieve the position as is in the database
            db_position = d['protein_pos']

            # update the positions to abide the users' expectation (start at 1, not zero)
            d.update((k, v+1) for k, v in d.items() if k == "protein_pos")

            # add domain and meta domain information per position
            d['domains'] = {}
            for domain in Pfam_domains:
                if d['protein_pos'] >= domain["start"] and d['protein_pos'] <= domain["stop"]:
                    # add the domain id for this position
                    d['domains'][domain['ID']] = None
                    if not meta_domains[domain['ID']] is None:
                        # retrieve the context for this protein
                        consensus_positions = meta_domains[domain['ID']].get_consensus_positions_for_uniprot_position(uniprot_ac=gene_region.uniprot_ac, uniprot_position=db_position)

                        if domain["metadomain"] and len(consensus_positions)>0:
                            d['domains'][domain['ID']] = create_meta_domain_entry(gene_region, meta_domains[domain['ID']], consensus_positions, db_position)
        result = {"transcript_id":transcript_id, "refseq_ids":refseq_ids, "protein_ac":gene_region.uniprot_ac, "gene_name":gene_region.gene_name, "positional_annotation":region_positional_annotation, "domains":Pfam_domains}
    else:
        result = {'error': 'No gene region could be build for transcript '+str(transcript_id)}

    return result

def create_meta_domain_entry(gene_region, metadomain, consensus_positions, protein_pos):
    metadom_entry = {}
    # update the consensus positions to abide the users' expectation (start at 1, not zero)
    metadom_entry['consensus_pos'] = [int(x) for x in np.array(consensus_positions)+1]
    metadom_entry['normal_missense_variant_count'] = 0
    metadom_entry['normal_variant_count'] = 0
    metadom_entry['pathogenic_missense_variant_count'] = 0
    metadom_entry['pathogenic_variant_count'] = 0
    for consensus_position in consensus_positions:
        # Retrieve the meta codons for this position
        meta_snvs = metadomain.get_annotated_SNVs_for_consensus_position(consensus_position)
        current_codon = gene_region.retrieve_codon_for_protein_position(protein_pos)
        # iterate over meta_codons and add to metadom_entry
        for meta_snv_repr in meta_snvs.keys():
            if not current_codon.unique_str_representation() in meta_snv_repr:
                # unique variant at homologous position, can just take the first from the list
                meta_snv = meta_snvs[meta_snv_repr][0]
                if meta_snv['variant_source'] == 'gnomAD':
                    # Update the variant counts
                    metadom_entry['normal_variant_count'] += 1
                    if meta_snv['variant_type'] == 'missense': metadom_entry['normal_missense_variant_count'] += 1
                elif meta_snv['variant_source'] == 'ClinVar':
                    # Update the variant counts
                    metadom_entry['pathogenic_variant_count'] += 1
                    if meta_snv['variant_type'] == 'missense': metadom_entry['pathogenic_missense_variant_count'] += 1
    return metadom_entry
