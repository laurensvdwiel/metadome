from metadome.domain.wrappers.hmmer import FoundNoPfamHMMException,\
    FoundMoreThanOnePfamHMMException, convert_pfam_fasta_alignment_to_original_aligned_sequence,\
    map_sequence_to_aligned_sequence, convert_pfam_fasta_alignment_to_strict_fasta,\
    create_alignment_of_sequences_according_to_PFAM_HMM, interpret_hmm_alignment_file, convert_pfam_fasta_alignment_to_strict_sequence
from metadome.domain.data_generation.mapping.Protein2ProteinMapping import createAlignedSequenceMapping    
from metadome.domain.repositories import MappingRepository, SequenceRepository, InterproRepository, ProteinRepository, \
    GeneRepository, MetaDomainRepository
from metadome.default_settings import METADOMAIN_DIR, METADOMAIN_ALIGNMENT_FILE_NAME
from metadome.domain.models.entities.codon import Codon, MalformedCodonException
import numpy as np
import time
import logging

_log = logging.getLogger(__name__)

class MalformedMappingsForAlignedCodonsPosition(Exception):
    pass

def generate_pfam_alignments(pfam_id):
    """
    Generates a Pfam alignment for each sequence in 
    the human genome, which corresponds to pfam sequence
    """
    _log.info("Started creating an alignment of all '"+pfam_id+"' Pfam domains in the human genome")
    start_time = time.clock()
    
    # retrieve all domain occurrences for the domain_id
    domain_of_interest_occurrences = InterproRepository.get_domains_for_ext_domain_id(pfam_id)
    
    # First: retrieve all protein ids for this domain
    protein_ids = [int(y) for y  in np.unique([x.protein_id for x in domain_of_interest_occurrences])]
    
    # Retrieve all acs for these proteins
    protein_acs = ProteinRepository.retrieve_protein_ac_for_multiple_protein_ids(protein_ids)
    
    # Second, get all mappings for these proteins
    protein_mappings = MappingRepository.get_mappings_for_multiple_protein_ids(protein_ids)

    # Third: cut out the sequences from these mappings to Retrieve all the sequences of the domain of interest
    domain_of_interest_sequences = [{'sequence':SequenceRepository.get_aa_region(
        sequence=SequenceRepository.get_aa_sequence(mappings=protein_mappings[domain_occurrence.protein_id], skip_asterix_at_end=True), 
        region_start=domain_occurrence.uniprot_start, 
        region_stop=domain_occurrence.uniprot_stop), 
        'uniprot_ac':protein_acs[domain_occurrence.protein_id], 
        'start':domain_occurrence.uniprot_start,
        'stop':domain_occurrence.uniprot_stop} for domain_occurrence in domain_of_interest_occurrences]
    _log.debug("Starting HMM based alignment on for domain '"+pfam_id+"' for '"+str(len(domain_of_interest_occurrences))+"' occurrences across HG19")
    # Run the HMMERAlign algorithm based on the Pfam HMM
    try:
        create_alignment_of_sequences_according_to_PFAM_HMM(sequences=domain_of_interest_sequences, pfam_ac=pfam_id, target_directory=METADOMAIN_DIR, target_file_alignments=METADOMAIN_ALIGNMENT_FILE_NAME)
    except (FoundNoPfamHMMException, FoundMoreThanOnePfamHMMException) as e:
        _log.error(e)
        time_step = time.clock()
        _log.error("Prematurely stopped creating the '"+pfam_id+"' 'meta'-domain in "+str(time_step-start_time)+" seconds")
        return None
    _log.debug("Finished HMM based alignment on for domain '"+pfam_id+"'")
    
    time_step = time.clock()
    _log.info("Finished the mappings for '"+str(len(domain_of_interest_occurrences)) +"' '"+pfam_id+"' domain occurrences in "+str(time_step-start_time)+" seconds")   

def retrieve_pfam_aligned_codons(pfam_id, genome_build):
    """
    Queries the database to generate a list of dictionaries (meta_codons_per_consensus_pos)
    wherein all aligned codons per domain consensus positions are located
    Also provides the consensus_length of the domain and the n_instances
    """
    _log.info("Started a meta-domain based on the alignment of all '" + pfam_id + "' Pfam domains in the human genome")
    start_time = time.perf_counter()

    # the consensus length
    consensus_length = 0
    # the meta_domain that is to be returned
    meta_codons_per_consensus_pos = {}
    # the mapping of the protein {protein_id: {protein_posistion: consensus_position}}
    consensus_pos_per_protein = {}
    # the amount of domain occurrences found
    n_instances = 0

    # Try and retrieve the meta-domain position mappings
    meta_domain_results = MetaDomainRepository.retrieve_all_mappings_for_meta_domain(pfam_id, genome_build)

    if not len(meta_domain_results) == 0:
        # Calculate consensus length and count unique interpro domains
        consensus_positions = set()
        unique_interpro_ids = set()

        # Dict for gene ids to transcription ids
        gene_ids_to_transcription_ids = {}

        # Dict for protein_id to uniprot_ac
        protein_id_to_uniprot_ac = {}

        # Dict for mappings per uniprot ac, all mappings for a multiple Protein objects as {protein_id: [ Mapping ]}
        protein_mappings = {}

        # Dict to convert uniprot_acs to protein ids
        uniprot_acs_to_ids = {}

        # first pass over results to retrieve consensus length and unique interpro ids and
        for mapping, gene, mdm, mdp, protein in meta_domain_results:
            consensus_positions.add(mdp.consensus_position)
            if mdm.interpro_id:
                unique_interpro_ids.add(mdm.interpro_id)

            # retrieve the transcription id for this gene_id
            if not mapping.gene_id in gene_ids_to_transcription_ids.keys():
                gene_ids_to_transcription_ids[mapping.gene_id] = gene.gencode_transcription_id

            # retrieve the uniprot_ac for this protein_id
            if not mapping.protein_id in protein_id_to_uniprot_ac.keys():
                protein_id_to_uniprot_ac[mapping.protein_id] = protein.uniprot_ac

            # retrieve the position in the domain consensus
            domain_consensus_pos = mdp.consensus_position

            # retrieve the uniprot position of the meta-domain position
            uniprot_pos = mapping.uniprot_position

            # Add the consensus pos to the protein
            if not protein.uniprot_ac in consensus_pos_per_protein.keys():
                consensus_pos_per_protein[protein.uniprot_ac] = {}
            if not uniprot_pos in consensus_pos_per_protein[protein.uniprot_ac].keys():
                consensus_pos_per_protein[protein.uniprot_ac][uniprot_pos] = []
            consensus_pos_per_protein[protein.uniprot_ac][uniprot_pos].append(domain_consensus_pos)

            # Add the mapping to the protein_mappings
            if not mapping.protein_id in protein_mappings.keys():
                protein_mappings[mapping.protein_id] = []
            protein_mappings[mapping.protein_id].append(mapping)

            # Add the uniprot_ac to protein_id mapping
            if not protein.uniprot_ac in uniprot_acs_to_ids.keys():
                uniprot_acs_to_ids[protein.uniprot_ac] = mapping.protein_id

        # the consensus length
        consensus_length = max(consensus_positions) if consensus_positions else 0

        # the amount of domain occurrences found
        n_instances = len(unique_interpro_ids)

        # create all aligned codons
        for uniprot_ac in consensus_pos_per_protein.keys():
            for uniprot_pos in consensus_pos_per_protein[uniprot_ac].keys():
                for domain_consensus_pos in consensus_pos_per_protein[uniprot_ac][uniprot_pos]:
                    # Retrieve the mapping for the corresponding uniprot_position
                    mappings_for_uniprot_pos = [x for x in protein_mappings[uniprot_acs_to_ids[uniprot_ac]] if x.uniprot_position == uniprot_pos]

                    # Seperate the mappings per gene_id
                    mapping_per_gene_id = {}
                    for mapping in mappings_for_uniprot_pos:
                        if not mapping.gene_id in mapping_per_gene_id.keys():
                            mapping_per_gene_id[mapping.gene_id] = []
                        mapping_per_gene_id[mapping.gene_id].append(mapping)

                    for gene_id in mapping_per_gene_id.keys():
                        # Obtain the mappings for this position
                        mappings = mapping_per_gene_id[gene_id]

                        # Group mappings by genomic position to handle duplicate entries
                        position_mappings = {}
                        for mapping in mappings:
                            position = mapping.chromosome_position
                            if position not in position_mappings:
                                position_mappings[position] = mapping

                        unique_mappings = list(position_mappings.values())

                        try:
                            # create a codon
                            codon = Codon.initializeFromMapping(unique_mappings, gene_ids_to_transcription_ids[gene_id], uniprot_ac)

                            # Add the codon to the consensus positions
                            if not domain_consensus_pos in meta_codons_per_consensus_pos.keys():
                                meta_codons_per_consensus_pos[domain_consensus_pos] = []

                            meta_codons_per_consensus_pos[domain_consensus_pos].append(codon)
                        except MalformedCodonException as e:
                            logging.info("All mappings used to create the codon: "+str([str(m) for m in unique_mappings]))
                            raise MalformedMappingsForAlignedCodonsPosition(
                                "Encountered a malformed codon mapping for domain '"
                                + str(pfam_id) + "' in gene '" + str(gene_id)
                                + "', at amino_acid_position '" + str(uniprot_pos)
                                + "':" + str(e))
    else:
        _log.info("No mappings found for domain '" + str(pfam_id) + "' in genome build '" + str(genome_build) + "'")

    time_step = time.perf_counter()
    _log.info("Finished the alignment of mappings for '" + str(
        n_instances) + "' instances '" + pfam_id + "' domain occurrences in " + str(
        time_step - start_time) + " seconds")
    return meta_codons_per_consensus_pos, consensus_length, n_instances

def generate_pfam_aligned_codons(pfam_id): #@todo mark for deletion
    """
    Utilizes the Pfam HMMER alignment files to generate a list of dictionaries (meta_codons_per_consensus_pos)
    wherein all aligned codons per domain consensus positions are located
    Also provides the consensus_length of the domain and the n_instances
    """
    _log.info("Started a meta-domain based on the alignment of all '"+pfam_id+"' Pfam domains in the human genome")
    start_time = time.perf_counter()
     
    # the consensus length 
    consensus_length = 0
    # the meta_domain that is to be returned
    meta_codons_per_consensus_pos = {}
    # the mapping of the protein {protein_id: {protein_posistion: consensus_position}}
    consensus_pos_per_protein = {}
    # the amount of domain occurrences found
    n_instances = 0 
    
    # retrieve the alignment
    hmmeralign_output = interpret_hmm_alignment_file(METADOMAIN_DIR+pfam_id+'/'+METADOMAIN_ALIGNMENT_FILE_NAME) 
    if not len (hmmeralign_output) == 0:
        #update the consensus length
        consensus_length = len(hmmeralign_output['consensus']['sequence'])
        
        # update the number of instances
        n_instances = len(hmmeralign_output['alignments'])
        _log.debug("Creating the alignment of mappings for '"+str(n_instances) +"' '"+pfam_id+"' domain occurrences based on the HMM alignment to consensus and original domain sequence")
        
        # ensure we can map consensus residues back to consensus positions
        hmmeralign_output['consensus']['aligned_sequence'] = convert_pfam_fasta_alignment_to_original_aligned_sequence(hmmeralign_output['consensus']['alignment'])
        hmmeralign_output['consensus']['mapping_consensus_alignment_to_positions'] = map_sequence_to_aligned_sequence(hmmeralign_output['consensus']['sequence'], hmmeralign_output['consensus']['aligned_sequence'])
         
        # create mappings between domain occurrences and the domain consensus sequence
        for _alignment in hmmeralign_output['alignments']:
            # retrieve current aligned domain
                
            # Create a mapping from the aligned domain sequence to the domain sequence
            aligned_sequence = convert_pfam_fasta_alignment_to_original_aligned_sequence(_alignment['alignment'])
            original_sequence = convert_pfam_fasta_alignment_to_strict_sequence(aligned_sequence)
            mapping_domain_alignment_to_sequence_positions = map_sequence_to_aligned_sequence(original_sequence, aligned_sequence)
                
            # Generate the strict sequence for this domain; leaving only residues that were aligned to the domain consensus
            strict_aligned_sequence = convert_pfam_fasta_alignment_to_strict_fasta(_alignment['alignment'])
                
            # create the mapping between the strict alignments and the original consensus sequence
            mapping_aligned_domain_to_domain_consensus = createAlignedSequenceMapping(strict_aligned_sequence, hmmeralign_output['consensus']['aligned_sequence'], False)
                
            # create a list of mapping positions that includes insertions
            mapping_positions = list(mapping_domain_alignment_to_sequence_positions.keys()) + list(set(mapping_aligned_domain_to_domain_consensus.keys()) - set(mapping_domain_alignment_to_sequence_positions.keys()))
                
            # Second add each aligned residue mapping
            for mapping_pos in sorted(mapping_positions):
                # retrieve the residue at the consensus position and the residue at the domain position
                consensus_domain_residue = hmmeralign_output['consensus']['aligned_sequence'][mapping_pos]
                    
                if consensus_domain_residue == '-':
                    # Set the default values for the insertion
                    continue
                else:
                    # retrieve the position in the domain consensus
                    domain_consensus_pos = hmmeralign_output['consensus']['mapping_consensus_alignment_to_positions'][mapping_pos]
                    
                # retrieve the position in the domain sequence
                ref_pos = mapping_domain_alignment_to_sequence_positions[mapping_pos]
                # convert the position in the domain sequence to the uniprot position and genomic position
                uniprot_pos = int(_alignment['start_pos']) + ref_pos -1
                 
                # Add the consensus pos to the protein
                if not _alignment['uniprot_ac'] in consensus_pos_per_protein.keys():
                    consensus_pos_per_protein[_alignment['uniprot_ac']] = {}
                if not uniprot_pos in consensus_pos_per_protein[_alignment['uniprot_ac']].keys():
                    consensus_pos_per_protein[_alignment['uniprot_ac']][uniprot_pos] = []
                consensus_pos_per_protein[_alignment['uniprot_ac']][uniprot_pos].append(domain_consensus_pos)                
                
        # now incorporate the alignment data into our domain model in form of mappings
        # First get the protein ids for the uniprot acs
        uniprot_acs_to_ids = ProteinRepository.retrieve_protein_id_for_multiple_protein_acs([x for x in consensus_pos_per_protein.keys()])
        protein_ids = [int(y) for y  in np.unique([x for x in uniprot_acs_to_ids.values()])]
        
        # Second, get all mappings for these proteins
        protein_mappings = MappingRepository.get_mappings_for_multiple_protein_ids(protein_ids)
        
        # retrieve all transcripts mapped to these protein_ids
        gene_ids = GeneRepository.retrieve_transcript_id_for_multiple_protein_ids(protein_ids)
        
        # create all aligned codons
        meta_codons_per_consensus_pos = {}
        for uniprot_ac in consensus_pos_per_protein.keys():
            for uniprot_pos in consensus_pos_per_protein[uniprot_ac].keys():
                for domain_consensus_pos in consensus_pos_per_protein[uniprot_ac][uniprot_pos]:
                    # Retrieve the mapping for the corresponding uniprot_position
                    mappings_for_uniprot_pos = [x for x in protein_mappings[uniprot_acs_to_ids[uniprot_ac]] if x.uniprot_position == uniprot_pos]
                     
                    # Seperate the mappings per gene_id
                    mapping_per_gene_id = {}
                    for mapping in mappings_for_uniprot_pos:
                        if not mapping.gene_id in mapping_per_gene_id.keys():
                            mapping_per_gene_id[mapping.gene_id] = []
                        mapping_per_gene_id[mapping.gene_id].append(mapping)
                   
                    for gene_id in mapping_per_gene_id.keys():
                        # Obtain the mappings for this position
                        mappings = mapping_per_gene_id[gene_id]

                        try:
                            # create a codon
                            codon = Codon.initializeFromMapping(mappings, gene_ids[gene_id], uniprot_ac)
                            
                            # Add the codon to the consensus positions
                            if not domain_consensus_pos in meta_codons_per_consensus_pos.keys():
                                meta_codons_per_consensus_pos[domain_consensus_pos] = []
                            
                            meta_codons_per_consensus_pos[domain_consensus_pos].append(codon)
                        except MalformedCodonException as e:
                            raise MalformedMappingsForAlignedCodonsPosition("Encountered a malformed codon mapping for domain '"
                                                                         +str(pfam_id)+"' in gene '"+str(gene_id)
                                                                         +"', at amino_acid_position '"+str(uniprot_pos)
                                                                         +"':" + str(e))
   
    time_step = time.perf_counter()
    _log.info("Finished the alignment of mappings for '"+str(n_instances) +"' instances '"+pfam_id+"' domain occurrences in "+str(time_step-start_time)+" seconds")
    return meta_codons_per_consensus_pos, consensus_length, n_instances
