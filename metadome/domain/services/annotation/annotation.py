from metadome.domain.models.gene import Strand

import logging

_log = logging.getLogger(__name__)

class NucleotideConversionException(Exception):
    pass    

def convertNucleotide(nucleotide):
    if nucleotide == 'A': return 'T'
    elif nucleotide == 'T': return 'A'
    elif nucleotide == 'C': return 'G'
    elif nucleotide == 'G': return 'C'
    else: raise NucleotideConversionException("Received nucleotide '"+str(nucleotide)+"', which could not be converted")

def annotateSNVs(annotateTranscriptFunction, mappings_per_chr_pos, strand, chromosome, regions):
    SNV_correct = 0
    SNV_incorrect = 0
    
    Annotation_mapping = {}

    # annotate the region with provided data
    annotation_data = annotateTranscriptFunction(chromosome, regions)

    try:
        for annotation in annotation_data:
            # Check if we are dealing with a SNV
            if annotation['REF'] is None or len(annotation['REF']) != 1: continue
            if annotation['ALT'] is None or len(annotation['ALT']) != 1: continue

            # retrieve the reference
            ref = str(annotation['REF'])
            # ensure we have correct strand, otherwise convert
            if strand == Strand.minus: ref = convertNucleotide(annotation['REF'])

            # Check if the reference in the annotation_date is correct, compared to our transcript
            if mappings_per_chr_pos[annotation['POS']]['base_pair'] == ref:
                SNV_correct = SNV_correct + 1
            else:
                SNV_incorrect = SNV_incorrect + 1

            # retrieve the alternate variation
            alt = annotation['ALT']

            # convert to string
            alt = str(alt)

            # ensure we have correct strand, otherwise convert
            if strand == Strand.minus: alt = convertNucleotide(alt)

            # construct SNV entry
            SNV_entry = {key:annotation['INFO'][key] for key in annotation['INFO'].keys()}
            SNV_entry['REF'] = ref
            SNV_entry['ALT'] = alt
            SNV_entry['CHROM'] = annotation['CHROM']
            SNV_entry['POS'] = annotation['POS']

            # Genome key
            genomic_pos = int(annotation['POS'])

            # Annotate
            if genomic_pos not in Annotation_mapping.keys():
                Annotation_mapping[genomic_pos] = []

            Annotation_mapping[genomic_pos].append(SNV_entry)

        if SNV_incorrect > 0:
            _log.error("SNV ERROR: "+str(SNV_incorrect)+' incorrect vs '+str(SNV_correct)+' correct')
    except Exception as ex:
        raise Exception("Caught exception {}, with function {} and parameters chromosome {} and regions {}".format(str(ex), annotateTranscriptFunction.__name__, chromosome, regions))

    return Annotation_mapping
        