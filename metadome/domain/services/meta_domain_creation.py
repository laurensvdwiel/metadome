# TODO(dead-data-generation): unreachable at runtime.
#  Reached only through database_creation and install.py, which no longer exists.
from metadome.domain.repositories import InterproRepository
from metadome.domain.data_generation.mapping.meta_domain_mapping import generate_pfam_alignments

import logging
import os
from metadome.default_settings import METADOMAIN_DIR

_log = logging.getLogger(__name__)

def create_metadomains(reconstruct):
    # initialize custom logging framework
    _log.info("Starting creation of all meta-domain alignments")
 
    # the genes that are to be checked
    domains_of_interest = [x for x in InterproRepository.get_all_Pfam_identifiers_suitable_for_metadomains()]
     
    # check for which domains there have already been made alignments
    finished_meta_domains = os.listdir(METADOMAIN_DIR)
     
    # (re-) construct the meta domain alignments    
    domains_processed = 1
    domains_already_processed = 0
 
    for domain_id in domains_of_interest:
        if domain_id in finished_meta_domains and not reconstruct:
            domains_already_processed+=1
            continue
         
        _log.info("Starting creation of alignment for '"+str(domains_processed)+"' out of '"+str(len(domains_of_interest))+"'")
        generate_pfam_alignments(domain_id)
        _log.info("Finished creation of alignment for '"+str(domains_processed)+"' out of '"+str(len(domains_of_interest))+"'")
        domains_processed+=1
         
    _log.info("Finished the creation of meta-domain alignments of '"+str(len(domains_of_interest))+"' domains, resulting in '"+str(domains_processed)+"' successful meta-domain alignments, previously processed (='"+str(domains_already_processed)+"')")