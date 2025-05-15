from metadome.database import db
from metadome.domain.models.protein import Protein
from metadome.domain.models.interpro import Interpro
from metadome.domain.repositories import MappingRepository, SequenceRepository
from metadome.domain.services.multi_threading import CalculateNumberOfActiveThreads
from metadome.domain.data_generation.mapping.mapping_generator import generate_gene_to_swissprot_mapping
from metadome.domain.infrastructure import add_gene_mapping_to_database,\
    filter_gene_names_present_in_database
from metadome.domain.wrappers.gencode import retrieve_all_protein_coding_gene_names
from metadome.domain.wrappers.interpro import retrieve_interpro_entries
from joblib.parallel import Parallel, delayed

import logging

_log = logging.getLogger(__name__)

def create_db():
    # initialize custom logging framework
    _log.info("Starting creation of new database")

    # the genes that are to be checked
    genes_of_interest = retrieve_all_protein_coding_gene_names()
         
    # (re-) construct the mapping database  => GENE2PROTEIN_MAPPING_DB
    generate_mappings_for_genes(genes_of_interest, batch_size=10, use_parallel=True)
      
    # annotate all the proteins with interpro_domains
    for protein in Protein.query.filter(Protein.evaluated_interpro_domains == False).all():
        # generate all pfam domain to swissprot mappings
        annotate_interpro_domains_to_proteins(protein)
    
def generate_mappings_for_genes(genes_of_interest, batch_size, use_parallel):
    # filter gene names alreay present in the database
    genes_of_interest = filter_gene_names_present_in_database(genes_of_interest)
    
    # Create batches
    genes_of_interest_batches = [genes_of_interest[i:i+batch_size] for i in range(0, len(genes_of_interest), batch_size)]
    n_batches = len(genes_of_interest_batches)
    n_genes = len(genes_of_interest)
     
    # Annotate the genes in batches
    _log.info("Starting the mapping of batched analysis of '"+str(n_genes)+"' over '"+str(n_batches)+"' batches")
    succeeded_genes = 0
    for batch_counter, gene_batch in enumerate(genes_of_interest_batches):
        _log.info("Starting batch '"+str(batch_counter+1)+"' out of '"+str(n_batches)+"', with '"+str(len(gene_batch))+"' genes")
     
        gene_mappings = []
        if use_parallel:
            gene_mappings = Parallel(n_jobs=CalculateNumberOfActiveThreads(batch_size))(delayed(generate_gene_to_swissprot_mapping)(gene) for gene in gene_batch)            
        else:
            gene_mappings = [generate_gene_to_swissprot_mapping(gene) for gene in gene_batch]

        # add the batches to the database
        for gene_mapping in gene_mappings:
            add_gene_mapping_to_database(gene_mapping)
            succeeded_genes +=1
         
        _log.info("Finished batch '"+str(batch_counter+1)+"' out of '"+str(n_batches)+"'")
    _log.info("Finished the mapping of batched analysis of '"+str(n_genes)+"' over '"+str(n_batches)+"' batches, resulting in '"+str(succeeded_genes)+"' successful gene mappings")

def annotate_interpro_domains_to_proteins(protein):
    _log.info("Protein "+str(protein.uniprot_ac)+" will now be annotated by Interpro...")
    # Annotate the interpro ids
    mappings = MappingRepository.get_mappings_for_protein(protein)
    aa_sequence = SequenceRepository.get_aa_sequence(mappings, skip_asterix_at_end=True)
    
    # Query the sequence to interpro
    interpro_results = retrieve_interpro_entries(protein.uniprot_ac, aa_sequence)
    
    # save the results to the database
    with db.session.no_autoflush:
        interpro_domains = []
        for interpro_result in interpro_results:
            if interpro_result['interpro_id'] is None and interpro_result['region_name'] == '':
                # skip non-informative results
                continue
                            
            # create a new interpro domain
            interpro_domain = Interpro(_interpro_id=interpro_result['interpro_id'],\
                                  _ext_db_id=interpro_result['ext_db_id'],\
                                  _region_name=interpro_result['region_name'],\
                                  _start_pos=interpro_result['start_pos'],\
                                  _end_pos=interpro_result['end_pos'])
             
            # Solve the required foreign key
            protein.interpro_domains.append(interpro_domain)
            
            # Add the interpro_domain to the database
            interpro_domains.append(interpro_domain)
        
        # Add the interpro domains to the database
        db.session.add_all(interpro_domains)
        
        # Adjust the value of the protein that it has already been evaluated
        protein.evaluated_interpro_domains = True
    
    _log.info("Protein "+str(protein.uniprot_ac)+" was annotated with '"+str(len(interpro_results))+"' interpro domains.")
    
    # Commit this session
    db.session.commit()