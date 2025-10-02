import csv
import logging
import traceback
import enum 
import argparse
import os
import sys
import gzip
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.exc import IntegrityError
from metadome.default_settings import SQLALCHEMY_DATABASE_URI

# --- Flask App and SQLAlchemy Setup ---
# This would typically be in your app.py or a config file
flask_app = Flask(__name__)

# IMPORTANT: Configure your actual database URI
flask_app.config['SQLALCHEMY_DATABASE_URI'] = SQLALCHEMY_DATABASE_URI #f"postgresql://{DB_USER}:{DB_PWD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
flask_app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # Optional: silence a warning

db = SQLAlchemy() # Initialize your db object globally or pass around

# --- Model Definitions (gene.py, protein.py, etc.) ---
# (These files would import 'db' from wherever you defined it above)
from metadome.database import db
from metadome.domain.models.gene import Gene
from metadome.domain.models.protein import Protein, ProteinSource 
from metadome.domain.models.mapping import Mapping
from metadome.domain.models.interpro import Interpro
from metadome.domain.models.meta_domain_mapping import MetaDomainMapping

# --- Helper Functions (safe_int_conversion, safe_bool_conversion, get_cleaned_str) ---
# (These remain unchanged)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def safe_int_conversion(value_str, default=None):
    if value_str is None or value_str.strip() == '' or value_str.lower() == 'na' or value_str.lower() == '-':
        return default
    try:
        return int(value_str)
    except ValueError:
        logging.warning(f"Could not convert '{value_str}' to int, using default: {default}")
        return default

def safe_bool_conversion(value_str, default=False):
    if value_str is None or value_str.strip() == '':
        return default
    val_lower = value_str.lower()
    if val_lower in ['true', 't', '1', 'yes', 'y']:
        return True
    if val_lower in ['false', 'f', '0', 'no', 'n', 'na', '-']:
        return False
    logging.warning(f"Could not convert '{value_str}' to bool, using default: {default}")
    return default

def get_cleaned_str(row, key, default=None):
    val = row.get(key)
    if val is None: return default
    val_stripped = val.strip()
    return val_stripped if val_stripped else default

# --- Batch Load Function (uses the passed sqlalchemy_session, which will be db.session) ---
def batch_load_data(csv_filepath, sqlalchemy_session, batch_size=1000):
    protein_cache = {}
    gene_cache = {}
    interpro_cache = {}
    # The rest of the batch_load_data function uses 'sqlalchemy_session'
    #      for all database operations like .query, .add, .commit, .rollback
    try:
        with gzip.open(csv_filepath, mode='rt') as infile:
            reader = csv.DictReader(infile)
            current_batch_count = 0
            total_processed_count = 0

            for i, row in enumerate(reader):
                line_num = i + 2 
                total_processed_count += 1
                current_protein = None
                current_gene = None

                try:
                    # 1. Process Protein
                    uniprot_ac = get_cleaned_str(row, 'uniprot_ac')
                    if uniprot_ac:
                        if uniprot_ac in protein_cache:
                            current_protein = protein_cache[uniprot_ac]
                        else:
                            current_protein = sqlalchemy_session.query(Protein).filter_by(uniprot_ac=uniprot_ac).one_or_none()
                            if current_protein is None:
                                uniprot_name_val = get_cleaned_str(row, 'uniprot_name')
                                source_val_str = get_cleaned_str(row, 'source').lower()
                                eval_interpro_val = safe_bool_conversion(row.get('evaluated_interpro_domains'), False)
                                if not source_val_str:
                                    logging.warning(f"L{line_num}: Missing source for Protein '{uniprot_ac}'. Skipping Protein.")
                                else:
                                    current_protein = Protein(_uniprot_ac=uniprot_ac, 
                                                              _uniprot_name=uniprot_name_val, 
                                                              _source=source_val_str)
                                    current_protein.evaluated_interpro_domains = eval_interpro_val
                                    sqlalchemy_session.add(current_protein)
                            protein_cache[uniprot_ac] = current_protein
                    
                    # 2. Process Gene
                    gencode_tr_id = get_cleaned_str(row, 'gencode_transcription_id')
                    gencode_tr_version = get_cleaned_str(row, 'GENCODE_version')
                    current_gene_key = gencode_tr_id+"_"+gencode_tr_version
                    if gencode_tr_id:
                        if gencode_tr_id in gene_cache:
                            current_gene = gene_cache[current_gene_key]
                        else:
                            current_gene = sqlalchemy_session.query(Gene).filter_by(gencode_transcription_id=gencode_tr_id, gencode_version=gencode_tr_version).one_or_none()
                            if current_gene is None:
                                strand_str = get_cleaned_str(row, 'strand')
                                gene_name_val = get_cleaned_str(row, 'gene_name')
                                gencode_transl_name_val = get_cleaned_str(row, 'gencode_translation_name')
                                if not strand_str or strand_str not in ['+', '-']:
                                    logging.warning(f"L{line_num}: Invalid strand for Gene '{gencode_tr_id}'. Skipping Gene.")
                                elif not gencode_transl_name_val:
                                    logging.warning(f"L{line_num}: Missing gencode_translation_name for Gene '{gencode_tr_id}'. Skipping Gene.")
                                else:
                                    current_gene = Gene(
                                        _strand=strand_str, _gene_name=gene_name_val,
                                        _gencode_transcription_id=gencode_tr_id, _gencode_translation_name=gencode_transl_name_val,
                                        _gencode_gene_id=get_cleaned_str(row, 'gencode_gene_id'),
                                        _gencode_version=gencode_tr_version,
                                        _gencode_basic=safe_bool_conversion(row.get('GencodeBasic'), False),
                                        _genome_build=get_cleaned_str(row, 'genome_build'),
                                        _refseq_transcript_id=get_cleaned_str(row, 'refseq_transcript_id'),
                                        _havana_gene_id=get_cleaned_str(row, 'havana_gene_id'),
                                        _havana_translation_id=get_cleaned_str(row, 'havana_translation_id'),
                                        _mane_transcript_type=get_cleaned_str(row, 'MANE'),
                                        _sequence_length=safe_int_conversion(row.get('sequence_length')))
                                    if current_protein: current_gene.protein = current_protein
                                    sqlalchemy_session.add(current_gene)
                            gene_cache[current_gene_key] = current_gene
                    
                    # 3. Process Interpro
                    ext_db_id_val = get_cleaned_str(row, 'ext_db_id')
                    ext_db_version_val = get_cleaned_str(row, 'PFAM_version')
                    uniprot_start_val = safe_int_conversion(row.get('uniprot_start'))
                    uniprot_stop_val = safe_int_conversion(row.get('uniprot_stop'))
                    current_interpro = None
                    if current_protein is not None and ext_db_id_val is not None and ext_db_version_val is not None and uniprot_start_val is not None and uniprot_stop_val is not None:
                        interpro_cache_key = (current_protein.uniprot_ac, ext_db_id_val, ext_db_version_val, uniprot_start_val, uniprot_stop_val)
                        if interpro_cache_key not in interpro_cache:
                            current_interpro = None
                            if current_protein.id:
                                current_interpro = sqlalchemy_session.query(Interpro).filter_by(
                                    protein_id=current_protein.id, ext_db_id=ext_db_id_val, ext_db_version=ext_db_version_val,
                                    uniprot_start=uniprot_start_val, uniprot_stop=uniprot_stop_val).one_or_none()
                            if not current_interpro:
                                current_interpro = Interpro(
                                    _ext_db_id=ext_db_id_val, _ext_db_version=ext_db_version_val, _start_pos=uniprot_start_val, _end_pos=uniprot_stop_val,
                                    _interpro_id=get_cleaned_str(row, 'interpro_id'), _region_name=get_cleaned_str(row, 'region_name'))
                                current_interpro.protein = current_protein
                                sqlalchemy_session.add(current_interpro)
                                interpro_cache[interpro_cache_key] = current_interpro
                            else:
                                interpro_cache[interpro_cache_key] = current_interpro
                        else:
                            current_interpro = interpro_cache[interpro_cache_key]
                    
                    # 4. Process Mapping
                    if current_gene is not None and current_protein is not None:
                        chromosome_val = get_cleaned_str(row, 'chromosome')
                        chromosome_pos_val = safe_int_conversion(row.get('chromosome_position'))
                        map_strand_enum = current_gene.strand
                        if not chromosome_val or chromosome_pos_val is None:
                            logging.warning(f"L{line_num}: Missing/invalid critical mapping info for gene '{current_gene.gencode_transcription_id}'. Skipping.")
                        else:
                            current_mapping = Mapping(
                                chromosome=chromosome_val,
                                chromosome_position=chromosome_pos_val,
                                strand=map_strand_enum,
                                base_pair=get_cleaned_str(row, 'base_pair'),
                                codon=get_cleaned_str(row, 'codon'),
                                codon_base_pair_position=safe_int_conversion(row.get('codon_base_pair_position')),
                                amino_acid_residue=get_cleaned_str(row, 'amino_acid_residue'),
                                amino_acid_position=safe_int_conversion(row.get('amino_acid_position')),
                                cDNA_position=safe_int_conversion(row.get('cDNA_position')),
                                uniprot_residue=get_cleaned_str(row, 'uniprot_residue'),
                                uniprot_position=safe_int_conversion(row.get('uniprot_position')),
                                exon_number=safe_int_conversion(row.get('exon_number'))
                            )
                            current_mapping.gene = current_gene
                            current_mapping.protein = current_protein
                            sqlalchemy_session.add(current_mapping)
                    else:
                        logging.warning(f"L{line_num}: Missing Gene or Protein for mapping. Skipping mapping for gene '{gencode_tr_id}' and protein '{uniprot_ac}'.")

                    if current_interpro is not None and current_mapping is not None:
                        # Check if there is a meta-domain mapping
                        PFAM_consensus_pos_val = safe_int_conversion(row.get('PFAM_consensus_pos'))
                        if PFAM_consensus_pos_val is None:
                            current_meta_domain_mapping = MetaDomainMapping(consensus_position=PFAM_consensus_pos_val, ext_db_id=current_interpro.ext_db_id)
                            current_meta_domain_mapping.interpro_domain = current_interpro
                            current_meta_domain_mapping.mapping = current_mapping
                            sqlalchemy_session.add(current_meta_domain_mapping)
                        else:
                            logging.info(f"L{line_num}: Missing consensus position for mapping that does contain a protein domain at this position. Skipping meta-domain mapping for gene '{gencode_tr_id}' and protein '{uniprot_ac}'.")

                    current_batch_count += 1
                    if current_batch_count >= batch_size:
                        logging.info(f"Processed {total_processed_count} rows. Committing batch of {current_batch_count}.")
                        sqlalchemy_session.commit()
                        current_batch_count = 0
                except ValueError as ve:
                    logging.error(f"L{line_num}: Data validation error: {ve}. Row: {row}. Skipping.")
                    sqlalchemy_session.rollback()
                except IntegrityError as ie:
                    logging.error(f"L{line_num}: Database integrity error: {ie}. Row: {row}. Rolling back.")
                    sqlalchemy_session.rollback()
            
            if current_batch_count > 0:
                logging.info(f"Processed {total_processed_count} rows. Committing final batch of {current_batch_count}.")
                sqlalchemy_session.commit()
            logging.info(f"Successfully processed {total_processed_count} rows from {csv_filepath}")
    except FileNotFoundError:
        logging.error(f"CSV file not found: {csv_filepath}")
    except Exception as e:
        error_type = type(e).__name__
        error_traceback = traceback.format_exc()
        logging.error(f"L{line_num}: Unexpected error: {error_type}: {e}")
        logging.error(f"Traceback: {error_traceback}")
        logging.error(f"Row: {row}")
        logging.error(f"An unrecoverable error occurred: {e}")
        if sqlalchemy_session: sqlalchemy_session.rollback()
    # finally:
        # If running within app_context, session cleanup is typically handled.
        # If you manually obtained db.session, consider db.session.remove() or db.session.close()
        # but often Flask-SQLAlchemy manages this well.
        # sqlalchemy_session.close() # Generally not needed with Flask-SQLAlchemy's default session handling

def verify_csv_data_completeness(csv_filepath):
    """
    Verify that the CSV file contents are complete, so that:
    All mappings have a valid gene and protein, and the protein
     has all required columns.

     Finally return the results to be checked for database completeness check.
    :param csv_filepath: Path to the CSV file to verify.
    :return: A dictionary with counts of valid and invalid entries.
    """
    #@todo: Implement this function to check the CSV file for completeness.

def verify_database_completeness(db_session, csv_completeness_results):
    """
    Verify that the database contains all expected entries.
    This function should check that all proteins, genes, and mappings
    are correctly linked and have the required fields populated.

    :param db_session: The SQLAlchemy session to use for querying the database.
    :return: A dictionary with counts of valid and invalid entries.
    """


# --- Main Execution Block ---
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Batch load data from CSV into database.")
    parser.add_argument('--csv', type=str, required=True, help='Path to the CSV file to load.')
    
    args = parser.parse_args()
    
    csv_file = args.csv # Path to your CSV
 
    if not csv_file or not os.path.isfile(csv_file):
        logging.error("CSV file path is required. Use --csv to specify the path.")
        exit(1)

    # init the db with the flask app
    db.init_app(flask_app) 
   
    # This is crucial for standalone scripts using Flask-SQLAlchemy
    with flask_app.app_context():
        # Create tables if they don't exist
        # In a real app, this is often done once or via migrations (e.g., Flask-Migrate)
        db.create_all()

        # Verify CSV data completeness
        # logging.info("--- Verifying CSV Data Completeness ---") # @todo: Implement this function to check the CSV file for completeness.
        # completeness_results = verify_csv_data_completeness(csv_file)

        # Perform the batch load using db.session
        batch_load_data(csv_file, db.session)

        # Verification of database completeness
        # logging.info("--- Verifying Database Completeness ---") # @todo: Implement this function to check the database for completeness.
        # verify_database_completeness(db.session, completeness_results)

        # Cleanup session through explicit removal
        db.session.remove() # Flask-SQLAlchemy often handles this at context teardown, but for long-running scripts, explicit removal might be good.
