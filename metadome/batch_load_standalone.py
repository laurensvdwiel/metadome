import csv
import logging
import enum 
import argparse
import os
import sys
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.exc import IntegrityError
from dotenv import load_dotenv

# --- Environment Variables ---
# Load environment variables from .env file
load_dotenv()

# --- Flask App and SQLAlchemy Setup ---
# This would typically be in your app.py or a config file
flask_app = Flask(__name__)

# IMPORTANT: Configure your actual database URI
DB_HOST = os.getenv('DB_HOST', 'localhost')
DB_PORT = os.getenv('DB_PORT', 5432)
DB_USER = os.getenv('DB_USER')
DB_PWD = os.getenv('DB_PWD')
DB_NAME = os.getenv('DB_NAME', 'metadome')

flask_app.config['SQLALCHEMY_DATABASE_URI'] =  f"postgresql://{DB_USER}:{DB_PWD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
flask_app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # Optional: silence a warning

db = SQLAlchemy() # Initialize your db object globally or pass around

# --- Model Definitions (gene.py, protein.py, etc.) ---
# (These files would import 'db' from wherever you defined it above)

# gene.py:
class Strand(enum.Enum):
    plus = '+'
    minus = '-'

class Gene(db.Model): # inherits from the Flask-SQLAlchemy 'db.Model'
    __tablename__ = 'genes'
    id = db.Column(db.Integer, primary_key=True)
    strand = db.Column(db.Enum(Strand), nullable=False)
    gene_name = db.Column(db.String(50))
    gencode_transcription_id = db.Column(db.String(50), unique=True, nullable=False)
    gencode_translation_name = db.Column(db.String(50), nullable=False)
    gencode_gene_id = db.Column(db.String(50))
    havana_gene_id = db.Column(db.String(50))
    havana_translation_id = db.Column(db.String(50))
    sequence_length = db.Column(db.Integer)
    protein_id = db.Column(db.Integer, db.ForeignKey('proteins.id'))
    
    protein = db.relationship("Protein", back_populates="genes")
    mappings = db.relationship('Mapping', back_populates="gene")
    
    def __init__(self, _strand_str, _gene_name, _gencode_transcription_id, 
                 _gencode_translation_name, _gencode_gene_id, _havana_gene_id, 
                 _havana_translation_id, _sequence_length):
        # __init__ logic remains the same
        if _strand_str == '-':
            self.strand = Strand.minus
        elif _strand_str == '+':
            self.strand = Strand.plus
        else:
            raise ValueError(f"Invalid strand value: '{_strand_str}' for gene '{_gencode_transcription_id}'")
        self.gene_name = _gene_name if _gene_name and _gene_name.strip() else None
        self.gencode_transcription_id = _gencode_transcription_id
        self.gencode_translation_name = _gencode_translation_name
        self.gencode_gene_id = _gencode_gene_id if _gencode_gene_id and _gencode_gene_id.strip() else None
        self.havana_gene_id = None if _havana_gene_id == '-' or not (_havana_gene_id and _havana_gene_id.strip()) else _havana_gene_id
        self.havana_translation_id = None if _havana_translation_id == '-' or not (_havana_translation_id and _havana_translation_id.strip()) else _havana_translation_id
        self.sequence_length = _sequence_length
    
    def __repr__(self):
        return (f"<Gene(id={self.id}, strand='{self.strand}', gene_name='{self.gene_name}', "
                f"gencode_transcription_id='{self.gencode_transcription_id}')>")

# protein.py:
class ProteinSource(enum.Enum):
    uniprot = 'uniprot'
    swissprot = 'swissprot'

class Protein(db.Model):
    __tablename__ = 'proteins'
    id = db.Column(db.Integer, primary_key=True)
    uniprot_ac = db.Column(db.String(12), unique=True, nullable=False)
    uniprot_name = db.Column(db.String(20))
    source = db.Column(db.Enum(ProteinSource), nullable=False)
    evaluated_interpro_domains = db.Column(db.Boolean, nullable=False, default=False)
    
    genes = db.relationship('Gene', back_populates="protein")
    mappings = db.relationship('Mapping', back_populates="protein")
    interpro_domains = db.relationship("Interpro", back_populates="protein")
    
    def __init__(self, _uniprot_ac, _uniprot_name, _source_str, _evaluated_interpro_domains):
        try:
            self.source = ProteinSource(_source_str.lower())
        except ValueError:
             raise ValueError(f"Invalid source database: '{_source_str}' for protein '{_uniprot_ac}'. Must be 'uniprot' or 'swissprot'.")
        self.uniprot_ac = _uniprot_ac
        self.uniprot_name = _uniprot_name if _uniprot_name and _uniprot_name.strip() else None
        self.evaluated_interpro_domains = _evaluated_interpro_domains

    def __repr__(self):
        return f"<Protein(id={self.id}, uniprot_ac='{self.uniprot_ac}', source='{self.source}')>"

# mapping.py:
class Mapping(db.Model):
    __tablename__ = 'mappings'
    id = db.Column(db.Integer, primary_key=True)
    base_pair = db.Column(db.String(1))
    codon = db.Column(db.String(3))
    codon_base_pair_position = db.Column(db.Integer)
    strand = db.Column(db.Enum(Strand), nullable=False)
    amino_acid_residue = db.Column(db.String(1))
    amino_acid_position = db.Column(db.Integer)
    cDNA_position = db.Column(db.Integer)
    uniprot_residue = db.Column(db.String(1))
    uniprot_position = db.Column(db.Integer)
    chromosome = db.Column(db.String(5), nullable=False)
    chromosome_position = db.Column(db.Integer, nullable=False)    
    gene_id = db.Column(db.Integer, db.ForeignKey('genes.id'), nullable=False)
    protein_id = db.Column(db.Integer, db.ForeignKey('proteins.id'))
    
    gene = db.relationship("Gene", back_populates="mappings")
    protein = db.relationship("Protein", back_populates="mappings")
    
    def __init__(self, chromosome, chromosome_position, strand_str, gene_obj,
                 base_pair=None, codon=None, codon_base_pair_position=None, 
                 amino_acid_residue=None, amino_acid_position=None, 
                 cDNA_position=None, 
                 uniprot_residue=None, uniprot_position=None,
                 protein_obj=None):
        self.chromosome = chromosome
        self.chromosome_position = chromosome_position
        if strand_str == '+': self.strand = Strand.plus
        elif strand_str == '-': self.strand = Strand.minus
        else: raise ValueError(f"Invalid strand value for mapping: '{strand_str}' at {chromosome}:{chromosome_position}")
        self.gene = gene_obj
        self.base_pair = base_pair if base_pair and base_pair.strip() else None
        self.codon = codon if codon and codon.strip() else None
        self.codon_base_pair_position = codon_base_pair_position
        self.amino_acid_residue = amino_acid_residue if amino_acid_residue and amino_acid_residue.strip() else None
        self.amino_acid_position = amino_acid_position
        self.cDNA_position = cDNA_position
        self.uniprot_residue = uniprot_residue if uniprot_residue and uniprot_residue.strip() else None
        self.uniprot_position = uniprot_position
        if protein_obj: self.protein = protein_obj

    def __repr__(self):
        return (f"<Mapping(id={self.id}, chr='{self.chromosome}', chr_pos='{self.chromosome_position}', "
                f"gene_id='{self.gene_id}', protein_id='{self.protein_id}')>")

# interpro.py:
class Interpro(db.Model):
    __tablename__ = 'interpro_domains'
    id = db.Column(db.Integer, primary_key=True)
    ext_db_id = db.Column(db.String, nullable=False)
    region_name = db.Column(db.String)
    interpro_id = db.Column(db.String(12))
    uniprot_start = db.Column(db.Integer, nullable=False)
    uniprot_stop = db.Column(db.Integer, nullable=False)
    protein_id = db.Column(db.Integer, db.ForeignKey('proteins.id'), nullable=False)
    
    protein = db.relationship("Protein", back_populates="interpro_domains")
    __table_args__ = (db.UniqueConstraint('protein_id', 'ext_db_id', 'uniprot_start', 'uniprot_stop', name='_unique_protein_region'),)
    
    def __init__(self, _ext_db_id, _uniprot_start, _uniprot_stop, _interpro_id=None, _region_name=None):
        self.ext_db_id = _ext_db_id
        self.region_name = _region_name if _region_name and _region_name.strip() else None
        self.interpro_id = _interpro_id if _interpro_id and _interpro_id.strip() else None
        self.uniprot_start = _uniprot_start
        self.uniprot_stop = _uniprot_stop
    
    def __repr__(self):
        return (f"<Interpro(id={self.id}, ext_db_id='{self.ext_db_id}', protein_id='{self.protein_id}', "
                f"start='{self.uniprot_start}', stop='{self.uniprot_stop}')>")

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
        with open(csv_filepath, mode='r', encoding='utf-8') as infile:
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
                            if not current_protein:
                                uniprot_name_val = get_cleaned_str(row, 'uniprot_name')
                                source_val_str = get_cleaned_str(row, 'source')
                                eval_interpro_val = safe_bool_conversion(row.get('evaluated_interpro_domains'), False)
                                if not source_val_str:
                                    logging.warning(f"L{line_num}: Missing source for Protein '{uniprot_ac}'. Skipping Protein.")
                                else:
                                    current_protein = Protein(_uniprot_ac=uniprot_ac, 
                                                              _uniprot_name=uniprot_name_val, 
                                                              _source_str=source_val_str,
                                                              _evaluated_interpro_domains=eval_interpro_val)
                                    sqlalchemy_session.add(current_protein)
                            protein_cache[uniprot_ac] = current_protein
                    
                    # 2. Process Gene
                    gencode_tr_id = get_cleaned_str(row, 'gencode_transcription_id')
                    if gencode_tr_id:
                        if gencode_tr_id in gene_cache:
                            current_gene = gene_cache[gencode_tr_id]
                        else:
                            current_gene = sqlalchemy_session.query(Gene).filter_by(gencode_transcription_id=gencode_tr_id).one_or_none()
                            if not current_gene:
                                strand_str = get_cleaned_str(row, 'strand')
                                gene_name_val = get_cleaned_str(row, 'gene_name')
                                gencode_transl_name_val = get_cleaned_str(row, 'gencode_translation_name')
                                if not strand_str or strand_str not in ['+', '-']:
                                    logging.warning(f"L{line_num}: Invalid strand for Gene '{gencode_tr_id}'. Skipping Gene.")
                                elif not gencode_transl_name_val:
                                    logging.warning(f"L{line_num}: Missing gencode_translation_name for Gene '{gencode_tr_id}'. Skipping Gene.")
                                else:
                                    current_gene = Gene(
                                        _strand_str=strand_str, _gene_name=gene_name_val,
                                        _gencode_transcription_id=gencode_tr_id, _gencode_translation_name=gencode_transl_name_val,
                                        _gencode_gene_id=get_cleaned_str(row, 'gencode_gene_id'),
                                        _havana_gene_id=get_cleaned_str(row, 'havana_gene_id'),
                                        _havana_translation_id=get_cleaned_str(row, 'havana_translation_id'),
                                        _sequence_length=safe_int_conversion(row.get('sequence_length')))
                                    if current_protein: current_gene.protein = current_protein
                                    sqlalchemy_session.add(current_gene)
                            gene_cache[gencode_tr_id] = current_gene
                    
                    # 3. Process Interpro
                    ext_db_id_val = get_cleaned_str(row, 'ext_db_id')
                    uniprot_start_val = safe_int_conversion(row.get('uniprot_start'))
                    uniprot_stop_val = safe_int_conversion(row.get('uniprot_stop'))
                    if current_protein and ext_db_id_val and uniprot_start_val is not None and uniprot_stop_val is not None:
                        interpro_cache_key = (current_protein.uniprot_ac, ext_db_id_val, uniprot_start_val, uniprot_stop_val)
                        if interpro_cache_key not in interpro_cache:
                            existing_interpro = None
                            if current_protein.id:
                                existing_interpro = sqlalchemy_session.query(Interpro).filter_by(
                                    protein_id=current_protein.id, ext_db_id=ext_db_id_val,
                                    uniprot_start=uniprot_start_val, uniprot_stop=uniprot_stop_val).one_or_none()
                            if not existing_interpro:
                                interpro_obj = Interpro(
                                    _ext_db_id=ext_db_id_val, _uniprot_start=uniprot_start_val, _uniprot_stop=uniprot_stop_val,
                                    _interpro_id=get_cleaned_str(row, 'interpro_id'), _region_name=get_cleaned_str(row, 'name'))
                                interpro_obj.protein = current_protein
                                sqlalchemy_session.add(interpro_obj)
                                interpro_cache[interpro_cache_key] = interpro_obj
                            else:
                                interpro_cache[interpro_cache_key] = existing_interpro
                    
                    # 4. Process Mapping
                    if current_gene:
                        chromosome_val = get_cleaned_str(row, 'chromosome')
                        chromosome_pos_val = safe_int_conversion(row.get('chromosome_position'))
                        map_strand_str = get_cleaned_str(row, 'strand')
                        if not chromosome_val or chromosome_pos_val is None or not map_strand_str or map_strand_str not in ['+', '-']:
                            logging.warning(f"L{line_num}: Missing/invalid critical mapping info for gene '{current_gene.gencode_transcription_id}'. Skipping.")
                        else:
                            mapping_obj = Mapping(
                                chromosome=chromosome_val, chromosome_position=chromosome_pos_val, strand_str=map_strand_str,
                                gene_obj=current_gene, base_pair=get_cleaned_str(row, 'base_pair'),
                                codon=get_cleaned_str(row, 'codon'), codon_base_pair_position=safe_int_conversion(row.get('codon_base_pair_position')),
                                amino_acid_residue=get_cleaned_str(row, 'amino_acid_residue'), amino_acid_position=safe_int_conversion(row.get('amino_acid_position')),
                                cDNA_position=safe_int_conversion(row.get('cDNA_position')),
                                uniprot_residue=get_cleaned_str(row, 'uniprot_residue'), uniprot_position=safe_int_conversion(row.get('uniprot_position')),
                                protein_obj=current_protein)
                            sqlalchemy_session.add(mapping_obj)
                    
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
                except Exception as e:
                    logging.error(f"L{line_num}: Unexpected error: {e}. Row: {row}. Skipping.")
                    sqlalchemy_session.rollback()
            
            if current_batch_count > 0:
                logging.info(f"Processed {total_processed_count} rows. Committing final batch of {current_batch_count}.")
                sqlalchemy_session.commit()
            logging.info(f"Successfully processed {total_processed_count} rows from {csv_filepath}")
    except FileNotFoundError:
        logging.error(f"CSV file not found: {csv_filepath}")
    except Exception as e:
        logging.error(f"An unrecoverable error occurred: {e}")
        if sqlalchemy_session: sqlalchemy_session.rollback()
    # finally:
        # If running within app_context, session cleanup is typically handled.
        # If you manually obtained db.session, consider db.session.remove() or db.session.close()
        # but often Flask-SQLAlchemy manages this well.
        # sqlalchemy_session.close() # Generally not needed with Flask-SQLAlchemy's default session handling


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

        # Perform the batch load using db.session
        batch_load_data(csv_file, db.session, batch_size=2)

        # Verification (optional)
        logging.info("--- Verification ---")
        for protein in db.session.query(Protein).all():
            logging.info(f"Loaded Protein: {protein} (Eval Interpro: {protein.evaluated_interpro_domains}) with {len(protein.genes)} genes and {len(protein.interpro_domains)} domains.")
        for gene in db.session.query(Gene).all():
            logging.info(f"Loaded Gene: {gene} (Protein ID: {gene.protein_id}) with {len(gene.mappings)} mappings.")
        for mapping_count, mapping in enumerate(db.session.query(Mapping).all()):
            logging.info(f"Loaded Mapping {mapping_count+1}: {mapping}")
        
        # db.session.remove() # Or db.session.close() - Flask-SQLAlchemy often handles this at context teardown
                          # For long-running scripts, explicit removal might be good.
