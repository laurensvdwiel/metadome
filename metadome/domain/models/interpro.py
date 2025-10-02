from metadome.database import db

class Interpro(db.Model):
    """
    Table: interpro_domains
    Contains the representation of Interpro domain occurrences
    in uniprot sequences
    
    Fields
    id                        identifier
    ext_db_id                 External domain database identifier code
    ext_db_version            External database version
    name                      Name of the Interpro domain
    domain_length             Length of the Interpro domain
    interpro_id               Interpro identifier
    uniprot_start             0 <= uniprot_start < uniprot_stop
    uniprot_stop              uniprot_stop <= uniprot_end
    protein_id                Foreign key
    
    Relationships
    many to one               protein
    """
    # Table configuration
    __tablename__ = 'interpro_domains'
    
    # Fields
    id = db.Column(db.Integer, primary_key=True)
    ext_db_id = db.Column(db.String, nullable=False)
    ext_db_version = db.Column(db.String, nullable=False)
    region_name = db.Column(db.String)
    interpro_id = db.Column(db.String(12))
    uniprot_start = db.Column(db.Integer, nullable=False)
    uniprot_stop = db.Column(db.Integer, nullable=False)
    protein_id = db.Column(db.Integer, db.ForeignKey('proteins.id'), nullable=False)
    
    # Relationships
    protein = db.relationship("Protein", back_populates="interpro_domains")
    meta_domain_associations = db.relationship("MetaDomainMappingAssociation",
                                               back_populates="interpro_domain")
    
    # Constraints
    __table_args__ = (db.UniqueConstraint('protein_id', 'ext_db_id', 'ext_db_version', 'uniprot_start', 'uniprot_stop', name='_unique_protein_region'),
                     )
    
    def __init__(self, _interpro_id, _ext_db_id, _ext_db_version, _region_name, _start_pos, _end_pos):
        self.interpro_id = _interpro_id
        self.ext_db_id = _ext_db_id
        self.ext_db_version = _ext_db_version
        self.region_name = _region_name
        self.uniprot_start = _start_pos
        self.uniprot_stop = _end_pos
    
    def __repr__(self):
        return "<Interpro(ext_db_id='%s', ext_db_version='%s', region_name='%s', interpro_id='%s', uniprot_start='%s', uniprot_stop='%s')>" % (
                            self.ext_db_id, self.ext_db_version, self.region_name, self.interpro_id, self.uniprot_start, self.uniprot_stop)