from metadome.database import db

class MetaDomainMapping(db.Model):
    """
    Table: meta_domain_mapping
    Representation of a single meta-domain mapping between a mapping and gencode gene
    translation and a uniprot protein sequence
    
    Fields
    id                        identifier
    consensus_position        Integer indicating the consensus position in the protein domain
    ext_db_id                 External domain database identifier code
    interpro_id               Foreign key
    mapping_id                Foreign key
    
    Relationships
    many to many              interpro
    many to many              mapping
    """
    # Table configuration
    __tablename__ = 'meta_domain_mappings'
    
    # Fields
    id = db.Column(db.Integer, primary_key=True)
    consensus_position = db.Column(db.Integer, nullable=False)
    ext_db_id = db.Column(db.String, nullable=False)
    mapping_id = db.Column(db.Integer, db.ForeignKey('mappings.id'), nullable=False)
    interpro_id = db.Column(db.Integer, db.ForeignKey('interpro_domains.id'), nullable=False)
    
    # Relationships
    interpro_domain = db.relationship("Interpro", back_populates="meta_domain_mappings")
    mapping = db.relationship("Mapping", back_populates="meta_domain_mapping")

    # Constraints
    __table_args__ = (db.UniqueConstraint('ext_db_id', 'consensus_position', 'mapping_id',
                                          name='_unique_meta_domain_position'),
                      )

    def __repr__(self):
        return "<MetaDomainMapping(ext_db_id='%s', consensus_position='%s')>" % (
                        self.ext_db_id, self.consensus_position)