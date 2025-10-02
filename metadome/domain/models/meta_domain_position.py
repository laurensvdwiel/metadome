from metadome.database import db

class MetaDomainPosition(db.Model):
    """
    Table: meta_domain_positions
    Representation of a single meta-domain position in a protein domain
    
    Fields
    id                        identifier
    consensus_position        Integer indicating the consensus position in the protein domain
    ext_db_id                 External domain database identifier code
    
    Relationships
    one to many               meta_domain_mapping_association
    """
    # Table configuration
    __tablename__ = 'meta_domain_positions'
    
    # Fields
    id = db.Column(db.Integer, primary_key=True)
    consensus_position = db.Column(db.Integer, nullable=False)
    ext_db_id = db.Column(db.String, nullable=False)
    
    # Relationships
    meta_domain_mapping_associations = db.relationship("MetaDomainMappingAssociation",
                                              back_populates="meta_domain_position")

    # Constraints
    __table_args__ = (db.UniqueConstraint('ext_db_id', 'consensus_position',
                                          name='_unique_meta_domain_position'),
                      )

    def __repr__(self):
        return "<MetaDomainMapping(ext_db_id='%s', consensus_position='%s')>" % (
                        self.ext_db_id, self.consensus_position)