from metadome.database import db

class MetaDomainMapping(db.Model):
    """
    Table: meta_domain_mapping
    Association table to represent the many-to-one relationship between
    mappings, interpro_domains, and meta_domain_positions

    Fields
    id                      identifier
    mapping_id              Foreign key to mappings table
    interpro_id             Foreign key to interpro_domains table
    meta_domain_position_id Foreign key to meta_domain_positions table

    Relationships
    many to one             mapping
    many to one             meta_domain_position
    many to one             interpro_domain
    """

    # Table configuration
    __tablename__ = 'meta_domain_mapping'

    # Fields
    id = db.Column(db.Integer, primary_key=True)
    mapping_id = db.Column(db.Integer, db.ForeignKey('mappings.id'), nullable=False)
    interpro_id = db.Column(db.Integer, db.ForeignKey('interpro_domains.id'), nullable=False)
    meta_domain_position_id = db.Column(db.Integer, db.ForeignKey('meta_domain_positions.id'), nullable=False)

    # Relationships
    mapping = db.relationship("Mapping", back_populates="meta_domain_mappings")
    meta_domain_position = db.relationship("MetaDomainPosition", back_populates="meta_domain_mappings")
    interpro_domain = db.relationship("Interpro", back_populates="meta_domain_mappings")

    # Constraints
    __table_args__ = (
        db.UniqueConstraint('mapping_id', 'meta_domain_position_id', 'interpro_id', name='_unique_meta_domain_mapping_association'),
        db.Index('idx_mdm_position_id', 'meta_domain_position_id'),
    )