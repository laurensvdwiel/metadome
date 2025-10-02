from metadome.database import db

class MetaDomainMappingAssociation(db.Model):
    """
    Table: meta_domain_mapping_association
    Association model to represent the many-to-many relationship between
    mappings, interpro_domains, and meta_domain_mappings

    Fields
    id                      identifier
    mapping_id              Foreign key to mappings table
    interpro_id             Foreign key to interpro_domains table
    meta_domain_mapping_id  Foreign key to meta_domain_mappings table

    Relationships
    many to one             mapping
    many to one             meta_domain_mapping
    many to one             interpro_domain
    """

    # Table configuration
    __tablename__ = 'meta_domain_mapping_association'

    # Fields
    id = db.Column(db.Integer, primary_key=True)
    mapping_id = db.Column(db.Integer, db.ForeignKey('mappings.id'), nullable=False)
    interpro_id = db.Column(db.Integer, db.ForeignKey('interpro_domains.id'), nullable=False)
    meta_domain_mapping_id = db.Column(db.Integer, db.ForeignKey('meta_domain_mappings.id'), nullable=False)

    # Relationships
    mapping = db.relationship("Mapping", back_populates="meta_domain_associations")
    meta_domain_mapping = db.relationship("MetaDomainMapping", back_populates="meta_domain_mapping_associations")
    interpro_domain = db.relationship("Interpro", back_populates="meta_domain_associations")

    # Constraints
    __table_args__ = (
        db.UniqueConstraint('mapping_id', 'meta_domain_mapping_id', 'interpro_id', name='_unique_meta_domain_mapping_association'),
    )