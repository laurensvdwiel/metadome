from metadome.database import db
import enum

class ProteinSource(enum.Enum):
    uniprot = 'uniprot'
    swissprot = 'swissprot'

class Protein(db.Model):
    """
    Table: proteins
    
    Fields
    id                         identifier
    uniprot_ac                 uniprot accession code
    uniprot_name               uniprot name
    uniprot_gn                 gene name from UniProt header GN=
    source                     'swissprot' or 'uniprot'
    evaluated_interpro_domains True if interpro domains have been annotated to this protein
    
    Relationships
    one to many                genes
    one to many                mappings
    one to many                interpro_domains
    """
    # Table configuration
    __tablename__ = 'proteins'
    
    # Fields
    id = db.Column(db.Integer, primary_key=True)
    uniprot_ac = db.Column(db.String(12), unique=True, nullable=False)
    uniprot_name = db.Column(db.String(20))
    uniprot_gn = db.Column(db.String(50))
    source = db.Column(db.Enum(ProteinSource), nullable=False)
    evaluated_interpro_domains = db.Column(db.Boolean)
    
    # Relationships
    genes = db.relationship('Gene', back_populates="protein")
    mappings = db.relationship('Mapping', back_populates="protein")
    interpro_domains = db.relationship("Interpro", back_populates="protein")
    
    def __init__(self, _uniprot_ac, _uniprot_name, _source, _uniprot_gn=None):
        if _source == 'swissprot':
            self.source = ProteinSource.swissprot
        elif _source == 'uniprot':
            self.source = ProteinSource.uniprot
        else:
            raise Exception('no source database defined for protein')
        self.uniprot_ac = _uniprot_ac
        self.uniprot_name = _uniprot_name
        self.uniprot_gn = _uniprot_gn
        self.evaluated_interpro_domains = False

    def __repr__(self):
        return "<Protein(uniprot_ac='%s', uniprot_name='%s', source='%s')>" % (
                            self.uniprot_ac, self.uniprot_name, self.source)