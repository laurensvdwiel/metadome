from metadome.database import db
from metadome.domain.models.gene import Strand

class Mapping(db.Model):
    """
    Table: mapping
    Representation of a single positional mapping between a gencode gene 
    translation and a uniprot protein sequence
    
    Fields
    id                        identifier
    base_pair                    'A', 'T', 'C', 'G'
    codon                     all triplet combinations of the alleles 
    codon_base_pair_position     0, 1, 2, or None
    amino_acid_residue        one of the 20 amino acids in the gene translation,
                              * or None
    amino_acid_position       the position in the amino acid sequence of the 
                              gene translation
    cDNA_position             the position in the cDNA
    uniprot_residue           one of the 20 amino acids in the uniprot sequence,
                              * or None
    uniprot_position          the position in the protein
    chromosome                the chromosome this mapping points to
    chromosome_position       the chromosomal position of this mapping
    exon_number               the numerical position of an exon within a gene transcript
    gene_id                   Foreign key
    protein_id                Foreign key
    
    Relationships
    many to one               gene
    many to one               protein
    """
    # Table configuration
    __tablename__ = 'mappings'
    
    # Fields
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
    exon_number = db.Column(db.Integer, nullable=False)
    gene_id = db.Column(db.Integer, db.ForeignKey('genes.id'), nullable=False)
    protein_id = db.Column(db.Integer, db.ForeignKey('proteins.id'))
    
    # Relationships
    gene = db.relationship("Gene", back_populates="mappings")
    protein = db.relationship("Protein", back_populates="mappings")
    meta_domain_associations = db.relationship("MetaDomainMappingAssociation",
                                               back_populates="mapping")
    
    def __repr__(self):
        return "<Mapping(chr='%s', chr_pos='%s' base_pair='%s', codon='%s', codon_base_pair_position='%s', amino_acid_residue='%s', cDNA_position='%s', uniprot_position='%s')>" % (
                            self.chromosome, self.chromosome_position, self.base_pair, self.codon, self.codon_base_pair_position, self.amino_acid_residue, self.cDNA_position, self.uniprot_position)
    