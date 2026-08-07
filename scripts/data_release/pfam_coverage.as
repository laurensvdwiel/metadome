table metaDomainPfamCoverage
"MetaDome: Pfam domain coverage of the coding genome. Each feature is a codon that aligns to a position in a Pfam domain consensus. Codons sharing a Pfam ID and consensus position are evolutionarily equivalent, so any variant set can be lifted between them."
    (
    string  chrom;              "Reference sequence chromosome"
    uint    chromStart;         "Start position of the codon"
    uint    chromEnd;           "End position of the codon"
    string  name;               "UniProt accession/position : Pfam ID : domain consensus position"
    uint    score;              "Unused; always 0"
    char[1] strand;             "Strand of the reading frame"
    string  uniprot_ac;         "UniProtKB/Swiss-Prot accession, isoform-suffixed where applicable"
    uint    uniprot_pos;        "Residue position, 1-based from the initiator methionine"
    string  pfam_id;            "Pfam domain accession"
    uint    consensus_pos;      "Position in the Pfam domain consensus; equal positions across genes are evolutionarily equivalent"
    lstring metadome_url;       "Link to this genomic position in MetaDome"
    )
