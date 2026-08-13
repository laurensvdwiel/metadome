table metaDomeTolerance
"MetaDome: missense over synonymous ratio computed over a sliding window of 21 codons centred on each position, corrected for codon composition. One feature per genomic codon; where several transcripts report different values at the same codon, the median is used."
    (
    string  chrom;        "Reference sequence chromosome"
    uint    chromStart;   "Start position of the codon"
    uint    chromEnd;     "End position of the codon"
    string  name;         "Unused"
    uint    score;        "Unused"
    char[1] strand;       "Strand of the reading frame"
    lstring sw_dn_ds;     "Missense over synonymous ratio; lower values indicate greater intolerance to missense variation"
    lstring sw_coverage;  "Fraction of the sliding window available; below 1 near the start and end of a coding region"
    )
