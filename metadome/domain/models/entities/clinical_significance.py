import enum

class ClinicalSignificance(enum.Enum):
    benign = 'Benign'
    likely_benign = 'Likely_benign'
    uncertain_significance = 'Uncertain_significance'
    likely_pathogenic = 'Likely_pathogenic'
    pathogenic = 'Pathogenic'


CLINVAR_CLINSIG_TO_SIGNIFICANCE = {
    'Benign': ClinicalSignificance.benign,
    'Likely_benign': ClinicalSignificance.likely_benign,
    'Benign/Likely_benign': ClinicalSignificance.likely_benign,
    'Uncertain_significance': ClinicalSignificance.uncertain_significance,
    'Likely_pathogenic': ClinicalSignificance.likely_pathogenic,
    'Pathogenic': ClinicalSignificance.pathogenic,
    'Pathogenic/Likely_pathogenic': ClinicalSignificance.likely_pathogenic,
}