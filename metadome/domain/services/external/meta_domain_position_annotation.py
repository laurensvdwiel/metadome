import logging

from metadome.domain.repositories import MetaDomainRepository

_log = logging.getLogger(__name__)


def annotate_variants_with_metadomain(variants):
    """
    Annotate normalized variant dictionaries.

    Expected input:
        [
            {
                "chr": "chr1",
                "pos": 123456,
                "ref": "A",
                "gene_id": "ENSG00000000001",
                "genome_build": "GRCh38.p14",
            },
            ...
        ]

    Returns:
        [
            {
                "chr": ...,
                "pos": ...,
                "ref": ...,
                "gene_id": ...,
                "genome_build": ...,
                "MetaDomainPositions": ...,
                "MetaDomainStatus": ...,
                "RefMatchStatus": ...,
            },
            ...
        ]
    """
    normalized_variants = []

    for variant in variants:
        chr_val = variant.get("chr")
        pos_val = variant.get("pos")
        ref_val = variant.get("ref")
        gene_id = variant.get("gene_id")
        genome_build = variant.get("genome_build")

        if not all([chr_val, pos_val, gene_id, genome_build]):
            normalized_variants.append({
                **variant,
                "MetaDomainPositions": "",
                "MetaDomainStatus": "invalid_input",
                "RefMatchStatus": "not_checked",
            })
            continue

        normalized_variants.append({
            "chr": chr_val,
            "pos": int(pos_val),
            "ref": ref_val,
            "gene_id": gene_id,
            "genome_build": genome_build,
        })

    valid_variants = [
        variant for variant in normalized_variants
        if variant.get("MetaDomainStatus") is None
        or "MetaDomainStatus" not in variant
    ]

    grouped_hits = MetaDomainRepository.get_meta_domain_annotation_for_variants(valid_variants)

    results = []
    for variant in normalized_variants:
        if variant.get("MetaDomainStatus") == "invalid_input":
            results.append(variant)
            continue

        normalized_gene_id = variant["gene_id"].split(".", 1)[0]
        variant_key = (
            variant["chr"],
            variant["pos"],
            variant.get("ref"),
            normalized_gene_id,
            variant["genome_build"],
        )

        annotation = grouped_hits.get(variant_key, {
            "MetaDomainPositions": "",
            "MetaDomainStatus": "no_mapping",
            "RefMatchStatus": "not_checked",
        })

        results.append({
            **variant,
            "MetaDomainPositions": annotation["MetaDomainPositions"],
            "MetaDomainStatus": annotation["MetaDomainStatus"],
            "RefMatchStatus": annotation["RefMatchStatus"],
        })

    return results