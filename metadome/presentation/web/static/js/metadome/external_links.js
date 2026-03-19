function fillTemplate(template, values) {
    return template
        .replaceAll("__GB__", encodeURIComponent(values.genome_build || ""))
        .replaceAll("__CHR__", encodeURIComponent(values.chromosome || ""))
        .replaceAll("__POS__", encodeURIComponent(values.position || ""))
        .replaceAll("__REF__", encodeURIComponent(values.ref || ""))
        .replaceAll("__ALT__", encodeURIComponent(values.alt || ""))
        .replaceAll("__CLINVAR_ID__", encodeURIComponent(values.clinvar_id || ""))
        .replaceAll("__REFSEQ__", encodeURIComponent(values.refseq || ""))
        .replaceAll("__PFAM__", encodeURIComponent(values.pfam || ""))
        .replaceAll("__UNIPROT__", encodeURIComponent(values.uniprot || ""))
        .replaceAll("__TRANSCRIPT__", encodeURIComponent(values.transcript_id || ""))
        .replaceAll("__GENE__", encodeURIComponent(values.gene_name || ""))
        .replaceAll("__AA_POSITION__", encodeURIComponent(values.aa_position || ""));
}

function isGrch37(genomeBuild) {
    return String(genomeBuild || "").toUpperCase().startsWith("GRCH37");
}

window.METADOME_LINKS = {
    ensemblTranscript(genomeBuild, transcriptId) {
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__"
            : "https://ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__";
        return fillTemplate(template, { transcript_id: transcriptId });
    },

    ensemblGenomicPosition(genomeBuild, chromosome, position) {
        const chr = String(chromosome || "").replace(/^chr/i, "");
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;r=__CHR__:__POS__-__POS__"
            : "https://ensembl.org/Homo_sapiens/Location/View?r=__CHR__:__POS__-__POS__";
        return fillTemplate(template, { chromosome: chr, position });
    },

    gnomadVariant(genomeBuild, chromosome, position, ref, alt) {
        const chr = String(chromosome || "").replace(/^chr/i, "");
        const template = isGrch37(genomeBuild)
            ? "https://gnomad.broadinstitute.org/variant/__CHR__-__POS__-__REF__-__ALT__?dataset=gnomad_r2_1"
            : "https://gnomad.broadinstitute.org/variant/__CHR__-__POS__-__REF__-__ALT__?dataset=gnomad_r4";
        return fillTemplate(template, { chromosome: chr, position, ref, alt });
    },

    clinvarVariant(clinvarId) {
        return fillTemplate("https://www.ncbi.nlm.nih.gov/clinvar/variation/__CLINVAR_ID__/", {
            clinvar_id: clinvarId
        });
    },

    refseq(refseqId) {
        return fillTemplate("https://www.ncbi.nlm.nih.gov/nuccore/__REFSEQ__", {
            refseq: refseqId
        });
    },

    pfam(pfamId) {
        return fillTemplate("https://www.ebi.ac.uk/interpro/entry/pfam/__PFAM__/", {
            pfam: pfamId
        });
    },

    uniprot(uniprotAc) {
        return fillTemplate("https://www.uniprot.org/uniprotkb/__UNIPROT__", {
            uniprot: uniprotAc
        });
    }
};