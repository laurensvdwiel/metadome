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

function isValidGenomeBuild(genomeBuild) {
    const gb = String(genomeBuild || "").toUpperCase();
    return gb.startsWith("GRCH37") || gb.startsWith("GRCH38");
}

function isNonEmptyString(value) {
    return typeof value === "string" && value.trim().length > 0;
}

function isPositiveIntegerLike(value) {
    if (value === null || typeof value === "undefined") return false;
    const n = Number(value);
    return Number.isInteger(n) && n > 0;
}

function normalizeChromosome(chromosome) {
    const raw = String(chromosome || "").trim();
    if (raw.length === 0) return null;
    const noChr = raw.replace(/^chr/i, "");
    if (/^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/i.test(noChr)) {
        return noChr.toUpperCase() === "MT" ? "M" : noChr.toUpperCase();
    }
    return null;
}

window.METADOME_LINKS = {
    ensemblTranscript(genomeBuild, transcriptId) {
        if (!isValidGenomeBuild(genomeBuild) || !isNonEmptyString(String(transcriptId || ""))) {
            return null;
        }
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__"
            : "https://ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__";
        return fillTemplate(template, { transcript_id: String(transcriptId).trim() });
    },

    ensemblGenomicPosition(genomeBuild, chromosome, position) {
        const chr = normalizeChromosome(chromosome);
        if (!isValidGenomeBuild(genomeBuild) || !chr || !isPositiveIntegerLike(position)) {
            return null;
        }
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;r=__CHR__:__POS__-__POS__"
            : "https://ensembl.org/Homo_sapiens/Location/View?r=__CHR__:__POS__-__POS__";
        return fillTemplate(template, { chromosome: chr, position: Number(position) });
    },

    gnomadVariant(genomeBuild, chromosome, position, ref, alt) {
        const chr = normalizeChromosome(chromosome);
        const refOk = isNonEmptyString(String(ref || ""));
        const altOk = isNonEmptyString(String(alt || ""));
        if (!isValidGenomeBuild(genomeBuild) || !chr || !isPositiveIntegerLike(position) || !refOk || !altOk) {
            return null;
        }
        const template = isGrch37(genomeBuild)
            ? "https://gnomad.broadinstitute.org/variant/__CHR__-__POS__-__REF__-__ALT__?dataset=gnomad_r2_1"
            : "https://gnomad.broadinstitute.org/variant/__CHR__-__POS__-__REF__-__ALT__?dataset=gnomad_r4";
        return fillTemplate(template, { chromosome: chr, position: Number(position), ref: ref, alt: alt });
    },

    clinvarVariant(clinvarId) {
        if (!isPositiveIntegerLike(clinvarId) && !isNonEmptyString(String(clinvarId || ""))) {
            return null;
        }
        return fillTemplate("https://www.ncbi.nlm.nih.gov/clinvar/variation/__CLINVAR_ID__/", {
            clinvar_id: String(clinvarId).trim()
        });
    },

    refseq(refseqId) {
        if (!isNonEmptyString(String(refseqId || ""))) {
            return null;
        }
        return fillTemplate("https://www.ncbi.nlm.nih.gov/nuccore/__REFSEQ__", {
            refseq: String(refseqId).trim()
        });
    },

    pfam(pfamId) {
        if (!isNonEmptyString(String(pfamId || ""))) {
            return null;
        }
        return fillTemplate("https://www.ebi.ac.uk/interpro/entry/pfam/__PFAM__/", {
            pfam: String(pfamId).trim()
        });
    },

    uniprot(uniprotAc) {
        if (!isNonEmptyString(String(uniprotAc || ""))) {
            return null;
        }
        return fillTemplate("https://www.uniprot.org/uniprotkb/__UNIPROT__", {
            uniprot: String(uniprotAc).trim()
        });
    }
};