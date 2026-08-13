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
    if (/^(?:[1-9]|1[0-9]|2[0-2]|X|Y)$/i.test(noChr)) {
        return "chr" + noChr.toUpperCase();
    }
    return null;
}

const STATIC_EXTERNAL_LINKS = {
    MetaDome: {
        app: "https://www.metadome.app",
        legacy: "https://stuart.radboudumc.nl/metadome/",
        github: "https://github.com/laurensvdwiel/metadome",
        githubReleases: "https://github.com/laurensvdwiel/metadome/releases",
        githubIssues: "https://github.com/laurensvdwiel/metadome/issues",
        zenodoV1: "https://zenodo.org/records/6625251",
        zenodoV2: "https://zenodo.org/records/19376150"
    },
    "MetaDome.papers": {
        humu2017: "https://doi.org/10.1002/humu.23313",
        humu2019: "https://doi.org/10.1002/humu.23798",
        humu2019Cover: "https://doi.org/10.1002/humu.23892",
        ajhg2023: "https://doi.org/10.1016/j.ajhg.2022.12.001",
        nature2020: "https://doi.org/10.1038/s41586-020-2832-5"
    },
    "External.resources": {
        gnomad: "http://gnomad.broadinstitute.org/",
        clinvar: "https://www.ncbi.nlm.nih.gov/clinvar/",
        pfam: "https://www.ebi.ac.uk/interpro/entry/pfam/",
        interpro: "https://www.ebi.ac.uk/interpro/",
        gencode: "https://www.gencodegenes.org/",
        uniprot: "http://www.uniprot.org/",
        docker: "https://www.docker.com/",
        bulma: "https://bulma.io/",
        d3: "https://d3js.org/",
        flask: "http://flask.pocoo.org/",
        hmmer: "http://hmmer.org/",
        postgresql: "https://www.postgresql.org/",
        redis: "https://redis.io/",
        celery: "https://docs.celeryq.dev/",
        rabbitmq: "https://www.rabbitmq.com/",
        chrome: "https://www.google.com/chrome/",
        firefox: "https://www.mozilla.org/en-US/firefox/new/",
        jalview: "http://www.jalview.org/",
        ddd: "https://en.wikipedia.org/wiki/Domain-driven_design",
        ucscConstraintTrack: "https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1735422030_sAYmaQ0c58lMHr1JWvkyZ4gtGBGI&db=hg38&c=chr1&g=constraintSuper",
        mobidetails: "https://mobidetails.iurc.montp.inserm.fr/MD/about",
        opensourceMit: "http://opensource.org/licenses/mit-license.php",
        ccByNcSa40: "http://creativecommons.org/licenses/by-nc-sa/4.0/"
    },
    "External.privacy": {
        googlePrivacy: "https://policies.google.com/privacy",
        googleTerms: "https://policies.google.com/terms",
        gaIpPrivacy: "https://support.google.com/analytics/answer/12017362"
    },
    "External.affiliates": {
        stanford: "https://www.stanford.edu/",
        stanfordMedicine: "https://med.stanford.edu/",
        radboudumcResearch: "https://www.radboudumc.nl/en/research",
        radboudHumanGenetics: "https://www.radboudumc.nl/en/research/research-groups/genome-bioinformatics",
        cmbi: "https://en.wikipedia.org/wiki/Center_for_Molecular_and_Biomolecular_Informatics"
    }
};

function getStaticExternalLink(group, key) {
    if (!isNonEmptyString(group) || !isNonEmptyString(key)) {
        return null;
    }

    const groupLinks = STATIC_EXTERNAL_LINKS[group] || {};
    return groupLinks[key] || null;
}

function applyStaticExternalLinks(root) {
    const scope = root || document;

    scope.querySelectorAll("[data-external-link-group][data-external-link-key]").forEach(function (element) {
        const href = getStaticExternalLink(
            element.dataset.externalLinkGroup,
            element.dataset.externalLinkKey
        );

        if (href) {
            element.href = href;
        }
    });
}

window.METADOME_LINKS = {
    ensemblTranscript(genomeBuild, transcriptId) {
        if (!isValidGenomeBuild(genomeBuild) || !isNonEmptyString(String(transcriptId || ""))) {
            return null;
        }
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__"
            : "https://ensembl.org/Homo_sapiens/Transcript/Summary?t=__TRANSCRIPT__";
        return fillTemplate(template, {transcript_id: String(transcriptId).trim()});
    },

    ensemblGenomicPosition(genomeBuild, chromosome, position) {
        const chr = normalizeChromosome(chromosome);
        if (!isValidGenomeBuild(genomeBuild) || !chr || !isPositiveIntegerLike(position)) {
            return null;
        }
        const template = isGrch37(genomeBuild)
            ? "https://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;r=__CHR__:__POS__-__POS__"
            : "https://ensembl.org/Homo_sapiens/Location/View?r=__CHR__:__POS__-__POS__";
        return fillTemplate(template, {chromosome: chr, position: Number(position)});
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
        return fillTemplate(template, {chromosome: chr, position: Number(position), ref: ref, alt: alt});
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
    },

    get(group, key) {
        return getStaticExternalLink(group, key);
    },

    apply(root) {
        applyStaticExternalLinks(root);
    }
};

document.addEventListener("DOMContentLoaded", function () {
    applyStaticExternalLinks(document);
});