document.addEventListener("DOMContentLoaded", function() {
    if (!window.METADOME_LINKS) {
        return;
    }

    function buildExternalAnchorOrFallback(href, label) {
        if (!href || typeof href !== "string" || href.trim().length === 0) {
            return '<span class="has-text-grey" title="External link unavailable">' + label + '</span>';
        }
        return '<a href="' + href + '" target="_blank" rel="noopener">' + label + '</a>';
    }

    function normalizeStrandSymbol(strandValue) {
        const raw = String(strandValue || "").trim();
        if (raw === "+" || raw === "-") {
            return raw;
        }
        if (raw.toLowerCase() === "plus") {
            return "+";
        }
        if (raw.toLowerCase() === "minus") {
            return "-";
        }
        return "";
    }

    function formatCodonWithHighlight(codon, codonPosition, strandSymbol) {
        const clean = String(codon || "").trim();
        const pos = parseInt(codonPosition, 10);
        const strand = String(strandSymbol || "").trim();

        if (clean.length !== 3 || ![1, 2, 3].includes(pos)) {
            return strand ? "(" + clean + " " + strand + ")" : "(" + clean + ")";
        }

        const chars = clean.split("");
        const index = pos - 1;
        chars[index] = "<strong><u>" + chars[index] + "</u></strong>";
        return strand ? "(" + chars.join("") + " " + strand + ")" : "(" + chars.join("") + ")";
    }

    const queryElement = document.getElementById("position-query");
    if (queryElement) {
        const genomeBuild = queryElement.dataset.genomeBuild || window.METADOME_CONFIG?.genomeBuild || "";
        const chr = queryElement.dataset.chr || "";
        const pos = queryElement.dataset.pos || "";
        const label = chr + ":" + pos + " (" + genomeBuild + ")";
        const href = window.METADOME_LINKS.ensemblGenomicPosition(genomeBuild, chr, pos);
        const container = queryElement.querySelector(".query-link-container");

        if (container) {
            container.innerHTML = buildExternalAnchorOrFallback(href, label);
        }
    }

    const table = document.getElementById("position-results-table");
    if (!table) {
        return;
    }

    const rows = table.querySelectorAll("tbody tr[data-gene]");
    rows.forEach(function(row) {
        const genomeBuild = row.dataset.build || window.METADOME_CONFIG?.genomeBuild || "GRCh38";
        const gencodeTranscript = row.dataset.gencodeTranscript || "";
        const refseqTranscript = row.dataset.refseqTranscript || "";
        const uniprotAc = row.dataset.uniprotAc || "";
        const mane = row.dataset.mane || "";
        const codon = row.dataset.codon || "";
        const codonPosition = row.dataset.codonPosition || "";
        const strand = normalizeStrandSymbol(row.dataset.strand);

        const proteinCell = row.querySelector(".protein-position-cell");

        const aaLabel = row.dataset.aminoAcid || "";
        const proteinPos = row.dataset.proteinPos || "";
        const codonLabel = formatCodonWithHighlight(codon, codonPosition, strand);

        proteinCell.innerHTML =
            "p." + aaLabel + proteinPos + (codon.trim() ? " " + codonLabel : "");

        const ensemblTranscriptHref = window.METADOME_LINKS.ensemblTranscript(genomeBuild, gencodeTranscript);
        const refseqHref = refseqTranscript ? window.METADOME_LINKS.refseq(refseqTranscript) : null;

        row.cells[1].innerHTML =
            buildExternalAnchorOrFallback(ensemblTranscriptHref, gencodeTranscript) +
            (refseqTranscript ? " / " + buildExternalAnchorOrFallback(refseqHref, refseqTranscript) : "") +
            (mane === "MANE_Select" ? ' <span class="tag is-info is-light">[MANE Select]</span>' : '');

        const uniprotHref = uniprotAc ? window.METADOME_LINKS.uniprot(uniprotAc) : null;
        row.cells[2].innerHTML = uniprotAc
            ? buildExternalAnchorOrFallback(uniprotHref, uniprotAc)
            : "";

        const pfamCell = row.querySelector(".pfam-domain-cell");
        if (pfamCell) {
            const pfamEntries = Array.from(pfamCell.querySelectorAll(".pfam-domain-entry"));
            if (pfamEntries.length > 0) {
                const renderedEntries = pfamEntries.map(function(entry) {
                    const pfamId = entry.dataset.pfamId || "";
                    const label = entry.dataset.pfamLabel || entry.textContent.trim();
                    const pfamHref = window.METADOME_LINKS.pfam(pfamId);
                    return buildExternalAnchorOrFallback(pfamHref, label);
                });
                pfamCell.innerHTML = renderedEntries.join(", ");
            }
        }
    });

    const headers = table.querySelectorAll("th.sortable");
    let currentSort = { col: null, asc: true };

    headers.forEach(function(th) {
        th.style.cursor = "pointer";
        th.addEventListener("click", function() {
            const sortKey = th.dataset.sort;

            if (currentSort.col === sortKey) {
                currentSort.asc = !currentSort.asc;
            } else {
                currentSort.col = sortKey;
                currentSort.asc = true;
            }

            const tbody = table.querySelector("tbody");
            const rows = Array.from(tbody.querySelectorAll("tr[data-gene]"));

            rows.sort(function(a, b) {
                let valA, valB;

                if (sortKey === "gene") {
                    valA = a.dataset.gene.toLowerCase();
                    valB = b.dataset.gene.toLowerCase();
                } else if (sortKey === "protein") {
                    valA = parseInt(a.dataset.proteinPos, 10) || 0;
                    valB = parseInt(b.dataset.proteinPos, 10) || 0;
                    return currentSort.asc ? valA - valB : valB - valA;
                } else if (sortKey === "transcript") {
                    valA = a.cells[1].textContent.toLowerCase();
                    valB = b.cells[1].textContent.toLowerCase();
                } else {
                    return 0;
                }

                if (valA < valB) return currentSort.asc ? -1 : 1;
                if (valA > valB) return currentSort.asc ? 1 : -1;
                return 0;
            });

            rows.forEach(function(row) {
                tbody.appendChild(row);
            });

            headers.forEach(function(h) {
                h.textContent = h.textContent.replace(/ [▲▼]$/, "");
            });
            th.textContent += currentSort.asc ? " ▲" : " ▼";
        });
    });
});