document.addEventListener("DOMContentLoaded", function() {
    const table = document.getElementById("position-results-table");
    if (!table || !window.METADOME_LINKS) {
        return;
    }

    const rows = table.querySelectorAll("tbody tr[data-gene]");
    rows.forEach(function(row) {
        const genomeBuild = row.dataset.build || window.METADOME_CONFIG.genomeBuild || "";
        const chr = row.dataset.chr || "";
        const genomicPos = row.dataset.pos || "";
        const proteinPos = row.dataset.proteinPos || "";
        const geneName = row.dataset.gene || "";
        const gencodeTranscript = row.dataset.gencodeTranscript || "";
        const refseqTranscript = row.dataset.refseqTranscript || "";
        const uniprotAc = row.dataset.uniprotAc || "";
        const mane = row.dataset.mane || "";
        const basePair = row.dataset.basePair || "";
        const strand = row.dataset.strand || "";

        const genomicCell = row.cells[2];
        const transcriptCell = row.cells[3];
        const uniprotCell = row.cells[4];
        const metadomeCell = row.cells[6];

        genomicCell.innerHTML =
            '<a href="' + window.METADOME_LINKS.ensemblGenomicPosition(genomeBuild, chr, genomicPos) + '" target="_blank" rel="noopener">' +
            chr + ':' + genomicPos +
            '</a> ' + basePair + ' (' + (strand === "plus" ? "+" : "-") + ')';

        transcriptCell.innerHTML =
            '<a href="' + window.METADOME_LINKS.ensemblTranscript(genomeBuild, gencodeTranscript) + '" target="_blank" rel="noopener">' +
            gencodeTranscript +
            '</a>' +
            (refseqTranscript ? ' / <a href="' + window.METADOME_LINKS.refseq(refseqTranscript) + '" target="_blank" rel="noopener">' + refseqTranscript + '</a>' : '') +
            (mane === "MANE_Select" ? ' <span class="tag is-info is-light">[MANE Select]</span>' : '');

        uniprotCell.innerHTML = uniprotAc
            ? '<a href="' + window.METADOME_LINKS.uniprot(uniprotAc) + '" target="_blank" rel="noopener">' + uniprotAc + '</a>'
            : '';
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
                } else if (sortKey === "build") {
                    valA = a.dataset.build.toLowerCase();
                    valB = b.dataset.build.toLowerCase();
                } else if (sortKey === "genomic") {
                    valA = a.dataset.chr + String(parseInt(a.dataset.pos)).padStart(12, "0");
                    valB = b.dataset.chr + String(parseInt(b.dataset.pos)).padStart(12, "0");
                } else if (sortKey === "protein") {
                    valA = parseInt(a.dataset.proteinPos) || 0;
                    valB = parseInt(b.dataset.proteinPos) || 0;
                    return currentSort.asc ? valA - valB : valB - valA;
                } else if (sortKey === "transcript") {
                    valA = a.cells[3].textContent.toLowerCase();
                    valB = b.cells[3].textContent.toLowerCase();
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