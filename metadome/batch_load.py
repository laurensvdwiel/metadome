import csv
import logging
import traceback
import argparse
import os
import gzip
from collections import Counter, defaultdict

from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func
from sqlalchemy.exc import IntegrityError

from metadome.default_settings import SQLALCHEMY_DATABASE_URI

# --- Flask App and SQLAlchemy Setup ---
# This would typically be in your app.py or a config file
flask_app = Flask(__name__)

# IMPORTANT: Configure your actual database URI
flask_app.config['SQLALCHEMY_DATABASE_URI'] = SQLALCHEMY_DATABASE_URI #f"postgresql://{DB_USER}:{DB_PWD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
flask_app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # Optional: silence a warning

db = SQLAlchemy() # Initialize your db object globally or pass around

# --- Model Definitions (gene.py, protein.py, etc.) ---
# (These files would import 'db' from wherever you defined it above)
from metadome.database import db
from metadome.domain.models.gene import Gene
from metadome.domain.models.protein import Protein, ProteinSource 
from metadome.domain.models.mapping import Mapping
from metadome.domain.models.interpro import Interpro
from metadome.domain.models.meta_domain_position import MetaDomainPosition
from metadome.domain.models.meta_domain_mapping import MetaDomainMapping

# --- Helper Functions (safe_int_conversion, safe_bool_conversion, get_cleaned_str) ---
# (These remain unchanged)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def safe_int_conversion(value_str, default=None):
    if value_str is None or value_str.strip() == '' or value_str.lower() == 'na' or value_str.lower() == '-':
        return default
    try:
        return int(value_str)
    except ValueError:
        logging.warning(f"Could not convert '{value_str}' to int, using default: {default}")
        return default

def safe_bool_conversion(value_str, default=False):
    if value_str is None or value_str.strip() == '':
        return default
    val_lower = value_str.lower()
    if val_lower in ['true', 't', '1', 'yes', 'y']:
        return True
    if val_lower in ['false', 'f', '0', 'no', 'n', 'na', '-']:
        return False
    logging.warning(f"Could not convert '{value_str}' to bool, using default: {default}")
    return default

def get_cleaned_str(row, key, default=None):
    val = row.get(key)
    if val is None: return default
    val_stripped = val.strip()
    return val_stripped if val_stripped else default

def _make_gene_key(row):
    return (
        get_cleaned_str(row, 'gencode_transcription_id'),
        get_cleaned_str(row, 'GENCODE_version'),
        get_cleaned_str(row, 'genome_build'),
    )

def _format_gene_key(gene_key):
    transcript_id, gencode_version, genome_build = gene_key
    return f"{transcript_id} | {gencode_version} | {genome_build}"

def _normalize_gene_symbol(value):
    if value is None:
        return None
    return str(value).strip().upper()

def _gene_names_match(gencode_gene_name, uniprot_gene_name):
    left = _normalize_gene_symbol(gencode_gene_name)
    right = _normalize_gene_symbol(uniprot_gene_name)
    if not left or not right:
        return False
    return left == right

def _new_report():
    return {
        'total_rows': 0,
        'rows_loaded': 0,
        'rows_skipped': 0,
        'proteins_created': 0,
        'genes_created': 0,
        'interpro_created': 0,
        'mappings_created': 0,
        'meta_domain_positions_created': 0,
        'meta_domain_mappings_created': 0,
        'protein_conflict_rows_skipped': 0,
        'non_canonical_uniprot_rows_skipped': 0,
        'invalid_rows': 0,
        'integrity_error_rows': 0,
        'value_error_rows': 0,
        'unexpected_error_rows': 0,
        'warnings': Counter(),
        'examples': defaultdict(list),
        'gene_conflict_details': defaultdict(lambda: {
            'selected_uniprot': None,
            'selected_gene_name': None,
            'conflicting_uniprots': Counter(),
        }),
        'gene_name_mismatch_rows': 0,
        'gene_name_mismatch_summary': defaultdict(lambda: {
            'count': 0,
            'uniprot_acs': Counter(),
        }),
    }

def _add_example(report, key, message, limit=10):
    if len(report['examples'][key]) < limit:
        report['examples'][key].append(message)

def print_final_report(csv_report, load_report, db_report):
    logging.info("========== BATCH LOAD FINAL REPORT ==========")

    logging.info("Conclusion:")
    logging.info("  loader = integrity-preserving + reporting")
    logging.info("  CSV generator = canonical transcript→Swiss-Prot resolver")
    logging.info("If needed, the conflict report below can be used to fix upstream CSV generation.")

    logging.info("--- CSV verification ---")
    for key in [
        'total_rows',
        'rows_with_gene_key',
        'rows_with_uniprot',
        'rows_with_mapping_minimum_fields',
        'rows_with_interpro_fields',
        'rows_with_metadomain_fields',
        'transcript_keys_seen',
        'transcript_keys_with_multiple_uniprots',
        'transcript_keys_with_unresolved_multiple_uniprots',
        'invalid_rows',
    ]:
        logging.info("%s: %s", key, csv_report.get(key))

    if csv_report.get('examples'):
        for category, examples in csv_report['examples'].items():
            logging.info("CSV example category: %s", category)
            for example in examples:
                logging.info("  %s", example)

    logging.info("--- Load report ---")
    for key in [
        'total_rows',
        'rows_loaded',
        'rows_skipped',
        'proteins_created',
        'genes_created',
        'interpro_created',
        'mappings_created',
        'meta_domain_positions_created',
        'meta_domain_mappings_created',
        'protein_conflict_rows_skipped',
        'non_canonical_uniprot_rows_skipped',
        'gene_name_mismatch_rows',
        'invalid_rows',
        'value_error_rows',
        'integrity_error_rows',
        'unexpected_error_rows',
    ]:
        logging.info("%s: %s", key, load_report.get(key))

    if load_report.get('warnings'):
        for warning_key, count in load_report['warnings'].most_common():
            logging.info("warning[%s]=%s", warning_key, count)

    if load_report.get('examples'):
        for category, examples in load_report['examples'].items():
            logging.info("LOAD example category: %s", category)
            for example in examples:
                logging.info("  %s", example)

    if load_report.get('gene_name_mismatch_summary'):
        logging.info("--- Aggregated GENCODE gene_name vs UniProt GN mismatches (first 50) ---")
        shown = 0
        sorted_items = sorted(
            load_report['gene_name_mismatch_summary'].items(),
            key=lambda item: item[1]['count'],
            reverse=True
        )
        for mismatch_key, detail in sorted_items:
            if shown >= 50:
                break
            transcript_id, gencode_gene_name, uniprot_gene_name = mismatch_key
            uniprot_ac_summary = ", ".join(
                f"{uniprot_ac} (rows={count})"
                for uniprot_ac, count in detail['uniprot_acs'].most_common()
            )
            logging.info(
                "  transcript=%s | gencode_gene_name=%s | uniprot_gene_name=%s | rows=%s | proteins=%s",
                transcript_id,
                gencode_gene_name,
                uniprot_gene_name,
                detail['count'],
                uniprot_ac_summary
            )
            shown += 1

    if load_report.get('gene_conflict_details'):
        logging.info("--- Actionable transcript→UniProt conflict summary (first 50) ---")
        shown = 0
        for gene_key, detail in load_report['gene_conflict_details'].items():
            if shown >= 50:
                break
            conflicts = ", ".join(
                f"{uniprot_ac} (skipped_rows={count})"
                for uniprot_ac, count in detail['conflicting_uniprots'].most_common()
            )
            logging.info(
                "  %s | gene=%s | selected=%s | skipped=%s",
                _format_gene_key(gene_key),
                detail['selected_gene_name'],
                detail['selected_uniprot'],
                conflicts
            )
            shown += 1

    logging.info("--- Database verification ---")
    for key in [
        'proteins_total',
        'genes_total',
        'mappings_total',
        'genes_without_protein',
        'mappings_without_gene',
        'mappings_without_protein',
        'genes_with_multiple_mapping_proteins',
        'transcript_keys_with_multiple_mapping_uniprots',
    ]:
        logging.info("%s: %s", key, db_report.get(key))

    if db_report.get('examples'):
        for category, examples in db_report['examples'].items():
            logging.info("DB example category: %s", category)
            for example in examples:
                logging.info("  %s", example)

    logging.info("========== END OF REPORT ==========")

# --- Batch Load Function (uses the passed sqlalchemy_session, which will be db.session) ---
def batch_load_data(csv_filepath, sqlalchemy_session, csv_report=None, batch_size=5000):
    protein_cache = {}
    gene_cache = {}
    interpro_cache = {}
    metadomain_position_cache = {}
    load_report = _new_report()
    canonical_uniprot_by_gene_key = {}
    canonical_uniprot_selection_reason = {}

    if csv_report is not None:
        canonical_uniprot_by_gene_key = csv_report.get('canonical_uniprot_by_gene_key', {})
        canonical_uniprot_selection_reason = csv_report.get('canonical_uniprot_selection_reason', {})

    try:
        with gzip.open(csv_filepath, mode='rt') as infile:
            reader = csv.DictReader(infile)
            current_batch_count = 0

            for i, row in enumerate(reader):
                line_num = i + 2
                load_report['total_rows'] += 1
                current_protein = None
                current_gene = None
                current_interpro = None
                current_mapping = None

                try:
                    # TODO(batch-integrity): a single bad row currently rolls back the whole uncommitted batch (up to batch_size good rows) and leaves counters + caches stale.
                    #  Fix: wrap this row body in `with sqlalchemy_session.begin_nested()` (per-row SAVEPOINT) so only the failing row unwinds; move the batch-commit block below out of the savepoint; then in each except handler: (a) restore CREATION_KEYS counters from a per-row snapshot, (b) _clear_caches(), (c) DROP the full rollback() call. See handlers at Exception Catches ~600/604/612.
                    gencode_tr_id = get_cleaned_str(row, 'gencode_transcription_id')
                    gencode_tr_version = get_cleaned_str(row, 'GENCODE_version')
                    genome_build = get_cleaned_str(row, 'genome_build')
                    current_gene_key = (gencode_tr_id, gencode_tr_version, genome_build)
                    uniprot_ac = get_cleaned_str(row, 'uniprot_ac')

                    selected_uniprot_ac = canonical_uniprot_by_gene_key.get(current_gene_key)
                    if selected_uniprot_ac is not None and uniprot_ac != selected_uniprot_ac:
                        load_report['rows_skipped'] += 1
                        load_report['protein_conflict_rows_skipped'] += 1
                        load_report['non_canonical_uniprot_rows_skipped'] += 1
                        load_report['warnings']['non_canonical_uniprot_skipped'] += 1

                        conflict_detail = load_report['gene_conflict_details'][current_gene_key]
                        conflict_detail['selected_uniprot'] = selected_uniprot_ac
                        conflict_detail['selected_gene_name'] = get_cleaned_str(row, 'gene_name')
                        conflict_detail['conflicting_uniprots'][uniprot_ac] += 1

                        _add_example(
                            load_report,
                            'non_canonical_uniprot_skipped',
                            f"L{line_num}: {_format_gene_key(current_gene_key)} selected "
                            f"'{selected_uniprot_ac}' reason="
                            f"'{canonical_uniprot_selection_reason.get(current_gene_key)}' "
                            f"skipping row protein '{uniprot_ac}'"
                        )
                        continue

                    # 1. Process Protein
                    uniprot_ac = get_cleaned_str(row, 'uniprot_ac')
                    if uniprot_ac:
                        if uniprot_ac in protein_cache:
                            current_protein = protein_cache[uniprot_ac]
                        else:
                            current_protein = sqlalchemy_session.query(Protein).filter_by(
                                uniprot_ac=uniprot_ac).one_or_none()
                            if current_protein is None:
                                uniprot_name_val = get_cleaned_str(row, 'uniprot_name')
                                uniprot_gn_val = get_cleaned_str(row, 'uniprot_gene_name')
                                source_val = get_cleaned_str(row, 'source')
                                source_val_str = source_val.lower() if source_val else None
                                eval_interpro_val = safe_bool_conversion(row.get('evaluated_interpro_domains'), False)

                                if not source_val_str:
                                    load_report['rows_skipped'] += 1
                                    load_report['invalid_rows'] += 1
                                    load_report['warnings']['missing_protein_source'] += 1
                                    _add_example(load_report, 'missing_protein_source',
                                                 f"L{line_num}: missing source for Protein '{uniprot_ac}'")
                                    protein_cache[uniprot_ac] = None
                                    continue

                                current_protein = Protein(
                                    _uniprot_ac=uniprot_ac,
                                    _uniprot_name=uniprot_name_val,
                                    _source=source_val_str,
                                    _uniprot_gn=uniprot_gn_val
                                )
                                current_protein.evaluated_interpro_domains = eval_interpro_val
                                sqlalchemy_session.add(current_protein)
                                load_report['proteins_created'] += 1
                            else:
                                if not current_protein.uniprot_gn:
                                    current_protein.uniprot_gn = get_cleaned_str(row, 'uniprot_gene_name')

                            protein_cache[uniprot_ac] = current_protein

                    # 2. Process Gene
                    gencode_tr_id = get_cleaned_str(row, 'gencode_transcription_id')
                    gencode_tr_version = get_cleaned_str(row, 'GENCODE_version')
                    genome_build = get_cleaned_str(row, 'genome_build')
                    current_gene_key = (gencode_tr_id, gencode_tr_version, genome_build)

                    if gencode_tr_id and gencode_tr_version and genome_build:
                        if current_gene_key in gene_cache:
                            current_gene = gene_cache[current_gene_key]
                        else:
                            current_gene = sqlalchemy_session.query(Gene).filter_by(
                                gencode_transcription_id=gencode_tr_id,
                                gencode_version=gencode_tr_version,
                                genome_build=genome_build
                            ).one_or_none()

                            if current_gene is None:
                                strand_str = get_cleaned_str(row, 'strand')
                                gene_name_val = get_cleaned_str(row, 'gene_name')
                                gencode_transl_name_val = get_cleaned_str(row, 'gencode_translation_name')

                                if not strand_str or strand_str not in ['+', '-']:
                                    load_report['rows_skipped'] += 1
                                    load_report['invalid_rows'] += 1
                                    load_report['warnings']['invalid_gene_strand'] += 1
                                    _add_example(load_report, 'invalid_gene_strand',
                                                 f"L{line_num}: invalid strand for Gene '{gencode_tr_id}'")
                                    gene_cache[current_gene_key] = None
                                    continue

                                if not gencode_transl_name_val:
                                    load_report['rows_skipped'] += 1
                                    load_report['invalid_rows'] += 1
                                    load_report['warnings']['missing_translation_name'] += 1
                                    _add_example(load_report, 'missing_translation_name',
                                                 f"L{line_num}: missing gencode_translation_name for Gene '{gencode_tr_id}'")
                                    gene_cache[current_gene_key] = None
                                    continue

                                current_gene = Gene(
                                    _strand=strand_str,
                                    _gene_name=gene_name_val,
                                    _gencode_transcription_id=gencode_tr_id,
                                    _gencode_translation_name=gencode_transl_name_val,
                                    _gencode_gene_id=get_cleaned_str(row, 'gencode_gene_id'),
                                    _gencode_version=gencode_tr_version,
                                    _gencode_basic=safe_bool_conversion(row.get('GencodeBasic'), False),
                                    _genome_build=genome_build,
                                    _refseq_transcript_id=get_cleaned_str(row, 'RefSeq'),
                                    _havana_gene_id=get_cleaned_str(row, 'havana_gene_id'),
                                    _havana_translation_id=get_cleaned_str(row, 'havana_translation_id'),
                                    _mane_transcript_type=get_cleaned_str(row, 'MANE'),
                                    _sequence_length=safe_int_conversion(row.get('sequence_length'))
                                )

                                if current_protein is not None:
                                    current_gene.protein = current_protein

                                sqlalchemy_session.add(current_gene)
                                load_report['genes_created'] += 1

                            gene_cache[current_gene_key] = current_gene
                    else:
                        load_report['rows_skipped'] += 1
                        load_report['invalid_rows'] += 1
                        load_report['warnings']['missing_gene_key'] += 1
                        _add_example(load_report, 'missing_gene_key',
                                     f"L{line_num}: missing gene key for transcript '{gencode_tr_id}'")
                        continue

                    if current_gene is None or current_protein is None:
                        load_report['rows_skipped'] += 1
                        load_report['warnings']['missing_gene_or_protein'] += 1
                        _add_example(load_report, 'missing_gene_or_protein',
                                     f"L{line_num}: missing Gene or Protein for transcript '{gencode_tr_id}', protein '{uniprot_ac}'")
                        continue

                    uniprot_gn_val = get_cleaned_str(row, 'uniprot_gene_name')
                    if uniprot_gn_val and not _gene_names_match(current_gene.gene_name, uniprot_gn_val):
                        load_report['gene_name_mismatch_rows'] += 1
                        mismatch_key = (
                            current_gene.gencode_transcription_id,
                            current_gene.gene_name,
                            uniprot_gn_val
                        )
                        load_report['gene_name_mismatch_summary'][mismatch_key]['count'] += 1
                        load_report['gene_name_mismatch_summary'][mismatch_key]['uniprot_acs'][
                            current_protein.uniprot_ac] += 1

                        _add_example(
                            load_report,
                            'gene_name_mismatch',
                            f"L{line_num}: transcript '{current_gene.gencode_transcription_id}' "
                            f"gencode_gene_name='{current_gene.gene_name}' "
                            f"uniprot_gene_name='{uniprot_gn_val}' "
                            f"uniprot_ac='{current_protein.uniprot_ac}'"
                        )

                    # Enforce 1 transcript/gene -> 1 protein strictly
                    if current_gene.protein is None:
                        current_gene.protein = current_protein
                    elif current_gene.protein.uniprot_ac != current_protein.uniprot_ac:
                        # Canonical selection has already happened during CSV verification.
                        if current_gene.protein is None:
                            current_gene.protein = current_protein
                        elif current_gene.protein.uniprot_ac != current_protein.uniprot_ac:
                            load_report['rows_skipped'] += 1
                            load_report['protein_conflict_rows_skipped'] += 1
                            load_report['warnings']['gene_protein_conflict_after_canonical_selection'] += 1

                            conflict_detail = load_report['gene_conflict_details'][current_gene_key]
                            conflict_detail['selected_uniprot'] = current_gene.protein.uniprot_ac
                            conflict_detail['selected_gene_name'] = current_gene.gene_name
                            conflict_detail['conflicting_uniprots'][current_protein.uniprot_ac] += 1

                            _add_example(
                                load_report,
                                'gene_protein_conflict_after_canonical_selection',
                                f"L{line_num}: {_format_gene_key(current_gene_key)} selected "
                                f"'{current_gene.protein.uniprot_ac}' conflicts with row protein "
                                f"'{current_protein.uniprot_ac}' even after canonical preselection"
                            )
                            continue

                    # 3. Process Interpro
                    ext_db_id_val = get_cleaned_str(row, 'ext_db_id')
                    ext_db_version_val = get_cleaned_str(row, 'PFAM_version')
                    uniprot_start_val = safe_int_conversion(row.get('uniprot_start'))
                    uniprot_stop_val = safe_int_conversion(row.get('uniprot_stop'))

                    if (
                        ext_db_id_val is not None and
                        ext_db_version_val is not None and
                        uniprot_start_val is not None and
                        uniprot_stop_val is not None
                    ):
                        interpro_cache_key = (
                            current_protein.uniprot_ac,
                            ext_db_id_val,
                            ext_db_version_val,
                            uniprot_start_val,
                            uniprot_stop_val
                        )
                        if interpro_cache_key not in interpro_cache:
                            current_interpro = None
                            if current_protein.id:
                                current_interpro = sqlalchemy_session.query(Interpro).filter_by(
                                    protein_id=current_protein.id,
                                    ext_db_id=ext_db_id_val,
                                    ext_db_version=ext_db_version_val,
                                    uniprot_start=uniprot_start_val,
                                    uniprot_stop=uniprot_stop_val
                                ).one_or_none()

                            if not current_interpro:
                                current_interpro = Interpro(
                                    _ext_db_id=ext_db_id_val,
                                    _ext_db_version=ext_db_version_val,
                                    _start_pos=uniprot_start_val,
                                    _end_pos=uniprot_stop_val,
                                    _interpro_id=get_cleaned_str(row, 'interpro_id'),
                                    _region_name=get_cleaned_str(row, 'region_name')
                                )
                                current_interpro.protein = current_protein
                                sqlalchemy_session.add(current_interpro)
                                load_report['interpro_created'] += 1

                            interpro_cache[interpro_cache_key] = current_interpro
                        else:
                            current_interpro = interpro_cache[interpro_cache_key]

                    # 4. Process Mapping
                    chromosome_val = get_cleaned_str(row, 'chromosome')
                    chromosome_pos_val = safe_int_conversion(row.get('chromosome_position'))
                    map_strand_enum = current_gene.strand

                    if not chromosome_val or chromosome_pos_val is None:
                        load_report['rows_skipped'] += 1
                        load_report['invalid_rows'] += 1
                        load_report['warnings']['missing_mapping_coordinates'] += 1
                        _add_example(load_report, 'missing_mapping_coordinates',
                                     f"L{line_num}: missing/invalid mapping coordinates for gene '{current_gene.gencode_transcription_id}'")
                        continue

                    current_mapping = Mapping(
                        chromosome=chromosome_val,
                        chromosome_position=chromosome_pos_val,
                        strand=map_strand_enum,
                        base_pair=get_cleaned_str(row, 'base_pair'),
                        codon=get_cleaned_str(row, 'codon'),
                        codon_base_pair_position=safe_int_conversion(row.get('codon_base_pair_position')),
                        cDNA_position=safe_int_conversion(row.get('cDNA_position')),
                        uniprot_residue=get_cleaned_str(row, 'uniprot_residue'),
                        uniprot_position=safe_int_conversion(row.get('uniprot_position')),
                        exon_number=safe_int_conversion(row.get('exon_number'))
                    )
                    current_mapping.gene = current_gene
                    current_mapping.protein = current_protein
                    sqlalchemy_session.add(current_mapping)
                    load_report['mappings_created'] += 1
                    load_report['rows_loaded'] += 1

                    if current_interpro is not None and current_mapping is not None:
                        pfam_consensus_pos_val = safe_int_conversion(row.get('PFAM_consensus_pos'))
                        pfam_consensus_length_val = safe_int_conversion(row.get('PFAM_consensus_length'))

                        if pfam_consensus_pos_val is not None:
                            metadomain_position_cache_key = (ext_db_id_val, pfam_consensus_pos_val)
                            current_meta_domain_position = None

                            if metadomain_position_cache_key not in metadomain_position_cache:
                                current_meta_domain_position = sqlalchemy_session.query(MetaDomainPosition).filter_by(
                                    consensus_position=pfam_consensus_pos_val,
                                    ext_db_id=ext_db_id_val
                                ).one_or_none()

                                if not current_meta_domain_position:
                                    current_meta_domain_position = MetaDomainPosition(
                                        consensus_position=pfam_consensus_pos_val,
                                        consensus_length=pfam_consensus_length_val,
                                        ext_db_id=current_interpro.ext_db_id
                                    )
                                    sqlalchemy_session.add(current_meta_domain_position)
                                    load_report['meta_domain_positions_created'] += 1

                                metadomain_position_cache[metadomain_position_cache_key] = current_meta_domain_position
                            else:
                                current_meta_domain_position = metadomain_position_cache[metadomain_position_cache_key]

                            if current_meta_domain_position is not None:
                                current_meta_domain_mapping = MetaDomainMapping()
                                current_meta_domain_mapping.meta_domain_position = current_meta_domain_position
                                current_meta_domain_mapping.mapping = current_mapping
                                current_meta_domain_mapping.interpro_domain = current_interpro
                                sqlalchemy_session.add(current_meta_domain_mapping)
                                load_report['meta_domain_mappings_created'] += 1

                    current_batch_count += 1 # TODO(batch-integrity): this bookkeeping must sit AFTER the begin_nested() block, not inside it.
                    if current_batch_count >= batch_size:
                        logging.info(
                            "Processed %s rows. Committing batch of %s.",
                            load_report['total_rows'],
                            current_batch_count
                        )
                        sqlalchemy_session.commit()
                        protein_cache.clear()
                        gene_cache.clear()
                        interpro_cache.clear()
                        metadomain_position_cache.clear()
                        current_batch_count = 0

                except ValueError as ve:
                    load_report['rows_skipped'] += 1
                    load_report['value_error_rows'] += 1
                    logging.error(f"L{line_num}: Data validation error: {ve}. Row: {row}. Skipping.")
                    sqlalchemy_session.rollback() # TODO(batch-integrity): remove this full rollback (savepoint handles it); instead restore counters_before + _clear_caches() so this row's counts and cache refs don't leak.
                except IntegrityError as ie:
                    load_report['rows_skipped'] += 1
                    load_report['integrity_error_rows'] += 1
                    logging.error(f"L{line_num}: Database integrity error: {ie}. Row: {row}. Rolling back.")
                    sqlalchemy_session.rollback() # TODO(batch-integrity): remove this full rollback (savepoint handles it); instead restore counters_before + _clear_caches() so this row's counts and cache refs don't leak.
                except Exception as e:
                    load_report['rows_skipped'] += 1
                    load_report['unexpected_error_rows'] += 1
                    logging.error(f"L{line_num}: Unexpected row error: {type(e).__name__}: {e}")
                    logging.error("Row: %s", row)
                    logging.error("Traceback: %s", traceback.format_exc())
                    sqlalchemy_session.rollback() # TODO(batch-integrity): remove this full rollback (savepoint handles it); instead restore counters_before + _clear_caches() so this row's counts and cache refs don't leak.

            if current_batch_count > 0:
                logging.info(
                    "Processed %s rows. Committing final batch of %s.",
                    load_report['total_rows'],
                    current_batch_count
                )
                sqlalchemy_session.commit()

            logging.info("Successfully processed %s rows from %s", load_report['total_rows'], csv_filepath)
            return load_report

    except FileNotFoundError:
        logging.error(f"CSV file not found: {csv_filepath}")
        raise
    except Exception as e:
        error_type = type(e).__name__
        error_traceback = traceback.format_exc()
        logging.error(f"Unexpected fatal error: {error_type}: {e}")
        logging.error(f"Traceback: {error_traceback}")
        if sqlalchemy_session:
            sqlalchemy_session.rollback()
        raise

def verify_csv_data_completeness(csv_filepath):
    """
    Verify that the CSV file contents are complete, so that:
    All mappings have a valid gene and protein, and the protein
     has all required columns.

     Finally return the results to be checked for database completeness check.
    :param csv_filepath: Path to the CSV file to verify.
    :return: A dictionary with counts of valid and invalid entries.
    """
    report = {
        'total_rows': 0,
        'rows_with_gene_key': 0,
        'rows_with_uniprot': 0,
        'rows_with_mapping_minimum_fields': 0,
        'rows_with_interpro_fields': 0,
        'rows_with_metadomain_fields': 0,
        'transcript_keys_seen': 0,
        'transcript_keys_with_multiple_uniprots': 0,
        'transcript_keys_with_unresolved_multiple_uniprots': 0,
        'invalid_rows': 0,
        'examples': defaultdict(list),
        'transcript_to_uniprots': defaultdict(set),
        'transcript_conflict_details': defaultdict(lambda: {
            'gencode_gene_names': Counter(),
            'uniprots': defaultdict(lambda: {
                'count': 0,
                'protein_names': set(),
                'uniprot_gene_names': Counter(),
                'first_line': None,
            })
        }),
        'canonical_uniprot_by_gene_key': {},
        'canonical_uniprot_selection_reason': {},
    }

    required_gene_fields = [
        'gencode_transcription_id',
        'GENCODE_version',
        'genome_build',
        'strand',
        'gencode_translation_name',
    ]

    with gzip.open(csv_filepath, mode='rt') as infile:
        reader = csv.DictReader(infile)

        for i, row in enumerate(reader):
            line_num = i + 2
            report['total_rows'] += 1

            gene_key = _make_gene_key(row)
            uniprot_ac = get_cleaned_str(row, 'uniprot_ac')
            uniprot_name = get_cleaned_str(row, 'uniprot_name')
            uniprot_gene_name = get_cleaned_str(row, 'uniprot_gene_name')
            gencode_gene_name = get_cleaned_str(row, 'gene_name')

            missing_gene_fields = [field for field in required_gene_fields if not get_cleaned_str(row, field)]
            if not missing_gene_fields:
                report['rows_with_gene_key'] += 1
            else:
                report['invalid_rows'] += 1
                if len(report['examples']['missing_gene_fields']) < 10:
                    report['examples']['missing_gene_fields'].append(
                        f"L{line_num}: missing gene fields {missing_gene_fields}"
                    )

            if uniprot_ac:
                report['rows_with_uniprot'] += 1

            chromosome = get_cleaned_str(row, 'chromosome')
            chromosome_position = safe_int_conversion(row.get('chromosome_position'))
            cdna_position = safe_int_conversion(row.get('cDNA_position'))
            exon_number = safe_int_conversion(row.get('exon_number'))

            if chromosome and chromosome_position is not None and cdna_position is not None and exon_number is not None:
                report['rows_with_mapping_minimum_fields'] += 1

            ext_db_id = get_cleaned_str(row, 'ext_db_id')
            pfam_version = get_cleaned_str(row, 'PFAM_version')
            uniprot_start = safe_int_conversion(row.get('uniprot_start'))
            uniprot_stop = safe_int_conversion(row.get('uniprot_stop'))

            if ext_db_id and pfam_version and uniprot_start is not None and uniprot_stop is not None:
                report['rows_with_interpro_fields'] += 1

            pfam_consensus_pos = safe_int_conversion(row.get('PFAM_consensus_pos'))
            pfam_consensus_length = safe_int_conversion(row.get('PFAM_consensus_length'))
            if pfam_consensus_pos is not None and pfam_consensus_length is not None:
                report['rows_with_metadomain_fields'] += 1

            if all(gene_key) and uniprot_ac:
                report['transcript_to_uniprots'][gene_key].add(uniprot_ac)

                conflict_details = report['transcript_conflict_details'][gene_key]
                if gencode_gene_name:
                    conflict_details['gencode_gene_names'][gencode_gene_name] += 1

                conflict_entry = conflict_details['uniprots'][uniprot_ac]
                conflict_entry['count'] += 1
                if uniprot_name:
                    conflict_entry['protein_names'].add(uniprot_name)
                if uniprot_gene_name:
                    conflict_entry['uniprot_gene_names'][uniprot_gene_name] += 1
                if conflict_entry['first_line'] is None:
                    conflict_entry['first_line'] = line_num

    report['transcript_keys_seen'] = len(report['transcript_to_uniprots'])

    unresolved_conflicts = []

    for gene_key, uniprots in report['transcript_to_uniprots'].items():
        details = report['transcript_conflict_details'][gene_key]
        gencode_gene_name = None
        if details['gencode_gene_names']:
            gencode_gene_name = details['gencode_gene_names'].most_common(1)[0][0]

        if len(uniprots) == 1:
            selected_uniprot = next(iter(uniprots))
            report['canonical_uniprot_by_gene_key'][gene_key] = selected_uniprot
            report['canonical_uniprot_selection_reason'][gene_key] = 'unique_uniprot'
            continue

        matching_uniprots = []
        for uniprot_ac in uniprots:
            uniprot_gene_names = details['uniprots'][uniprot_ac]['uniprot_gene_names']
            if any(_gene_names_match(gencode_gene_name, uniprot_gn) for uniprot_gn in uniprot_gene_names.keys()):
                matching_uniprots.append(uniprot_ac)

        if len(matching_uniprots) == 1:
            selected_uniprot = matching_uniprots[0]
            report['canonical_uniprot_by_gene_key'][gene_key] = selected_uniprot
            report['canonical_uniprot_selection_reason'][gene_key] = 'unique_uniprot_gn_match'
        else:
            sorted_candidates = sorted(
                uniprots,
                key=lambda ac: (
                    details['uniprots'][ac]['first_line'] if details['uniprots'][ac]['first_line'] is not None else 10 ** 18,
                    ac
                )
            )
            selected_uniprot = sorted_candidates[0]
            report['canonical_uniprot_by_gene_key'][gene_key] = selected_uniprot
            report['canonical_uniprot_selection_reason'][gene_key] = (
                'unresolved_multiple_uniprots_first_seen'
                if len(matching_uniprots) == 0
                else 'unresolved_multiple_uniprot_gn_matches_first_seen'
            )
            unresolved_conflicts.append((gene_key, sorted(uniprots), matching_uniprots))

    conflicting_transcripts = [
        (gene_key, sorted(uniprots))
        for gene_key, uniprots in report['transcript_to_uniprots'].items()
        if len(uniprots) > 1
    ]
    report['transcript_keys_with_multiple_uniprots'] = len(conflicting_transcripts)
    report['transcript_keys_with_unresolved_multiple_uniprots'] = len(unresolved_conflicts)

    for gene_key, uniprots in conflicting_transcripts[:20]:
        details = report['transcript_conflict_details'][gene_key]['uniprots']
        selected_uniprot = report['canonical_uniprot_by_gene_key'].get(gene_key)
        selection_reason = report['canonical_uniprot_selection_reason'].get(gene_key)

        formatted = []
        for uniprot_ac in uniprots:
            detail = details[uniprot_ac]
            protein_names = sorted(detail['protein_names'])
            uniprot_gene_names = [
                f"{gn} (rows={count})"
                for gn, count in detail['uniprot_gene_names'].most_common()
            ]
            formatted.append(
                f"{uniprot_ac} ({', '.join(protein_names) if protein_names else 'unknown'}), "
                f"GN={'; '.join(uniprot_gene_names) if uniprot_gene_names else 'unknown'}, "
                f"rows={detail['count']}, first_line={detail['first_line']}"
            )

        report['examples']['multiple_uniprots_per_transcript'].append(
            f"{_format_gene_key(gene_key)} -> selected={selected_uniprot} "
            f"reason={selection_reason}; candidates: {'; '.join(formatted)}"
        )

    for gene_key, uniprots, matching_uniprots in unresolved_conflicts[:20]:
        report['examples']['unresolved_multiple_uniprots'].append(
            f"{_format_gene_key(gene_key)} -> candidates={uniprots}, "
            f"gn_matching_candidates={matching_uniprots}, "
            f"selected={report['canonical_uniprot_by_gene_key'].get(gene_key)}, "
            f"reason={report['canonical_uniprot_selection_reason'].get(gene_key)}"
        )

    return report

def verify_database_completeness(db_session, csv_completeness_results):
    """
    Verify that the database contains all expected entries.
    This function should check that all proteins, genes, and mappings
    are correctly linked and have the required fields populated.

    :param db_session: The SQLAlchemy session to use for querying the database.
    :return: A dictionary with counts of valid and invalid entries.
    """
    report = {
        'proteins_total': db_session.query(func.count(Protein.id)).scalar(),
        'genes_total': db_session.query(func.count(Gene.id)).scalar(),
        'mappings_total': db_session.query(func.count(Mapping.id)).scalar(),
        'genes_without_protein': db_session.query(func.count(Gene.id)).filter(Gene.protein_id.is_(None)).scalar(),
        'mappings_without_gene': db_session.query(func.count(Mapping.id)).filter(Mapping.gene_id.is_(None)).scalar(),
        'mappings_without_protein': db_session.query(func.count(Mapping.id)).filter(Mapping.protein_id.is_(None)).scalar(),
        'genes_with_multiple_mapping_proteins': 0,
        'transcript_keys_with_multiple_mapping_uniprots': 0,
        'examples': defaultdict(list),
    }

    multi_mapping_protein_rows = db_session.query(
        Gene.id,
        Gene.gene_name,
        Gene.gencode_transcription_id,
        Gene.gencode_version,
        Gene.genome_build,
        func.count(func.distinct(Mapping.protein_id)).label('distinct_mapping_proteins')
    ).join(
        Mapping, Mapping.gene_id == Gene.id
    ).group_by(
        Gene.id,
        Gene.gene_name,
        Gene.gencode_transcription_id,
        Gene.gencode_version,
        Gene.genome_build
    ).having(
        func.count(func.distinct(Mapping.protein_id)) > 1
    ).all()

    report['genes_with_multiple_mapping_proteins'] = len(multi_mapping_protein_rows)

    for row in multi_mapping_protein_rows[:20]:
        report['examples']['genes_with_multiple_mapping_proteins'].append(
            f"{row.gencode_transcription_id} | {row.gencode_version} | {row.genome_build} -> {row.distinct_mapping_proteins}"
        )

    transcript_multi_uniprot_rows = db_session.query(
        Gene.gencode_transcription_id,
        Gene.gencode_version,
        Gene.genome_build,
        func.count(func.distinct(Protein.uniprot_ac)).label('distinct_uniprots')
    ).join(
        Mapping, Mapping.gene_id == Gene.id
    ).join(
        Protein, Protein.id == Mapping.protein_id
    ).group_by(
        Gene.gencode_transcription_id,
        Gene.gencode_version,
        Gene.genome_build
    ).having(
        func.count(func.distinct(Protein.uniprot_ac)) > 1
    ).all()

    report['transcript_keys_with_multiple_mapping_uniprots'] = len(transcript_multi_uniprot_rows)

    for row in transcript_multi_uniprot_rows[:20]:
        report['examples']['transcript_keys_with_multiple_mapping_uniprots'].append(
            f"{row.gencode_transcription_id} | {row.gencode_version} | {row.genome_build} -> {row.distinct_uniprots}"
        )

    if csv_completeness_results.get('transcript_keys_with_multiple_uniprots'):
        report['examples']['csv_precheck_summary'].append(
            f"CSV transcript keys with multiple UniProt accessions: "
            f"{csv_completeness_results['transcript_keys_with_multiple_uniprots']}"
        )

    return report


# --- Main Execution Block ---
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Batch load data from CSV into database.")
    parser.add_argument('--csv', type=str, required=True, help='Path to the CSV file to load.')
    
    args = parser.parse_args()
    
    csv_file = args.csv # Path to your CSV
 
    if not csv_file or not os.path.isfile(csv_file):
        logging.error("CSV file path is required. Use --csv to specify the path.")
        exit(1)

    # init the db with the flask app
    db.init_app(flask_app) 
   
    # This is crucial for standalone scripts using Flask-SQLAlchemy
    with flask_app.app_context():
        # Create tables if they don't exist
        # In a real app, this is often done once or via migrations (e.g., Flask-Migrate)
        db.create_all()

        # Verify CSV data completeness
        logging.info("--- Verifying CSV Data Completeness ---")
        csv_report = verify_csv_data_completeness(csv_file)

        # Perform the batch load using db.session
        load_report = batch_load_data(csv_file, db.session, csv_report=csv_report)

        # Verification of database completeness
        logging.info("--- Verifying Database Completeness ---")
        db_report = verify_database_completeness(db.session, csv_report)

        print_final_report(csv_report, load_report, db_report)

        # Cleanup session through explicit removal
        db.session.remove() # Flask-SQLAlchemy often handles this at context teardown, but for long-running scripts, explicit removal might be good.
