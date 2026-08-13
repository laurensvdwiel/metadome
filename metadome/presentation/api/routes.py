from metadome.domain.models.entities.gene_region import GenomeBuild
from metadome.domain.models.entities.single_nucleotide_variant import VariantSource
from metadome.domain.repositories import GeneRepository
from metadome.domain.services.helper_functions import is_transcript_id
from metadome.controllers.job import retrieve_error
from flask import Blueprint, jsonify, request, make_response
from builtins import Exception
import traceback
import logging
from metadome.controllers.job import (create_visualization_job_if_needed,
                                      get_visualization_status,
                                      retrieve_visualization)

_log = logging.getLogger(__name__)

bp = Blueprint('api', __name__)

##########################################
#########    Route end points    #########
##########################################

@bp.route('/get_transcripts/<string:genome_build>/<string:gene_name>', methods=['GET'])
def get_transcript_ids_for_gene(genome_build, gene_name):
    _log.debug('get_transcript_ids_for_gene')
    transcripts = GeneRepository.retrieve_all_transcript_ids(genome_build, gene_name)

    if len(transcripts) > 0:
        message = "Retrieved transcripts for gene '" + transcripts[0]['gene_name'] + "'"
    else:
        message = "No transcripts available in database for gene '" + gene_name + "'"

    transcript_results = []
    for t in transcripts:
        transcript_entry = {}
        transcript_entry['aa_length'] = t['sequence_length']
        transcript_entry['gencode_id'] = t['gencode_transcription_id']
        if t['refseq_transcript_id'] is None:
            transcript_entry['refseq_nm_numbers'] = ""
        else:
            transcript_entry['refseq_nm_numbers'] = ", ".join(
                nm_number for nm_number in t['refseq_transcript_id'].split(',')
            )
        transcript_entry['mane_transcript_type'] = t['mane_transcript_type'] if t['mane_transcript_type'] is not None else ""
        transcript_entry['has_protein_data'] = t['protein_id'] is not None
        transcript_results.append(transcript_entry)

    response = make_response(jsonify(
        transcript_ids=transcript_results,
        message=message,
        genome_build=genome_build,
        gene_name=gene_name
    ))
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response

@bp.route('/get_metadomain_annotation/', methods=['POST'])
def get_metadomain_annotation():
    data = request.get_json()

    _log.debug("data is {}".format(data))

    if not 'transcript_id' in data:
        return jsonify({"error: no transcript id"}), 400
    elif not 'protein_position' in data:
        return jsonify({"error: no protein position"}), 400
    elif not 'requested_domains' in data:
        return jsonify({"error: no requested domains"}), 400
    elif not 'genome_build' in data:
        return jsonify({"error: no genome build"}), 400

    variant_sources = data.get('variant_sources')
    if variant_sources is not None:
        try:
            variant_sources = [VariantSource(source) for source in variant_sources]
        except (TypeError, ValueError):
            return jsonify({'error': "invalid variant sources: {}".format(data.get('variant_sources'))}), 400

    transcript_id = data['transcript_id']
    protein_pos = data['protein_position']
    requested_domains = data['requested_domains']
    genome_build_full = data['genome_build']

    # set the genome build enum
    if genome_build_full.upper()[:6] == 'GRCH37':
        genome_build = GenomeBuild.GRCh37
    elif genome_build_full.upper()[:6] == 'GRCH38':
        genome_build = GenomeBuild.GRCh38
    else:
        return jsonify({"error: at meta domain annotation request: no genome build {}".format(genome_build_full)}), 400

    _log.debug("get_metadomain_annotation with transcript: {}, protein position: {}, requested_domains: {}, genome_build: {}"
               .format(transcript_id, protein_pos, requested_domains, genome_build))

    # attempt to retrieve the response for a metadomain position
    from metadome.tasks import retrieve_metadomain_annotation as rma
    response = rma(transcript_id, protein_pos, requested_domains, genome_build, variant_sources=variant_sources)

    return jsonify(response)


@bp.route('/submit_visualization/', methods=['POST'])
def submit_visualization_job_for_transcript():
    data = request.get_json()

    if not 'transcript_id' in data:
        return jsonify({'error': "no transcript id"}), 400
    if not 'genome_build' in data:
        return jsonify({'error': "no genome build"}), 400

    transcript_id = data['transcript_id']
    genome_build = data['genome_build']

    # check if genome build is valid
    if not genome_build in GeneRepository.retrieve_all_genome_builds_from_db():
        return jsonify({'error': "not a valid genome build: {}".format(genome_build)}), 400


    if not is_transcript_id(transcript_id):
        return jsonify({'error': "for genome build {}, not a valid transcript id: {}".format(genome_build, transcript_id)}), 400


    _log.debug("submitted submit_visualization_job_for_transcript with genome_build {}, transcript_id {}".format(genome_build, transcript_id))

    create_visualization_job_if_needed(transcript_id, genome_build)

    # It has to return something :(
    return jsonify({'transcript_id': transcript_id, 'genome_build': genome_build})


@bp.route('/status/<genome_build>/<transcript_id>', methods=['GET'])
def get_visualization_status_for_transcript(genome_build, transcript_id):
    status = get_visualization_status(transcript_id, genome_build)

    response = {'status': status}
    return jsonify(response)


@bp.route('/result/<genome_build>/<transcript_id>', methods=['GET'])
def get_visualization_result_for_transcript(genome_build, transcript_id):
    try:
        result = retrieve_visualization(transcript_id, genome_build)
        return jsonify(result)
    except FileNotFoundError:
        return jsonify({'error': "file not found"}), 404


@bp.route('/error/<genome_build>/<transcript_id>', methods=['GET'])
def get_visualization_error_for_transcript(genome_build, transcript_id):
    stacktrace, timestamp = retrieve_error(transcript_id, genome_build)
    error = "error running visualization job"
    return jsonify({'error': error, 'stacktrace': stacktrace, 'timestamp': timestamp})


@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    _log.error("Unhandled exception:{}\n{}"
               .format(error, traceback.format_exc()))
    return jsonify({'error': str(error),
                    'stacktrace': traceback.format_exc()}), 500
