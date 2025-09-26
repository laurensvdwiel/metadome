from metadome.domain.repositories import GeneRepository
from metadome.domain.services.helper_functions import is_transcript_id
from metadome.controllers.job import retrieve_error
from flask import Blueprint, jsonify, request
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
    # retrieve the transcript ids for this gene
    transcripts = GeneRepository.retrieve_all_transcript_ids(genome_build, gene_name)
    
    # check if there was any return value
    if len(transcripts) > 0:
        message = "Retrieved transcripts for gene '"+transcripts[0].gene_name+"'"
    else:
        message = "No transcripts available in database for gene '"+gene_name+"'"

    transcript_results = []
    for t in transcripts:
        transcript_entry = {}
        transcript_entry['aa_length'] = t.sequence_length
        transcript_entry['gencode_id'] = t.gencode_transcription_id
        transcript_entry['refseq_nm_numbers'] = ", ".join(nm_number for nm_number in t.refseq_transcript_id.split(','))
        transcript_entry['has_protein_data'] = not t.protein_id is None
        transcript_results.append(transcript_entry)

    return jsonify(transcript_ids=transcript_results, message=message, genome_build=genome_build, gene_name=gene_name)

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

    transcript_id = data['transcript_id']
    protein_pos = data['protein_position']
    requested_domains = data['requested_domains']

    _log.debug("get_metadomain_annotation with transcript: {}, protein position: {}, requested_domains: {}"
               .format(transcript_id, protein_pos, requested_domains))

    # attempt to retrieve the response for a metadomain position
    from metadome.tasks import retrieve_metadomain_annotation as rma
    response = rma(transcript_id, protein_pos, requested_domains)

    return jsonify(response)


@bp.route('/submit_visualization/', methods=['POST'])
def submit_visualization_job_for_transcript():
    data = request.get_json()

    if not 'transcript_id' in data:
        return jsonify({'error': "no transcript id"}), 400

    transcript_id = data['transcript_id']

    if not is_transcript_id(transcript_id):
        return jsonify({'error': "not a valid transcript id: {}".format(transcript_id)}), 400

    _log.debug("submitted {}".format(transcript_id))

    create_visualization_job_if_needed(transcript_id)

    # It has to return something :(
    return jsonify({'transcript_id': transcript_id})


@bp.route('/status/<transcript_id>/', methods=['GET'])
def get_visualization_status_for_transcript(transcript_id):
    status = get_visualization_status(transcript_id)

    response = {'status': status}
    return jsonify(response)


@bp.route('/result/<transcript_id>/', methods=['GET'])
def get_visualization_result_for_transcript(transcript_id):
    try:
        result = retrieve_visualization(transcript_id)
        return jsonify(result)
    except FileNotFoundError:
        return jsonify({'error': "file not found"}), 404


@bp.route('/error/<transcript_id>/', methods=['GET'])
def get_visualization_error_for_transcript(transcript_id):
    stacktrace = retrieve_error(transcript_id)
    error = "error running visualization job"
    return jsonify({'error': error, 'stacktrace': stacktrace})


@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    _log.error("Unhandled exception:{}\n{}"
               .format(error, traceback.format_exc()))
    return jsonify({'error': str(error),
                    'stacktrace': traceback.format_exc()}), 500
