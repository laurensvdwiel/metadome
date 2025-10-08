from flask import Blueprint, g, render_template, redirect, url_for, session, jsonify
from flask_mail import Message
import json
import traceback
from metadome import get_version
from metadome.controllers.job import retrieve_error
from metadome.domain.repositories import GeneRepository
from metadome.domain.services.mail.mail import mail
from metadome.presentation.api.routes import get_transcript_ids_for_gene
from metadome.presentation.web.forms import SupportForm
from metadome.default_settings import DEFAULT_RECIPIENT, MAIL_SERVER

import logging

_log = logging.getLogger(__name__)

bp = Blueprint('web', __name__)

@bp.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@bp.route('/dashboard_js')
def dashboard_js():
    # Renders the javascript used on the index page
    return render_template('/js/dashboard.js')

@bp.route('/dashboard/', methods=['GET'])
def dashboard():
    # Reset the dashboard session to ensure a clean start
    reset_dashboard_session()

    # Redirect to the GRCh38 dashboard by default
    return dashboard_grch38()

@bp.route('/dashboard/tour/', methods=['GET'])
def dashboard_tour():
    # Reset the dashboard session to ensure a clean start
    reset_dashboard_session()

    return render_template('dashboard.html')

@bp.route('/dashboard_grch37/', methods=['GET'])
def dashboard_grch37():
    genome_build_url_safe = dashboard_genome_build_safety_check('GRCh37')

    return redirect(url_for('web.dashboard_genome_build', genome_build=genome_build_url_safe))

@bp.route('/dashboard_grch38/', methods=['GET'])
def dashboard_grch38():
    genome_build_url_safe = dashboard_genome_build_safety_check('GRCh38')

    return redirect(url_for('web.dashboard_genome_build', genome_build=genome_build_url_safe))

@bp.route('/dashboard/<genome_build>/', methods=['GET'])
def dashboard_genome_build(genome_build):
    # Check if the genome build is valid, redirect to default dashboard if not
    if not validate_genome_build(genome_build):
        delete_selected_genome_build_from_session()
        _log.error("Genome build '{}' is not available in the database. Redirecting to default".format(genome_build))
        return dashboard()
    else:
        # add valid genome_build to the session
        session['genome_build'] = genome_build

    # Retrieve all gene names
    gene_names = GeneRepository.retrieve_all_gene_names_from_file()

    # Render and return the template
    return render_template('dashboard.html', data=map(json.dumps, gene_names), genome_build = genome_build) #@todo move genome build to session handling

@bp.route('/dashboard/<genome_build>/<gene_name>/', methods=['GET'])
def dashboard_gene_name(genome_build, gene_name):
    # We can add the gene name to the session, even if it does not exist in the database
    session['gene_name'] = gene_name

    # check of transcript_ids_for_gene is in the session
    transcript_ids_for_gene = get_transcript_ids_for_gene_from_session()
    query_required = True

    # check if the genome build and gene name match the session variables
    if transcript_ids_for_gene is not None:
        if 'genome_build' in transcript_ids_for_gene.keys() and transcript_ids_for_gene['genome_build'] == genome_build:
            if 'gene_name' in transcript_ids_for_gene.keys() and transcript_ids_for_gene['gene_name'] == gene_name:
                # no need to query the database again if both genome build and gene name match
                query_required = False
    else:
        # ensure a clean state if transcript_ids_for_gene is None
        delete_transcript_ids_for_gene_from_session()

    if query_required:
        # reset the session variables for gene name, geome build, and transcript ids
        session['transcript_ids_for_gene'] = json.loads(get_transcript_ids_for_gene(genome_build, gene_name).get_data().decode('utf-8'))

    # Check if there is a transcript id selected in the session
    transcript_id = get_selected_transcript_id_from_session()
    if transcript_id is not None:
        transcript_id_match = False

        # check if the transcript id is valid for the retrieved transcript ids for gene
        if 'transcript_ids' in transcript_ids_for_gene.keys():
            for transcript in transcript_ids_for_gene['transcript_ids']:
                if transcript['gencode_id'].startswith(transcript_id):
                    # check if the transcript contains a dot, if so, split and check only the main part
                    if not '.' in transcript_id and '.' in transcript['gencode_id']:
                        transcript_id_match = transcript['gencode_id'].split('.')[0] == transcript_id
                    else:
                        transcript_id_match = transcript['gencode_id'] == transcript_id

                    # break the loop if we found a match
                    if transcript_id_match:
                        break

        if not transcript_id_match:
            delete_selected_transcript_id_from_session()
            raise Exception("Transcript id '{}' is not valid. Redirected to gene_name endpoint.".format(transcript_id))

    else:
        # ensure a clean state if transcript_id is None
        delete_selected_transcript_id_from_session()

    # After adding session variables, handle further validation logic in the dashboard_genome_build() method
    return dashboard_genome_build(genome_build)

@bp.route('/dashboard/<genome_build>/<gene_name>/<transcript_id>/', methods=['GET'])
def transcript(genome_build, gene_name, transcript_id):
    # validate transcript id format
    valid_transcript_id_format = False

    # type check the transcript id
    transcript_id_corrected = transcript_id.strip().upper()

    # check if transcript_id_corrected contains ENST, followed by numericals and optionally a dot:
    if transcript_id_corrected.startswith('ENST'):
        # tokenize by 'ENST' and check the rest
        transcript_id_corrected_tail = transcript_id_corrected[4:]  # get the part after 'ENST
        if '.' in transcript_id_corrected_tail:
            main_part, dot_part = transcript_id_corrected_tail.split('.', 1)
            valid_transcript_id_format = main_part.isdigit() and dot_part.isdigit()
        else:
            valid_transcript_id_format = transcript_id_corrected_tail.isdigit()

    if not valid_transcript_id_format:
        delete_selected_transcript_id_from_session()
        return visualization_error("Transcript id '{}' is not valid. Redirected to gene_name endpoint.".format(transcript_id))

    # First make sure the string is upper case and has no leading or trailing spaces
    if transcript_id_corrected != transcript_id:
        delete_selected_transcript_id_from_session()
        return redirect(url_for('web.transcript', genome_build=genome_build, gene_name=gene_name, transcript_id=transcript_id_corrected))

    # Add the transcript id to the session after validation
    session['transcript_id'] = transcript_id

    return dashboard_gene_name(genome_build, gene_name)

@bp.route('/contact', methods=['GET', 'POST'])
def contact():
    if DEFAULT_RECIPIENT is None or MAIL_SERVER is None:
        return render_template('error.html', msg="Mail server is not configured, therefore there can be no mails sent to support."), 500
    
    form = SupportForm()
    if form.validate_on_submit():
        email = form.email.data
        body = form.body.data
        
        msg = Message(subject="Support Request",
                  sender=email,
                  recipients=[DEFAULT_RECIPIENT])
        
        msg.body = body
        
        mail.send(msg)

        return render_template('contact_sent.html')

    return render_template('contact.html', form=form)

@bp.route('/about', methods=['GET'])
def about():
    return render_template('about.html')

@bp.route('/method', methods=['GET'])
def method():
    return render_template('method.html')

@bp.route('/faq', methods=['GET'])
def help_page():
    return render_template('faq.html')

@bp.route('/visualization_error/<genome_build>/<transcript_id>/', methods=['GET'])
def visualization_error(genome_build, transcript_id):
    stacktrace = retrieve_error(transcript_id, genome_build)
    error = "error during visualization generation"
    return render_template('error.html', msg=error, stack_trace=stacktrace)

@bp.before_request
def before_request():
    g.metadome_version = get_version()

@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    _log.error("Unhandled exception:\n{}".format(error))
    _log.error(traceback.format_exc())
    return render_template('error.html', msg=error, stack_trace=traceback.format_exc()), 500

@bp.route('/get_dashboard_session')
def get_dashboard_session():
    result = {}

    # Check if genome_builds are in the session
    result['genome_builds'] = get_genome_builds_from_session()

    # retrieve the selected genome build from the session
    result['genome_build'] = get_selected_genome_build_from_session()
    if result['genome_build'] is None:
        # ensure the session remains clear of other variables
        delete_selected_genome_build_from_session()
        return jsonify(genome_builds = result['genome_builds'])

    # retrieve the selected gene name from the session
    result['gene_name'] = get_selected_gene_name_from_session()
    if result['gene_name'] is None:
        # ensure a partial session reset
        delete_selected_gene_name_from_session()
        # return jsonify(data = result['data'], genome_builds = result['genome_builds'], genome_build = result['genome_build'])
        return jsonify(genome_builds = result['genome_builds'], genome_build = result['genome_build'])

    # Check if transcript ids have already been retrieved for the selected gene
    transcript_ids_for_gene = get_transcript_ids_for_gene_from_session()
    if transcript_ids_for_gene is None:
        # ensure a partial session reset
        delete_transcript_ids_for_gene_from_session()
        return jsonify(genome_builds = result['genome_builds'], genome_build = result['genome_build'], gene_names = result['gene_name'])

    result['transcript_ids'] = transcript_ids_for_gene['transcript_ids']

    # check if a transcript id is selected
    result['transcript_id'] = get_selected_transcript_id_from_session()
    if result['transcript_id'] is None:
        # ensure a partial session reset
        delete_selected_transcript_id_from_session()
        return jsonify(genome_builds = result['genome_builds'], genome_build = result['genome_build'], gene_name = result['gene_name'], transcript_ids_for_gene = transcript_ids_for_gene)

    # if all variables are set, return the session data
    return jsonify(genome_builds = result['genome_builds'], genome_build = result['genome_build'], gene_name = result['gene_name'], transcript_ids_for_gene = transcript_ids_for_gene, transcript_id = result['transcript_id'])

### Helper functions & Session handling
def reset_dashboard_session():
    # Deletes user input and system input from the session
    session.clear()

def delete_selected_transcript_id_from_session():
    """ Delete the selected transcript id from the session."""
    session.pop('transcript_id', None)

def get_selected_transcript_id_from_session():
    """ Retrieve the selected transcript id from the session, if available."""
    return session.get('transcript_id', None)

def delete_transcript_ids_for_gene_from_session():
    """ Delete the transcript ids for a gene from the session."""
    session.pop('transcript_ids_for_gene', None)
    # Also delete the transcript id to ensure a clean state
    session.pop('transcript_id', None)

def get_transcript_ids_for_gene_from_session():
    """ Retrieve the transcript ids for a gene from the session, if available."""
    transcript_ids_for_gene = session.get('transcript_ids_for_gene', None)
    if transcript_ids_for_gene is not None:
        try:
            # check if the transcript_ids_for_gene variable holds the expected keys
            if 'transcript_ids' not in transcript_ids_for_gene.keys():
                raise Exception("Field 'transcript_ids' is missing in 'transcript_ids_for_gene' session variable.")
            if 'gene_name' not in transcript_ids_for_gene.keys():
                raise Exception("Field 'gene_name' is missing in 'transcript_ids_for_gene' session variable.")
            if 'genome_build' not in transcript_ids_for_gene.keys():
                raise Exception("Field 'genome_build' is missing in 'transcript_ids_for_gene' session variable.")
            if 'message' not in transcript_ids_for_gene.keys():
                raise Exception("Field 'message' is missing in 'transcript_ids_for_gene' session variable.")
        except Exception as e:
            raise Exception("Error in `get_transcript_ids_for_gene_from_session` with message: {}, and stack trace: {}".format(e, traceback.format_exc()))

    return transcript_ids_for_gene

def delete_selected_gene_name_from_session():
    """ Delete the selected gene name from the session."""
    session.pop('gene_name', None)
    # Also delete the transcript id and transcript ids for gene to ensure a clean state
    session.pop('transcript_ids_for_gene', None)
    session.pop('transcript_id', None)

def get_selected_gene_name_from_session():
    """ Retrieve the selected gene name from the session, if available."""
    return session.get('gene_name', None)

def delete_selected_genome_build_from_session():
    """ Delete all user input from the session, e.g. gene_name, transcript_id, transcript_ids_for_gene."""
    session.pop('genome_build', None)
    # Also delete the gene name, transcript ids for gene, and transcript id to ensure a clean state
    session.pop('gene_name', None)
    session.pop('transcript_ids_for_gene', None)
    session.pop('transcript_id', None)

def get_selected_genome_build_from_session():
    """ Retrieve the selected genome build from the session, if available."""
    genome_build = session.get('genome_build', None)
    if genome_build is not None:
        # safety check for the genome build
        try:
            genome_url_safe = dashboard_genome_build_safety_check(genome_build)
            if genome_url_safe != genome_build:
                raise Exception("Genome build in session '{}' does not match the available genome builds.".format(genome_build))
        except Exception as e:
            raise Exception("Error in `get_selected_genome_build_from_session` with message: {}, and stack trace: {}".format(e, traceback.format_exc()))

    return genome_build

def get_genome_builds_from_session():
    """ Retrieve the genome builds from the session, if available."""
    genome_builds = session.get('genome_builds', None)
    if genome_builds is None:
        genome_builds = GeneRepository.retrieve_all_genome_builds_from_file()
        session['genome_builds'] = genome_builds

    return genome_builds

def validate_genome_build(genome_build):
    """ Validate the genome build against the available genome builds in the session."""
    genome_url_safe = dashboard_genome_build_safety_check(genome_build)
    return genome_build == genome_url_safe

def dashboard_genome_build_safety_check(genome_build):
    """ Check if the genome build is available in the file extracted from the database."""
    genome_builds = get_genome_builds_from_session()

    # check which genome build is available for genome build
    result = [x for x in genome_builds if x.lower().startswith(genome_build.lower())]

    # redirect to the dashboard with the genome build
    if len(result) != 1:
        raise Exception("No genome build found for "+str(genome_build))

    return result[0]