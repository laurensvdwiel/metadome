from flask import Blueprint, g, render_template, redirect, url_for, session, jsonify, current_app, request, flash
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_mail import Message
from metadome import get_version
from metadome.controllers.job import retrieve_error
from metadome.domain.repositories import GeneRepository
from metadome.domain.services.mail.mail import mail
from metadome.presentation.api.routes import get_transcript_ids_for_gene
from metadome.presentation.web.forms import SupportForm
from metadome.default_settings import MAIL_SERVER, MAIL_DEFAULT_SENDER, SUPPORT_EMAIL, ERROR_EMAIL_NOTIFICATION_WINDOW, \
    ISSUES_EMAIL
import json
import traceback
import time
import requests
from datetime import datetime

import logging

_log = logging.getLogger(__name__)

bp = Blueprint('web', __name__)

# Initialize limiter at the top of your routes file
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["200 per day", "50 per hour"]
)

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
            _log.error("Transcript id '{}' in dashboard_gene_name() is not valid. Redirected to gene_name endpoint.".format(transcript_id))
            return redirect(url_for('web.dashboard_gene_name', genome_build=genome_build, gene_name=gene_name))

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
        _log.error("Transcript id '{}' in transcript() is not valid. Redirected to gene_name endpoint.".format(transcript_id))
        return redirect(url_for('web.dashboard_gene_name', genome_build=genome_build, gene_name=gene_name))

    # First make sure the string is upper case and has no leading or trailing spaces
    if transcript_id_corrected != transcript_id:
        delete_selected_transcript_id_from_session()
        return redirect(url_for('web.transcript', genome_build=genome_build, gene_name=gene_name, transcript_id=transcript_id_corrected))

    # Add the transcript id to the session after validation
    session['transcript_id'] = transcript_id

    return dashboard_gene_name(genome_build, gene_name)

@bp.route('/dashboard/<genome_build>/<gene_name>/<transcript_id>/p.<int:position>', methods=['GET'], strict_slashes=False)
def transcript_position(genome_build, gene_name, transcript_id, position):
    """Handle position-specific links - validate position against transcript length"""

    # Get transcript info from session (already validated in transcript() route)
    transcript_ids_for_gene = get_transcript_ids_for_gene_from_session()

    # Find the matching transcript to get its length
    valid_position = False
    if transcript_ids_for_gene and 'transcript_ids' in transcript_ids_for_gene:
        for transcript_data in transcript_ids_for_gene['transcript_ids']:
            # Match the transcript ID (handle with/without version)
            transcript_gencode_id = transcript_data['gencode_id']
            if transcript_gencode_id == transcript_id or transcript_gencode_id.split('.')[0] == transcript_id:
                # Validate position is within the actual protein length
                aa_length = transcript_data['aa_length']
                if 1 <= position <= aa_length:
                    valid_position = True
                    session['selected_position'] = position
                else:
                    _log.warning(f"Position {position} out of range for {transcript_id} (length: {aa_length})")
                break

    if not valid_position:
        _log.warning(f"Invalid position {position} for {transcript_id}, redirecting to transcript view")
        return redirect(
            url_for('web.transcript', genome_build=genome_build, gene_name=gene_name, transcript_id=transcript_id))

    # Call the transcript route function directly
    return transcript(genome_build, gene_name, transcript_id)

@bp.route('/get_selected_position')
def get_selected_position():
    """Retrieve the selected position from session"""
    position = session.get('selected_position', None)
    return jsonify(selected_position=position)

@bp.route('/clear_selected_position', methods=['POST'])
def clear_selected_position():
    """Clear the selected position from session after it's been used"""
    session.pop('selected_position', None)
    return jsonify(success=True)

@bp.route('/contact', methods=['GET', 'POST'])
@limiter.limit("5 per minute; 30 per hour")
def contact():
    if MAIL_SERVER is None or MAIL_DEFAULT_SENDER is None or SUPPORT_EMAIL is None:
        return render_template('error.html', msg="Mail server is not configured, therefore there can be no mails sent to support."), 500
    
    form = SupportForm()
    if request.method == 'POST':
        # 1. Check honeypot first
        if request.form.get('website'):
            _log.warning(f"Honeypot triggered from IP: {request.remote_addr}")
            # Fake success to confuse bots
            flash('Thank you for your message! We will get back to you soon.', 'success')
            return redirect(url_for('web.index'))

        # 2. Verify reCAPTCHA
        recaptcha_token = request.form.get('recaptcha_token')
        if current_app.config.get('RECAPTCHA_SECRET_KEY'):
            try:
                verify_response = requests.post(
                    'https://www.google.com/recaptcha/api/siteverify',
                    data={
                        'secret': current_app.config['RECAPTCHA_SECRET_KEY'],
                        'response': recaptcha_token
                    },
                    timeout=5
                )

                result = verify_response.json()
                score = result.get('score', 0)

                if not result.get('success') or score < current_app.config.get('RECAPTCHA_THRESHOLD', 0.5):
                    _log.warning(f"reCAPTCHA failed - IP: {request.remote_addr}, Score: {score}")
                    flash('Security verification failed. Please try again.', 'error')
                    return render_template('contact.html', form=form)

            except Exception as e:
                _log.error(f"reCAPTCHA error: {str(e)}")
                # Continue anyway but log the error

    if form.validate_on_submit():
        email = form.email.data
        body = form.body.data

        # 1. Send to support team
        support_msg = Message(
            subject="Support Request",
            sender=('MetaDome Support', MAIL_DEFAULT_SENDER),
            recipients=[SUPPORT_EMAIL],
            reply_to=email
        )
        support_msg.body = f"""New support request from: {email}

        {body}"""

        # 2. Send confirmation to user with HTML template
        user_msg = Message(
            subject="We received your message - MetaDome Support",
            sender=('MetaDome Support', MAIL_DEFAULT_SENDER),
            recipients=[email],
            reply_to=SUPPORT_EMAIL
        )

        # Render HTML template
        # For emails, we need absolute URLs for images
        logo_url = url_for('static', filename='img/metadome_logo.png', _external=True)

        user_msg.html = render_template(
            'emails/default_support_email.html',
            body=body,
            logo_url=logo_url
        )

        # Send both emails
        try:
            mail.send(support_msg)
            mail.send(user_msg)
        except Exception as e:
            error = f'Failed to send email: {str(e)}'
            render_template('error.html', msg=error, stack_trace=traceback.format_exc())
            return render_template('contact.html', form=form)

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
    stacktrace, error_timestamp = retrieve_error(transcript_id, genome_build)
    error = "error during visualization generation"

    # Handle sending error notification email if the error is recent
    if error_timestamp:
        send_error_notification(transcript_id, genome_build, stacktrace, error_timestamp)

    return render_template('error.html', msg=error, stack_trace=stacktrace)

def send_error_notification(transcript_id, genome_build, error_content, error_timestamp):
    """Send email notification for recent error"""
    if ISSUES_EMAIL is None or MAIL_SERVER is None:
        _log.error("Mail server not configured, cannot send error notification")
        return

    # Create a unique key for this error
    cache_key = f"error_email:{transcript_id}:{genome_build}:{int(error_timestamp)}"

    # Check if error is recent (within last 5 minutes)
    if (time.time() - error_timestamp) >= ERROR_EMAIL_NOTIFICATION_WINDOW:
        _log.error(f"Error '{cache_key}' occurred to long ago for mail to be send not configured, cannot send error notification")
        return

    # Get cache from current app
    cache = current_app.cache

    # Check if we've already sent this exact error
    if cache.get(cache_key) is not None:
        _log.info(f"Already sent email for {cache_key} recently")
        return

    try:
        # Create error notification email
        msg = Message(
            subject=f"MetaDome Error: {transcript_id} ({genome_build})",
            sender=('MetaDome Automated', MAIL_DEFAULT_SENDER),
            recipients=[ISSUES_EMAIL]
        )

        # Format timestamp
        error_time = datetime.fromtimestamp(error_timestamp).strftime('%Y-%m-%d %H:%M:%S UTC')

        # Get logo URL (absolute URL for email)
        logo_url = url_for('static', filename='img/metadome_logo.png', _external=True)

        # Render HTML template
        msg.html = render_template(
            'emails/default_issue_email.html',
            transcript_id=transcript_id,
            genome_build=genome_build,
            error_time=error_time,
            error_content=error_content,
            logo_url=logo_url
        )

        mail.send(msg)

        # Add the cache key to prevent sending duplicate emails for the same error within the notification window.
        # This will help avoid spamming the support team with multiple emails for the same issue.
        cache.set(cache_key, True, timeout=3600)  # Remember for 1 hour

        _log.info(f"Error notification sent for {transcript_id} ({genome_build})")

    except Exception as e:
        _log.error(f"Failed to send error notification: {str(e)}")

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