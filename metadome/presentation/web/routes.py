from flask import Blueprint, g, render_template, redirect, url_for, current_app, request, flash, make_response
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_mail import Message
from metadome import get_version
from metadome.controllers.job import retrieve_error
from metadome.domain.repositories import GeneRepository, MappingRepository
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
import re

_log = logging.getLogger(__name__)

bp = Blueprint('web', __name__)

# Initialize limiter at the top of your routes file
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["200 per day", "50 per hour"]
)

@bp.app_context_processor
def inject_navbar_search_context():
    genome_builds = GeneRepository.retrieve_all_genome_builds_from_db()
    gene_names = GeneRepository.retrieve_all_gene_names_from_db()
    gene_data = list(map(json.dumps, gene_names))

    return {
        "navbar_genome_builds": genome_builds,
        "navbar_selected_genome_build": genome_build_safety_check('GRCh38'),
        "navbar_gene_data": gene_data
    }

@bp.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@bp.route('/frontpage/search', methods=['GET'])
def frontpage_search():
    genome_build = genome_build_safety_check(request.args.get('genome_build', 'GRCh38'))
    query = (request.args.get('query', '') or '').strip()

    chromosome, position = parse_position_query(query)
    if chromosome and position is not None:
        return redirect(url_for(
            'web.position_lookup',
            genome_build=genome_build,
            chromosome=chromosome,
            position=position
        ))

    gene_name = parse_gene_query(query)
    if gene_name:
        return redirect(url_for(
            'web.dashboard_gene_name',
            genome_build=genome_build,
            gene_name=gene_name
        ))

    flash("Please enter a gene name or genomic position.", "warning")
    return redirect(url_for('web.index'))

@bp.route('/dashboard_js')
def dashboard_js():
    # Renders the javascript used on the index page
    response = make_response(render_template('/js/dashboard.js'))
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response

@bp.route('/dashboard/', methods=['GET'])
def dashboard():
    # Redirect to the GRCh38 dashboard by default
    return dashboard_grch38()

@bp.route('/dashboard/tour/', methods=['GET'])
def dashboard_tour():
    return render_template('dashboard.html')

@bp.route('/dashboard/tour/<path:extra_path>', methods=['GET'], strict_slashes=False)
def dashboard_tour_redirect(extra_path):
    return redirect(url_for('web.dashboard_tour'))

@bp.route('/dashboard_grch37/', methods=['GET'])
def dashboard_grch37():
    genome_build_url_safe = genome_build_safety_check('GRCh37')

    return redirect(url_for('web.dashboard_genome_build', genome_build=genome_build_url_safe))

@bp.route('/dashboard_grch38/', methods=['GET'])
def dashboard_grch38():
    genome_build_url_safe = genome_build_safety_check('GRCh38')

    return redirect(url_for('web.dashboard_genome_build', genome_build=genome_build_url_safe))

def _render_dashboard(genome_build, gene_name=None, transcript_id=None):
    """Render the dashboard from URL state only."""
    if not validate_genome_build(genome_build):
        _log.error("Genome build '%s' is not available in the database. Redirecting to default", genome_build)
        return dashboard()

    gene_names = GeneRepository.retrieve_all_gene_names_from_db()

    response = make_response(render_template(
        'dashboard.html',
        data=map(json.dumps, gene_names),
        genome_build=genome_build,
        gene_name=gene_name or "",
        transcript_id=transcript_id or ""
    ))
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response

@bp.route('/dashboard/<genome_build>/', methods=['GET'])
def dashboard_genome_build(genome_build):
    return _render_dashboard(genome_build)

@bp.route('/dashboard/<genome_build>/<gene_name>/', methods=['GET'])
def dashboard_gene_name(genome_build, gene_name):
    if not validate_genome_build(genome_build):
        _log.error("Genome build '%s' is not available in the database. Redirecting to default", genome_build)
        return dashboard()

    transcript_ids_for_gene = json.loads(
        get_transcript_ids_for_gene(genome_build, gene_name).get_data().decode('utf-8')
    )

    if not transcript_ids_for_gene or 'transcript_ids' not in transcript_ids_for_gene:
        _log.error("Failed to retrieve transcripts for gene '%s' on build '%s'", gene_name, genome_build)
        return redirect(url_for('web.dashboard_genome_build', genome_build=genome_build))

    return _render_dashboard(genome_build, gene_name=gene_name)

@bp.route('/dashboard/<genome_build>/<gene_name>/<transcript_id>/', methods=['GET'])
def transcript(genome_build, gene_name, transcript_id):
    transcript_context = _get_valid_transcript_context(genome_build, gene_name, transcript_id)

    if not transcript_context['format_valid']:
        _log.error("Transcript id '%s' in transcript() is not valid. Redirected to gene_name endpoint.", transcript_id)
        return _redirect_to_gene_dashboard(genome_build, gene_name)

    if transcript_context['normalized_transcript_id'] != transcript_id:
        return redirect(url_for(
            'web.transcript',
            genome_build=genome_build,
            gene_name=gene_name,
            transcript_id=transcript_context['normalized_transcript_id']
        ))

    if transcript_context['transcript_match'] is None:
        _log.error(
            "Transcript id '%s' in transcript() is not valid for gene '%s'. Redirected to gene_name endpoint.",
            transcript_context['normalized_transcript_id'],
            gene_name
        )
        return _redirect_to_gene_dashboard(genome_build, gene_name)

    if not transcript_context['has_protein_data']:
        _log.warning(
            "Transcript id '%s' for gene '%s' has no protein data. Redirected to gene_name endpoint.",
            transcript_context['normalized_transcript_id'],
            gene_name
        )
        return _redirect_to_gene_dashboard(genome_build, gene_name)

    if _dashboard_query_state_is_invalid(request.args, transcript_context['aa_length']):
        _log.warning(
            "Invalid dashboard query state for transcript '%s'. Redirecting to clean transcript URL.",
            transcript_context['normalized_transcript_id']
        )
        return _redirect_to_clean_transcript(
            genome_build,
            gene_name,
            transcript_context['normalized_transcript_id']
        )

    return _render_dashboard(
        genome_build,
        gene_name=gene_name,
        transcript_id=transcript_context['normalized_transcript_id']
    )

@bp.route('/dashboard/<genome_build>/<gene_name>/<transcript_id>/p.<int:position>', methods=['GET'], strict_slashes=False)
def transcript_position(genome_build, gene_name, transcript_id, position):
    """Legacy position-specific links redirect to the transcript route with query-string state."""

    transcript_context = _get_valid_transcript_context(genome_build, gene_name, transcript_id)

    if transcript_context['transcript_match'] is None or not transcript_context['has_protein_data']:
        _log.warning("Invalid transcript '%s' for legacy position link, redirecting to transcript view", transcript_id)
        return _redirect_to_clean_transcript(
            genome_build,
            gene_name,
            transcript_context['normalized_transcript_id']
        )

    if transcript_context['aa_length'] is None or not (1 <= position <= transcript_context['aa_length']):
        _log.warning(
            "Position %s out of range for %s (length: %s)",
            position,
            transcript_context['normalized_transcript_id'],
            transcript_context['aa_length']
        )
        return _redirect_to_clean_transcript(
            genome_build,
            gene_name,
            transcript_context['normalized_transcript_id']
        )

    return redirect(url_for(
        'web.transcript',
        genome_build=genome_build,
        gene_name=gene_name,
        transcript_id=transcript_context['normalized_transcript_id'],
        selected_positions=str(position),
        positional_information=position
    ))

@bp.route('/position/', methods=['GET'])
def position_lookup_home():
    genome_builds = GeneRepository.retrieve_all_genome_builds_from_db()
    return render_template(
        'position.html',
        genome_builds=genome_builds,
        selected_genome_build=genome_build_safety_check('GRCh38'),
        selected_query="",
        selected_chromosome="",
        selected_position="",
        results=[]
    )

@bp.route('/position/search', methods=['GET'])
def position_lookup_search():
    genome_build = genome_build_safety_check(request.args.get('genome_build'))
    query = request.args.get('query', '')
    chromosome, position = parse_position_query(query)

    if not chromosome or position is None:
        return redirect(url_for('web.position_lookup_home'))

    return redirect(url_for(
        'web.position_lookup',
        genome_build=genome_build,
        chromosome=chromosome,
        position=position
    ))

@bp.route('/position/<genome_build>/<chromosome>/<int:position>/', methods=['GET'])
def position_lookup(genome_build, chromosome, position):
    selected_genome_build = genome_build_safety_check(genome_build)
    normalized_chr = MappingRepository._normalize_chromosome_input(chromosome)

    raw_results = MappingRepository.lookup_position_results(
        chromosome=normalized_chr,
        position=position,
        genome_build=selected_genome_build
    )

    results = []
    for r in raw_results:
        metadome_link = url_for(
            'web.transcript_position',
            genome_build=r["genome_build"],
            gene_name=r["gene_name"],
            transcript_id=r["gencode_transcript"],
            position=r["protein_position"]
        )
        item = dict(r)
        item["metadome_link"] = metadome_link
        item["metadome_label"] = "Analyse"
        results.append(item)

    genome_builds = GeneRepository.retrieve_all_genome_builds_from_db()
    return render_template(
        'position.html',
        genome_builds=genome_builds,
        selected_genome_build=selected_genome_build,
        selected_query=f"{normalized_chr}:{position}",
        selected_chromosome=normalized_chr,
        selected_position=position,
        results=results
    )

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

@bp.route('/privacy', methods=['GET'])
def privacy():
    return render_template('privacy.html')


@bp.route('/visualization_error/<genome_build>/<transcript_id>/', methods=['GET'])
def visualization_error(genome_build, transcript_id):
    stacktrace, error_timestamp = retrieve_error(transcript_id, genome_build)
    error = "error during visualization generation"
    previous_url = request.referrer

    # Handle sending error notification email if the error is recent
    if error_timestamp:
        send_error_notification(
            transcript_id,
            genome_build,
            stacktrace,
            error_timestamp,
            previous_url=previous_url
        )

    return render_template('error.html', msg=error, stack_trace=stacktrace)

def send_error_notification(transcript_id, genome_build, error_content, error_timestamp, previous_url=None):
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
            previous_url=previous_url,
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

def _normalize_transcript_id(transcript_id):
    return (transcript_id or "").strip().upper()

def _is_valid_transcript_id_format(transcript_id):
    if not transcript_id.startswith('ENST'):
        return False

    transcript_id_tail = transcript_id[4:]
    if '.' in transcript_id_tail:
        main_part, dot_part = transcript_id_tail.split('.', 1)
        return main_part.isdigit() and dot_part.isdigit()

    return transcript_id_tail.isdigit()

def _redirect_to_gene_dashboard(genome_build, gene_name):
    return redirect(url_for('web.dashboard_gene_name', genome_build=genome_build, gene_name=gene_name))

def _redirect_to_clean_transcript(genome_build, gene_name, transcript_id):
    return redirect(url_for(
        'web.transcript',
        genome_build=genome_build,
        gene_name=gene_name,
        transcript_id=transcript_id
    ))

def _get_valid_transcript_context(genome_build, gene_name, transcript_id):
    normalized_transcript_id = _normalize_transcript_id(transcript_id)
    format_valid = _is_valid_transcript_id_format(normalized_transcript_id)

    transcript_match = None
    if format_valid:
        transcript_match = _get_transcript_match(genome_build, gene_name, normalized_transcript_id)

    has_protein_data = bool(transcript_match and transcript_match.get('has_protein_data', False))
    aa_length = transcript_match.get('aa_length') if transcript_match else None

    return {
        'normalized_transcript_id': normalized_transcript_id,
        'format_valid': format_valid,
        'transcript_match': transcript_match,
        'has_protein_data': has_protein_data,
        'aa_length': aa_length
    }

def _get_transcript_match(genome_build, gene_name, transcript_id):
    transcript_ids_for_gene = json.loads(
        get_transcript_ids_for_gene(genome_build, gene_name).get_data().decode('utf-8')
    )

    if 'transcript_ids' not in transcript_ids_for_gene:
        return None

    for transcript_data in transcript_ids_for_gene['transcript_ids']:
        gencode_id = transcript_data['gencode_id']
        if gencode_id == transcript_id or gencode_id.split('.')[0] == transcript_id:
            return transcript_data

    return None

def _parse_selected_positions_query(raw_selected_positions, aa_length):
    if raw_selected_positions is None or raw_selected_positions == "":
        return []

    if aa_length is None:
        return None

    if not re.fullmatch(r"\d+(,\d+)*", raw_selected_positions):
        return None

    positions = [int(value) for value in raw_selected_positions.split(',')]

    if len(set(positions)) != len(positions):
        return None

    if positions != sorted(positions):
        return None

    if any(position < 1 or position > aa_length for position in positions):
        return None

    return positions

def _parse_positional_information_query(raw_positional_information, aa_length):
    if raw_positional_information is None or raw_positional_information == "":
        return None

    if aa_length is None:
        return False

    if not str(raw_positional_information).isdigit():
        return False

    position = int(raw_positional_information)
    if not (1 <= position <= aa_length):
        return False

    return position

def _dashboard_query_state_is_invalid(args, aa_length):
    selected_positions = _parse_selected_positions_query(args.get('selected_positions'), aa_length)
    positional_information = _parse_positional_information_query(args.get('positional_information'), aa_length)

    if selected_positions is None:
        return True

    if positional_information is False:
        return True

    allowed_keys = {'selected_positions', 'positional_information'}
    query_keys = set(args.keys())

    if not query_keys.issubset(allowed_keys):
        return True

    return False

def validate_genome_build(genome_build):
    """ Validate the genome build against the available genome builds in the db."""
    genome_url_safe = genome_build_safety_check(genome_build)
    return genome_build == genome_url_safe

def genome_build_safety_check(genome_build):
    """ Check if the genome build is consistent with the database."""
    genome_builds = GeneRepository.retrieve_all_genome_builds_from_db()

    # check which genome build is available for genome build
    result = [x for x in genome_builds if x.lower().startswith(genome_build.lower())]

    # redirect to the dashboard with the genome build
    if len(result) != 1:
        raise Exception("No genome build found for "+str(genome_build))

    return result[0]

def parse_gene_query(query: str) -> str | None:
    raw = (query or "").strip()
    if not raw:
        return None

    if ":" in raw or any(ch.isspace() for ch in raw):
        return None

    return raw

def parse_position_query(query: str) -> tuple[str | None, int | None]:
    raw = (query or "").strip()
    if not raw:
        return None, None

    match = re.match(r"^\s*(chr[\w]+|[\w]+)\s*[:\s]\s*(\d+)\s*$", raw, re.IGNORECASE)
    if not match:
        return None, None

    chromosome = MappingRepository._normalize_chromosome_input(match.group(1))
    position = int(match.group(2)) if match.group(2).isdigit() else None
    return chromosome, position