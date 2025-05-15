from flask import Blueprint, g, render_template, redirect, url_for
from flask_mail import Message
import json
import traceback
from metadome import get_version
from metadome.controllers.job import retrieve_error
from metadome.domain.repositories import GeneRepository
from metadome.domain.services.mail.mail import mail
from metadome.presentation.web.forms import SupportForm
from metadome.default_settings import DEFAULT_RECIPIENT, MAIL_SERVER

import logging

_log = logging.getLogger(__name__)

bp = Blueprint('web', __name__)

@bp.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@bp.route('/dashboard/transcript/', methods=['GET'])
def transcript_empty():
    return redirect(url_for('web.dashboard'))

@bp.route('/dashboard/transcript/<transcript_id>/', methods=['GET'])
def transcript(transcript_id):
    return render_template('dashboard.html', transcript_id=transcript_id)

@bp.route('/dashboard', methods=['GET'])
def dashboard():
    gene_names = GeneRepository.retrieve_all_gene_names_from_file()
    return render_template('dashboard.html', data=map(json.dumps, gene_names))

@bp.route('/dashboard_js')
def dashboard_js():
    # Renders the javascript used on the index page
    return render_template('/js/dashboard.js')

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

@bp.route('/visualization_error/<transcript_id>/', methods=['GET'])
def visualization_error(transcript_id):
    stacktrace = retrieve_error(transcript_id)
    error = "error during visualization generation"
    return render_template('error.html', msg=error, stack_trace=stacktrace)


@bp.before_request
def before_request():
    g.metadom_version = get_version()


@bp.errorhandler(Exception)
def exception_error_handler(error):  # pragma: no cover
    _log.error("Unhandled exception:\n{}".format(error))
    _log.error(traceback.format_exc())
    return render_template('error.html', msg=error, stack_trace=traceback.format_exc()), 500
