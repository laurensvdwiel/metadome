import logging

from flask import Flask
import flask_jsglue
from metadome.domain.services.mail.mail import mail
from metadome.default_settings import MAIL_SERVER
from celery import Celery

_log = logging.getLogger(__name__)

def create_app(settings=None):
    _log.info("Creating flask app")

    app = Flask(__name__, static_folder='presentation/web/static', static_url_path='/metadome/static',
               template_folder='presentation/web/templates')
    app.config.from_object('metadome.default_settings')
    if settings:
        app.config.update(settings)

    # Ignore Flask's built-in logging
    # app.logger is accessed here so Flask tries to create it
    app.logger_name = "nowhere"
    app.logger

    # Configure logging.
    #
    # It is somewhat dubious to get _log from the root package, but I can't see
    # a better way. Having the email handler configured at the root means all
    # child loggers inherit it.
    from metadome import _log as metadom_logger

    # Only log to the console during development and production, but not during
    # testing.
    if app.testing:
        metadom_logger.setLevel(logging.DEBUG)
    else:
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        metadom_logger.addHandler(ch)

        if app.debug:
            metadom_logger.setLevel(logging.DEBUG)
        else:
            metadom_logger.setLevel(logging.INFO)


    # Add mail instance
    if not MAIL_SERVER is None:
        mail.init_app(app)
    else:
        _log.warn('MAILSERVER is not set in default_settings.py')

    # Initialize extensions
    from metadome import toolbar
    toolbar.init_app(app)


    # Specify the Blueprints
    from metadome.presentation.web.routes import bp as web_bp
    from metadome.presentation.api.routes import bp as api_bp

    # Register the Blueprints
    app.register_blueprint(api_bp, url_prefix='/metadome/api')
    app.register_blueprint(web_bp, url_prefix='/metadome')

    # Database
    from metadome.database import db
    db.init_app(app)
    with app.app_context():
        db.create_all()

    from flask_jsglue import JSGlue
    flask_jsglue.JSGLUE_JS_PATH = '/metadome/static/js/JSGlue.js'
    jsglue = JSGlue(app)

    return app

def make_celery(app):
    _log.info("Creating celery app")

    app = app or create_app()

    celery = Celery(__name__)
    celery.conf.update(
        broker_url=app.config['CELERY_BROKER_URL'],
        result_backend=app.config['CELERY_RESULT_BACKEND'],
        task_track_started=app.config['CELERY_TRACK_STARTED'],
        task_serializer=app.config['CELERY_TASK_SERIALIZER'],
        result_serializer=app.config['CELERY_RESULT_SERIALIZER'],
        accept_content=app.config['CELERY_ACCEPT_CONTENT'])
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask

    # needed here to register celery tasks, this should not be anywhere else
    import metadome.tasks

    return celery
