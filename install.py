from metadome.default_settings import RECONSTRUCT_METADOMAINS
from metadome.domain.services.meta_domain_creation import create_metadomains
from metadome.domain.services.database_creation import create_db
from metadome.database import db
from metadome.application import app


db.init_app(app)
with app.app_context():
    # Extensions like Flask-SQLAlchemy now know what the "current" app
    # is while within this block. Therefore, you can now run........
    db.create_all()
    create_db()

    # now create all meta_domains
    create_metadomains(reconstruct=RECONSTRUCT_METADOMAINS)

    # retrieve all gene names and write to disk
    write_all_genes_names_to_disk()
