from flask_sqlalchemy import SQLAlchemy
from sqlalchemy_schemadisplay import create_schema_graph

db = SQLAlchemy()

## Save the schema of the database
def write_db_schema_graph(schema_filename):
    # create the pydot graph object by autoloading all tables via a bound metadata object
    graph = create_schema_graph(db.engine,
        metadata=db.Model.metadata,
        show_datatypes=True,
        show_indexes=True,
        rankdir='LR',  # From left to right (instead of top to bottom)
        concentrate=False,  # Don't try to join the relation lines together
        relation_options={'fontsize': '11.0'}
    )
    graph.set('ranksep', '1.2')  # horizontal gap between table columns
    graph.set('nodesep', '0.6')  # vertical gap between tables
    graph.set('dpi', '150')  # print resolution
    graph.set('mclimit', '4')  # spend extra layout passes on minimizing line crossings
    # tables are transparent plaintext nodes drawn before the edges; make them opaque
    # and draw them last, so relation lines cannot strike through a table
    graph.set('outputorder', 'edgesfirst')
    for node in graph.get_nodes():
        label = node.get('label')
        if label:
            node.set('label', label.replace('<TABLE', '<TABLE BGCOLOR="white"', 1))
            node.set('fontsize', '11.0')  # the library hardcodes 7.0, too small in print
    # replace the labels at the arrow ends by a single mid-edge label; dot reserves
    # layout space for mid-edge labels, so they cannot overlap a table
    for edge in graph.get_edges():
        fk_column = (edge.get('taillabel') or '').lstrip('+ ')
        referenced_column = (edge.get('headlabel') or '').lstrip('+ ')
        edge.set('taillabel', '')
        edge.set('headlabel', '')
        edge.set('label', '%s -> %s' % (fk_column, referenced_column))
        # keep the genes->proteins relation short and straight, so the longer
        # relations arriving at proteins route around it instead of across it
        if edge.get_source() == 'genes' and edge.get_destination() == 'proteins':
            edge.set('weight', '3')
    graph.write_png(schema_filename)  # write out the file
