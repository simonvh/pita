import sys
import atexit
from pita.db_backend import *

def db_session(conn, new=True):
    if not hasattr(db_session, 'session') or not db_session.session:
        engine = create_engine(conn)
        db_session.engine = engine
        if new:
            sys.stderr.write("Dropping all tables")
            Base.metadata.drop_all(db_session.engine)
        sys.stderr.write("Creating tables")
        Base.metadata.create_all(engine)
        Base.metadata.bind = engine
        db_session.session = scoped_session(sessionmaker(bind=engine))
    elif new:
        db_session.session.commit()
        sys.stderr.write("Dropping all tables")
        Base.metadata.drop_all(db_session.engine)
        sys.stderr.write("Creating tables")
        Base.metadata.create_all(db_session.engine)
    
    return db_session.session
