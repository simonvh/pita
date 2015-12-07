import logging
import sys
from pita.config import DEBUG_LEVELS

def setup_logging(basename, debug_level):
    debug_level = debug_level.upper()
    
    if not debug_level in DEBUG_LEVELS:
        sys.stderr.write("Invalid debug level {0}\n".format(debug_level))
        sys.stderr.write("Valid values are {0}\n".format(",".join(DEBUG_LEVELS)))
        sys.exit(1)

    logger = logging.getLogger("pita")
    logger.setLevel(getattr(logging, debug_level))
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(getattr(logging, debug_level))
    fh = logging.FileHandler("{0}.log".format(basename))
    fh.setFormatter(formatter)
    fh.setLevel(getattr(logging, debug_level))
    logger.addHandler(handler)
    logger.addHandler(fh)

    return logger

