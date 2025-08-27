import logging
from unittest.mock import patch
from ppanggolin.main import main


# create a module-level logger
logger = logging.getLogger(__name__)

def run_ppanggolin_command(cmd):
    logger.info(cmd)
    with patch("sys.argv", cmd.split()):
        main()
