from calendar import c
from unittest.mock import patch
from ppanggolin.main import main


def run_ppanggolin_command(cmd):
    print(cmd)
    with patch("sys.argv", cmd.split()):
        main()
