from pathlib import Path
import shutil

import io
from contextlib import redirect_stdout
from unittest.mock import patch

from ppanggolin.main import main

GOLDEN_DIR = Path("testingDataset/expected_info_files")


def assert_or_update_file(file_path: Path, update_golden: bool):
    """Compare file content with golden file or update it if --update-golden is set."""
    golden_file = GOLDEN_DIR / file_path.name

    if update_golden:
        # Copy current file to golden
        shutil.copy(file_path, golden_file)
        print(f"[update-golden] Updated golden file for {file_path.name}")
    else:
        assert (
            golden_file.exists()
        ), f"No golden file for {file_path.name}. Run pytest with --update-golden first."
        content_actual = file_path.read_text()
        content_expected = golden_file.read_text()
        print(content_actual)
        assert content_actual == content_expected, (
            f"Content mismatch for {file_path.name}. "
            f"Use --update-golden to update the reference."
        )


def check_pangenome_info(pangenome, info_file, update_golden):
    cmd = f"ppanggolin info --pangenome {pangenome}"

    # Capture stdout
    f = io.StringIO()
    with patch("sys.argv", cmd.split()):
        with redirect_stdout(f):
            main()

    # Write captured output to file
    info_file.write_text(f.getvalue())

    assert info_file.stat().st_size > 0
    print(info_file)
    assert_or_update_file(info_file, update_golden)
