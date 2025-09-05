from pathlib import Path
import hashlib
import json
import gzip

GOLDEN_JSON = Path("testingDataset/expected_info_files/expected_hashes.json")


def compute_sha256(file_path: Path, decompress_gz: bool = False) -> str:
    """Compute sha256 hash of a file. If decompress_gz=True, hash the uncompressed content."""
    sha256 = hashlib.sha256()
    if decompress_gz and file_path.suffix == ".gz":
        with gzip.open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
    else:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
    return sha256.hexdigest()


def load_golden_hashes() -> dict:
    """Load golden hashes from JSON, return empty dict if file does not exist."""
    if GOLDEN_JSON.exists():
        with GOLDEN_JSON.open() as f:
            return json.load(f)
    return {}


def save_golden_hashes(hashes: dict):
    """Save golden hashes to JSON."""
    GOLDEN_JSON.parent.mkdir(parents=True, exist_ok=True)
    with GOLDEN_JSON.open("w") as f:
        json.dump(hashes, f, indent=4)


def assert_or_update_hash(file_path: Path, update_golden: bool):
    """Compare file hash with golden or update it if --update-golden is set.
    If file is gzipped, compute hash on the uncompressed content and drop the `.gz` extension.
    """
    is_gz = file_path.suffix == ".gz"
    actual_hash = compute_sha256(file_path, decompress_gz=is_gz)

    # Normalize filename (strip .gz if needed)
    fname = str(file_path.relative_to(file_path.parents[1]))
    if is_gz:
        fname = fname[:-3]  # remove ".gz"

    golden_hashes = load_golden_hashes()

    if update_golden:
        golden_hashes[fname] = actual_hash
        save_golden_hashes(golden_hashes)
        print(f"[update-golden] Updated hash for {fname}: {actual_hash}")
    else:
        expected_hash = golden_hashes.get(fname)
        assert (
            expected_hash is not None
        ), f"No golden hash for {fname}. Run pytest with --update-golden first."
        assert actual_hash == expected_hash, f"Checksum mismatch for {fname}"
