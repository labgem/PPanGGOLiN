from pathlib import Path
import hashlib
import json

GOLDEN_JSON = Path("testingDataset/expected_info_files/expected_hashes.json")


def compute_sha256(file_path: Path) -> str:
    """Compute SHA256 of a file."""
    sha256 = hashlib.sha256()
    with file_path.open("rb") as f:
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
    """Compare file hash with golden or update it if --update-golden is set."""
    actual_hash = compute_sha256(file_path)
    golden_hashes = load_golden_hashes()
    fname = str(file_path.relative_to(file_path.parents[1]))  # store relative path in JSON

    if update_golden:
        golden_hashes[fname] = actual_hash
        save_golden_hashes(golden_hashes)
        print(f"[update-golden] Updated hash for {fname}: {actual_hash}")
    else:
        expected_hash = golden_hashes.get(fname)
        assert expected_hash is not None, (
            f"No golden hash for {fname}. Run pytest with --update-golden first."
        )
        assert actual_hash == expected_hash, f"Checksum mismatch for {fname}"
