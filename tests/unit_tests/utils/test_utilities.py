#! /usr/bin/env python3

import pytest
from pathlib import Path
import gzip
import bz2
import zipfile
from typing import Generator

from ppanggolin.utils import (
    is_compressed,
    read_compressed_or_not,
    write_compressed_or_not,
    has_non_ascii,
    replace_non_ascii,
)


class TestCompressed:
    """
    Test cases for the is_compressed function.
    """

    @pytest.fixture
    def plain_file(self, tmp_path: Path) -> Generator[Path, None, None]:
        """
        Creates a temporary plain text file for testing.
        """
        file_path = tmp_path / "test.txt"
        with open(file_path, "wb") as f:
            f.write(b"Test data")
        yield file_path

    @pytest.fixture
    def gzip_file(self, tmp_path: Path) -> Generator[Path, None, None]:
        """
        Creates a temporary gzip file for testing.
        """
        file_path = tmp_path / "test.gz"
        with gzip.open(file_path, "wb") as f:
            f.write(b"Test data")
        yield file_path

    @pytest.fixture
    def bz2_file(self, tmp_path: Path) -> Generator[Path, None, None]:
        """
        Creates a temporary bz2 file for testing.
        """
        file_path = tmp_path / "test.bz2"
        with bz2.open(file_path, "wb") as f:
            f.write(b"Test data")
        yield file_path

    @pytest.fixture
    def zip_file(self, tmp_path: Path) -> Generator[Path, None, None]:
        """
        Creates a temporary zip file for testing.
        """
        file_path = tmp_path / "test.zip"
        with zipfile.ZipFile(file_path, "w") as z:
            z.writestr("test.txt", "Test data")
        yield file_path


class TestIsCompressed(TestCompressed):
    def test_is_compressed_with_plain_file(self, plain_file: Path) -> None:
        """Test is_compressed function with a plain text file."""
        assert is_compressed(plain_file) == (False, None)

    def test_is_compressed_with_gzip_file(self, gzip_file: Path) -> None:
        """
        Test is_compressed function with a gzip file.
        """
        assert is_compressed(gzip_file) == (True, "gzip")

    def test_is_compressed_with_bz2_file(self, bz2_file: Path) -> None:
        """
        Test is_compressed function with a bz2 file.
        """
        assert is_compressed(bz2_file) == (True, "bz2")

    def test_is_compressed_with_zip_file(self, zip_file: Path) -> None:
        """
        Test is_compressed function with a zip file.
        """
        assert is_compressed(zip_file) == (True, "zip")

    def test_is_compressed_with_unsupported_type(self) -> None:
        """
        Test is_compressed function with an unsupported file type.
        """
        with pytest.raises(TypeError):
            is_compressed(123)


class TestReadCompressedOrNot(TestCompressed):
    """
    Test cases for the read_compressed_or_not function.
    """

    def test_read_compressed_gzip(self, gzip_file: Path) -> None:
        """
        Test read_compressed_or_not function with a gzip file.
        """
        with read_compressed_or_not(gzip_file) as f:
            assert f.read() == "Test data"

    def test_read_compressed_bz2(self, bz2_file: Path) -> None:
        """
        Test read_compressed_or_not function with a bz2 file.
        """
        with read_compressed_or_not(bz2_file) as f:
            assert f.read() == "Test data"

    def test_read_compressed_zip(self, zip_file: Path) -> None:
        """
        Test read_compressed_or_not function with a zip file.
        """
        with read_compressed_or_not(zip_file) as f:
            assert f.read() == "Test data"

    def test_read_plain_file(self, plain_file: Path) -> None:
        """
        Test read_compressed_or_not function with a plain text file.
        """
        with read_compressed_or_not(plain_file) as f:
            assert f.read() == "Test data"

    def test_read_unsupported_type(self) -> None:
        """
        Test read_compressed_or_not function with an unsupported file type.
        """
        with pytest.raises(TypeError):
            read_compressed_or_not(123)


class TestWriteCompressedOrNot:
    """
    Test cases for the write_compressed_or_not function.
    """

    @pytest.fixture
    def plain_file_path(self, tmp_path: Path) -> Path:
        """
        Provides a temporary file path for testing.
        """
        return tmp_path / "test.txt"

    def test_write_compressed(self, plain_file_path: Path) -> None:
        """
        Test write_compressed_or_not function with compression enabled.
        """
        with write_compressed_or_not(plain_file_path, compress=True) as f:
            f.write("Test data")
        with gzip.open(plain_file_path.with_suffix(".txt.gz"), "rt") as f:
            assert f.read() == "Test data"

    def test_write_uncompressed(self, plain_file_path: Path) -> None:
        """
        Test write_compressed_or_not function without compression.
        """
        with write_compressed_or_not(plain_file_path, compress=False) as f:
            f.write("Test data")
        with open(plain_file_path, "r") as f:
            assert f.read() == "Test data"


# Test cases for has_non_ascii
@pytest.mark.parametrize(
    "input_string, expected",
    [
        ("Escherichia_coli", False),  # All ASCII characters
        ("Escherichia_colí", True),  # Contains non-ASCII character 'í'
        ("simple_string", False),  # Simple ASCII string
        ("Ωmega", True),  # Contains non-ASCII character 'Ω'
        ("", False),  # Empty string should return False
    ],
)
def test_has_non_ascii(input_string, expected):
    assert has_non_ascii(input_string) == expected


# Test cases for replace_non_ascii
@pytest.mark.parametrize(
    "input_string, replacement, expected",
    [
        (
            "Escherichia_coli",
            "_",
            "Escherichia_coli",
        ),  # All ASCII characters, no replacement needed
        ("Escherichia_colí", "_", "Escherichia_col_"),  # Replace 'í' with '_'
        ("Ωmega", "-", "-mega"),  # Replace 'Ω' with '-'
        ("Escherichia_Ωcoli", "X", "Escherichia_Xcoli"),  # Replace 'Ω' with 'X'
        ("", "_", ""),  # Empty string, no replacement
    ],
)
def test_replace_non_ascii(input_string, replacement, expected):
    assert replace_non_ascii(input_string, replacement) == expected
