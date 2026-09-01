import tables
from pathlib import Path
from dataclasses import dataclass
import typing as tp


@dataclass
class TableAttribute:
    """Describe how a table column should be read and converted.

    :param name: Name of the HDF5 column in the table.
    :param predicate: Optional filter used to decide whether a row should be kept.
    :param transform: Optional callable used to transform the raw value before instantiation.
    """

    name: str
    predicate: tp.Callable = lambda x: True
    transform: tp.Callable = lambda x: x


def get_table_attributes(cls):
    """Collect table column metadata annotations from a dataclass.

    :param cls: Dataclass type describing the HDF5 table schema.
    :return: Mapping of field names to their TableAttribute metadata.
    """
    attrs = {}
    for field_name, annotated_type in tp.get_type_hints(
        cls, include_extras=True
    ).items():
        if tp.get_origin(annotated_type) is tp.Annotated:
            _, table_attr = tp.get_args(annotated_type)
            attrs[field_name] = table_attr
    return attrs


class H5Reader:
    """Lightweight HDF5 table reader for PPanGGOLiN-style dataclass mappings."""

    def __init__(self, h5_file: Path) -> None:
        """Initialize the reader for an HDF5 file.

        :param h5_file: Path to the pangenome HDF5 file.
        :raises FileNotFoundError: If the HDF5 file does not exist.
        """
        if not h5_file.exists():
            raise FileNotFoundError(f"The file {h5_file} does not exist.")
        self.h5_file = h5_file
        self.handle: tables.File | None = None

    def open(self) -> None:
        """Open the underlying HDF5 file in read mode."""
        self.handle = tables.open_file(str(self.h5_file), mode="r")

    def close(self) -> None:
        """Close the HDF5 file if it is currently open."""
        if self.handle:
            self.handle.close()
            self.handle = None

    def __enter__(self) -> "H5Reader":
        """Open the reader when used as a context manager.

        :return: The initialized H5Reader instance.
        """
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Close the reader when leaving a context manager block."""
        self.close()

    def _check(self, obj: tp.Any) -> dict[str, TableAttribute]:
        """Validate that a dataclass exposes a readable HDF5 table definition.

        :param obj: Dataclass type describing the table schema.
        :raises ValueError: If the schema is missing required metadata or the table is absent.
        :return: The table name and column metadata mapping.
        """
        attrs = get_table_attributes(obj)
        if not attrs:
            raise ValueError(f"No TableAttributes found in {obj.__name__}.")
        if "_table" not in obj.__dict__:
            raise ValueError(
                f'Missing required TableAttribute "_table" in {obj.__name__}.'
            )

        table_name = obj._table

        if table_name not in self.handle:
            raise ValueError(
                f"Required table '{table_name}' is missing from the pangenome H5 file: {self.h5_file}"
            )

        return table_name, attrs

    def _apply_override(self, attrs, override):
        """Apply row-level overrides to a table attribute mapping.

        :param attrs: Mapping of fields to TableAttribute definitions.
        :param override: Optional per-field overrides with transform and/or predicate values.
        """
        if not override:
            return

        for field_name, attr in attrs.items():
            if field_name in override:
                if "transform" in override[field_name]:
                    attr.transform = override[field_name]["transform"]
                if "predicate" in override[field_name]:
                    attr.predicate = override[field_name]["predicate"]

    def fetch(
        self, obj: tp.Any, override: dict = {}
    ) -> tp.Generator[tp.Any, None, None]:
        """Yield dataclass instances generated from rows in the matching HDF5 table.

        :param obj: Dataclass type used to construct each row object.
        :param override: Optional field overrides for transform and predicate hooks.
        :raises RuntimeError: If the HDF5 file is not open.
        :return: Generator of instantiated dataclass objects.
        """
        if self.handle is None:
            raise RuntimeError(
                "H5 file is not opened. Use 'with' statement or call open() method."
            )

        table_name, attrs = self._check(obj)

        h5_table = self.handle.get_node(table_name)

        self._apply_override(attrs, override)

        for row in h5_table.iterrows():
            args = {}
            for field_name, attr in attrs.items():
                value = row[attr.name]
                if not attr.predicate(value):
                    break
                else:
                    args[field_name] = attr.transform(value)
            else:
                yield obj(**args)

    def fetch_raw(self, obj: tp.Any) -> tp.Generator[tp.Any, None, None]:
        """Yield raw table rows as dataclass instances without applying row filtering.

        :param obj: Dataclass type used to construct each object.
        :raises RuntimeError: If the HDF5 file is not open.
        :return: Generator of dataclass objects with raw column values.
        """
        if self.handle is None:
            raise RuntimeError(
                "H5 file is not opened. Use 'with' statement or call open() method."
            )

        table_name, attrs = self._check(obj)
        h5_table = self.handle.get_node(table_name)

        for row in h5_table.iterrows():
            obj(**{field_name: row[attr.name] for field_name, attr in attrs.items()})

    def fetch_rows(self, obj: tp.Any) -> tp.Generator[tp.Any, None, None]:
        """Yield the raw row objects stored in a table without object conversion.

        :param obj: Dataclass type describing the table schema.
        :raises RuntimeError: If the HDF5 file is not open.
        :return: Generator of raw rows from the selected HDF5 table.
        """
        if self.handle is None:
            raise RuntimeError(
                "H5 file is not opened. Use 'with' statement or call open() method."
            )
        table_name, _ = self._check(obj)
        yield from self.handle.get_node(table_name).iterrows()
