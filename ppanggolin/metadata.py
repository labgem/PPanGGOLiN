#!/usr/bin/env python3

# default libraries
import logging
from typing import Generator, List, Tuple, Union, Any, Dict
from collections import defaultdict

# installed libraries
from pandas import isna


class Metadata:
    """The Metadata class represents a metadata link to genes, gene families, organisms, regions, spot or modules.

    Methods:
        - number_of_attribute: Returns the number of attributes in the Metadata object.
        - get: Returns the value of a specific attribute, or None if the attribute does not exist.
        - fields: Returns a list of all the attributes in the Metadata object.

    Fields:
        - source: A string representing the source of the metadata.
        - kwargs: A dictionary of attributes and values representing the metadata. The attributes can be any string, and the values can be any type except None or NaN.
    """

    def __init__(self, source: str, **kwargs):
        """Constructor Method

        :param source: A string representing the source of the metadata.
        :param kwargs: A dictionary of attributes and values representing the metadata. The attributes can be any string, and the values can be any type except None or NaN.

        :raises TypeError: Source name is not a string
        :raises Exception: Source name is empty
        :raises Exception: Metadata is empty
        """
        if not isinstance(source, str):
            raise TypeError(
                f"Metadata source name must be a string. Given type {type(source)}"
            )
        if source == "":
            raise ValueError("Metadata source name should not be empty.")
        self.source = source
        if len(kwargs) == 0:
            raise Exception(f"No metadata given for source: {source}")
        for attr, value in kwargs.items():
            if isinstance(value, list):
                value = self._join_list(value)
            if value is not None and not isna(value):
                setattr(self, attr, value)

    def __repr__(self):
        return f"Metadata source: {self.source}, #attr: {len(self)}"

    def __len__(self) -> int:
        """Get the number of attribute links to the metadata object

        :return: Number of fields (attribute) of the metadata
        """
        return len(self.__dict__) - 1

    def __getattr__(self, attr: str) -> Any:
        """Get the value corresponding to the given attribute

        :return: Value of the attribute

        :raises AttributeError: The attribute does not exist in the metadata
        """
        if attr not in self.__dict__:
            raise AttributeError(f"{attr} is not an attribute of metadata")
        return self[attr]

    @property
    def fields(self) -> List[str]:
        """Get all the field of the metadata

        :return: List of the field in the metadata
        """
        fields = list(self.__dict__)
        fields.remove("source")
        return fields

    def to_dict(self) -> Dict[str, Any]:
        """
        Get metadata in dict format.
        """

        return self.__dict__.copy()

    @staticmethod
    def _join_list(attr_list: Union[str, List[str]]):
        return ",".join(attr_list)


class MetaFeatures:
    """
    The MetaFeatures class provides methods to access and manipulate metadata in all ppanggolin classes.

    Methods
    metadata: Generate all metadata from all sources.
    sources: Generate all metadata sources.
    get_metadata: Get metadata based on attribute values.
    max_metadata_by_source: Gets the source with the maximum number of metadata and the corresponding count.
    """

    def __init__(self):
        """Constructor method"""
        self._metadata_getter = defaultdict(dict)

    @property
    def number_of_metadata(self) -> int:
        """Get the number of metadata associated to feature"""
        return sum(len(meta_dict) for meta_dict in self._metadata_getter.values())

    @property
    def metadata(self) -> Generator[Metadata, None, None]:
        """Generate metadata in gene families

        :return: Metadata from all sources
        """

        for meta_dict in self._metadata_getter.values():
            yield from meta_dict.values()

    @property
    def sources(self) -> Generator[str, None, None]:
        """Get all metadata source in gene family

        :return: Metadata source
        """
        yield from self._metadata_getter.keys()

    def formatted_metadata_dict(self) -> Dict[str, List[str]]:
        """
        Format metadata by combining source and field values.

        Given an object with metadata, this function creates a new dictionary where the keys
        are formatted as 'source_field'.

        :return: A dictionary with formatted metadata.
        """

        source_field_2_values = defaultdict(list)
        for metadata in self.metadata:
            for field in metadata.fields:
                value = str(getattr(metadata, field))

                source_field_2_values[f"{metadata.source}_{field}"].append(value)

        return source_field_2_values

    def formatted_metadata_dict_to_string(self, separator: str = "|") -> Dict[str, str]:
        """
        Format metadata by combining source and field values.

        Given an object with metadata, this function creates a new dictionary where the keys
        are formatted as 'source_field'. In some cases, it is possible to have multiple values for the same field,
        in this situation, values are concatenated with the specified separator.

        :param separator: The separator used to join multiple values for the same field (default is '|').
        :return: A dictionary with formatted metadata.
        """

        source_field_2_concat_values = {}

        source_field_2_values = self.formatted_metadata_dict()

        for source_field, values in source_field_2_values.items():
            if len(values) > 1:
                for value in values:
                    if separator in value:
                        raise ValueError(
                            f"Metadata {source_field}={value} associated to {self} "
                            f"contains in its value the separator character '{separator}'. "
                            "Please change separator in order to be able to write the metadata."
                        )
            source_field_2_concat_values[source_field] = separator.join(values)

        return source_field_2_concat_values

    def add_metadata(self, metadata: Metadata, metadata_id: int = None) -> None:
        """
        Add metadata to metadata getter

        :param metadata: metadata value to add for the source
        :param metadata_id: metadata identifier

        :raises AssertionError: Source or metadata is not with the correct type
        """
        assert isinstance(
            metadata, Metadata
        ), f"Metadata is not with type Metadata but with {type(metadata)}"

        # Metadata_id should not already exist because the metadata are added from scratch to a new source,
        # or they are ridden
        if metadata_id is None:
            if self.has_source(metadata.source):
                metadata_id = max(self._metadata_getter[metadata.source].keys()) + 1
            else:
                # Set first as 1 for PANORAMA
                metadata_id = 1

        try:
            self.get_metadata(metadata.source, metadata_id)
        except KeyError:
            self._metadata_getter[metadata.source][metadata_id] = metadata
        else:
            raise KeyError(
                f"A metadata with ID {metadata_id} already exist "
                f"for source {metadata.source} in {str(self)}"
            )

    def get_metadata(self, source: str, metadata_id: int = None) -> Metadata:
        """Get metadata from metadata getter by its source and identifier

        :param source: source of the metadata
        :param metadata_id: metadata identifier

        :raises KeyError: No metadata with ID or source is found
        """
        try:
            metadata = self._metadata_getter[source][metadata_id]
        except KeyError:
            raise KeyError(
                f"No metadata exist with ID {metadata_id}"
                f"for source {source} in {str(self)}"
            )
        else:
            return metadata

    def get_metadata_by_source(self, source: str) -> Union[Dict[int, Metadata], None]:
        """Get all the metadata feature corresponding to the source

        :param source: Name of the source to get

        :return: List of metadata corresponding to the source

        :raises AssertionError: Source is not with the correct type
        """
        assert isinstance(
            source, str
        ), f"Source is not a string but with {type(source)}"
        return self._metadata_getter.get(
            source
        )  # if source in _metadata_getter return value else None

    def get_metadata_by_attribute(self, **kwargs) -> Generator[Metadata, None, None]:
        """Get metadata by one or more attribute

        :return: Metadata searched
        """
        for metadata in self.metadata:
            for attr, value in kwargs.items():
                if hasattr(metadata, attr):
                    # BUG If value is a list, the join block detection.
                    # It would be better to keep a list and change in writing and reading metadata to join the list
                    if (
                        getattr(metadata, attr, None) in value
                        or getattr(metadata, attr, None) == value
                    ):
                        yield metadata

    def del_metadata_by_source(self, source: str):
        """Remove a source from the feature

        :param source: Name of the source to delete

        :raises AssertionError: Source is not with the correct type
        :raises KeyError: Source does not belong in the MetaFeature
        """
        assert isinstance(
            source, str
        ), f"Source is not a string but with {type(source)}"
        if self._metadata_getter.pop(source, None) is None:
            logging.getLogger("PPanGGOLiN").warning(
                "The source to remove does not exist"
            )

    def del_metadata_by_attribute(self, **kwargs):
        """Remove a source from the feature"""
        for source, metadata_dict in self._metadata_getter.items():
            for attr, value in kwargs.items():
                for meta_id, metadata in metadata_dict.items():
                    if hasattr(metadata, attr):
                        # BUG If value is a list, the join block detection.
                        # It would be better to keep a list and change in writing and reading metadata to join the list
                        if (
                            getattr(metadata, attr, None) in value
                            or getattr(metadata, attr, None) == value
                        ):
                            del self._metadata_getter[source][meta_id]

    def max_metadata_by_source(self) -> Tuple[str, int]:
        """Get the maximum number of metadata for one source

        :return: Name of the source with the maximum annotation and the number of metadata corresponding
        """
        max_source, max_meta = max(
            self._metadata_getter.items(), key=lambda x: len(x[1])
        )
        return max_source, len(max_meta)

    def has_metadata(self) -> bool:
        """
        Does the feature has some metadata associated.

        :return: True if it has metadata else False
        """

        return self.number_of_metadata > 0

    def has_source(self, source: str) -> bool:
        """Check if the source is in the metadata feature

        :param source: name of the source

        :return: True if the source is in the metadata feature else False
        """
        return source in self._metadata_getter
