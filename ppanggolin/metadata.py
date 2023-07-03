#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Tuple, Union

# installed libraries
from pandas import isna


class Metadata:
    """The Metadata class represents a metadata link to genes, gene families, organisms, regions, spot or modules. It allows the creation of metadata objects with different attributes and values, and provides methods to access and manipulate these attributes. The class has a constructor method that initializes the Metadata object with a source and a dictionary of attributes and values. The class also has methods to get the value of a specific attribute, return a list of all the attributes, and join a list of strings into a single string separated by commas. The class has two fields: source, which represents the source of the metadata, and **kwargs, which is a dictionary of attributes and values representing the metadata.

        Methods:
        - __init__(self, source: str, **kwargs): Constructor method that initializes the Metadata object with a source and a dictionary of attributes and values.
        - number_of_attribute(self): Returns the number of attributes in the Metadata object.
        - get(self, name: str, skip_error: bool = False): Returns the value of a specific attribute in the Metadata object, or None if the attribute does not exist. If skip_error is True, it does not raise an AttributeError if the attribute does not exist.
        - fields(self) -> List[str]: Returns a list of all the attributes in the Metadata object.
        - _join_list(attr_list: Union[str, List[str]]): Joins a list of strings into a single string separated by commas.

        Fields:
        - source: A string representing the source of the metadata.
        - **kwargs: A dictionary of attributes and values representing the metadata. The attributes can be any string, and the values can be any type except None or NaN.
    """
    def __init__(self, source: str, **kwargs):
        """
        The Metadata class represents a metadata link to genes, gene families, organisms, regions, spot or modules.
        It allows the creation of metadata objects with different attributes and values, and provides methods to access
        and manipulate these attributes. Add attributes and values representing the metadata as mush as you want.
        The attributes can be any string, and the values can be any type except None or NaN.

        :param source: A string representing the source of the metadata.
        :param kwargs: A dictionary of attributes and values representing the metadata. The attributes can be any string, and the values can be any type except None or NaN.
        """
        self.source = source
        if len(kwargs) == 0:
            raise Exception(f"No metadata given for source: {source}")
        for attr, value in kwargs.items():
            if isinstance(value, list):
                value = self._join_list(value)
            if value is not None and not isna(value):
                setattr(self, attr, value)

    def number_of_attribute(self):
        return len(self.__dict__.keys())

    def get(self, name: str):
        try:
            value = self.__getattribute__(name)
        except AttributeError as attr_error:
                raise AttributeError(attr_error)
        else:
            return value

    @property
    def fields(self) -> List[str]:
        return list(self.__dict__.keys())

    @staticmethod
    def _join_list(attr_list: Union[str, List[str]]):
        return ','.join(attr_list)


class MetaFeatures:
    """
    This represents a methods to access metadata in genes, gene families, organisms, regions, spot or modules
    """
    def __init__(self):
        self._metadataGetter = {}

    @property
    def metadata(self) -> Generator[Metadata, None, None]:
        """Generate metadatas in gene families

        :return: Generator with all metadata from all sources
        """

        for meta_list in self._metadataGetter.values():
            for metadata in meta_list:
                yield metadata

    @property
    def sources(self) -> List[str]:
        """ Get all metadata source in gene family

        :return: List of metadata source
        """
        return list(self._metadataGetter.keys())

    def get_source(self, source: str) -> Union[List[Metadata], None]:
        """ Get the metadata for a specific source in gene family

        :param source: Name of the source

        :return: All the metadata from the source if exist else None
        """
        return self._metadataGetter[source] if source in self.sources else None

    def get_metadata(self, **kwargs) -> Generator[Metadata, None, None]:
        """Get metadata by one or more attribute

        :return: metadata searched
        """
        for metadata in self.metadata:
            for attr, value in kwargs.items():
                if hasattr(metadata, attr):
                    if metadata.__getattribute__(attr) in value or metadata.__getattribute__(attr) == value:
                        yield metadata

    def add_metadata(self, source: str, metadata: Metadata):
        """ Add metadata

        :param source: Name of database source
        :param metadata: Identifier of the metadata
        """
        assert isinstance(metadata, Metadata)
        source_annot = self.get_source(source)
        if source_annot is not None:
            self._metadataGetter[source].append(metadata)
        else:
            self._metadataGetter[source] = [metadata]

    def max_metadata_by_source(self) -> Tuple[str, int]:
        """Get the maximum number of annotation for one source
        :return: Name of the source with the maximum annotation and the number of annotation corresponding
        """
        max_meta = 0
        max_source = None
        for source, metadata in self._metadataGetter.items():
            if len(metadata) > max_meta:
                max_meta = len(metadata)
                max_source = source
        return max_source, max_meta