#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Tuple, Union

# installed libraries
from pandas import isna


class Metadata:
    """
    This represents a metadata link to genes, gene families, organisms, regions, spot or modules

    :param source: source of the metadata
    :param kwargs: all metadata name with there value
    """

    def __init__(self, source: str, **kwargs):
        """Constructor Method
        """
        self.source = source
        for attr, value in kwargs.items():
            if value is not None and not isna(value):
                if isinstance(value, list):
                    value = self._join_list(value)
                setattr(self, attr, value)

    def get(self, name: str):
        return self.__getattribute__(name)

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
        same_value = False
        if source_annot is not None:
            index_annot = 0
            while index_annot < len(source_annot):
                current_annot = source_annot[index_annot]
                for attr, value in metadata.__dict__.items():
                    if hasattr(current_annot, attr) and current_annot.__getattribute__(attr) == value:
                        same_value = True
                if not same_value:
                    index_annot += 1
                else:
                    break
            if not same_value:
                source_annot.append(metadata)
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