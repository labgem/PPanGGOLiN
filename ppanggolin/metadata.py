#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Tuple, Union

# installed libraries
from pandas import isna


class Metadata:
    """
    This represents a metadata lionk to gene families organisms or ...

    :param source: source of the metadata
    :param value: name of the annotation/function
    :param accession: accesion identifier
    :param secondary_names: Other possible name for the annotation
    :param description: description of the annotation
    :param e_val: E-value of the gene family/profile comparison
    :param score: Bit score of the gene family/profile comparison.
    :param bias: The biased composition score correction that was applied to the bit score
    """

    def __init__(self, source: str, value: Union[str, int, float], accession: str = None,
                 secondary_names: Union[str, List[str]] = None, description: str = None,
                 score: float = None, e_val: float = None, bias: float = None):
        """Constructor Method
        """
        self.source = source
        self.accession = accession
        self.value = value
        self.secondary_names = self._write_secondary_names(secondary_names)
        self.description = description
        self.score = score
        self.e_val = e_val
        self.bias = bias

    def __str__(self):
        return self.value

    @staticmethod
    def _write_secondary_names(secondary_names: Union[str, List[str]]):
        if isinstance(secondary_names, list):
            return ','.join(secondary_names)
        elif isinstance(secondary_names, str):
            return secondary_names
        elif secondary_names is None or isna(secondary_names):
            return ''


class MetaFeatures:

    def __init__(self):
        self._metadataGetter = {}

    @property
    def metadata(self) -> Generator[Metadata, None, None]:
        """Generate annotations in gene families

        :return: Gene family annotation"""
        for meta_list in self._metadataGetter.values():
            for metadata in meta_list:
                yield metadata

    @property
    def sources(self) -> List[str]:
        """ Get all metadata source in gene family

        :return: List of annotation source
        """
        return list(self._metadataGetter.keys())

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

    def get_source(self, name: str) -> Union[List[Metadata], None]:
        """ Get the annotation for a specific source in gene family

        :param name: Name of the source

        :return: All the annotation from the source if exist else None
        """
        return self._metadataGetter[name] if name in self.sources else None

    def get_metadata(self, value: Union[str, int, float, List[str], List[int], List[float]],
                     accession: Union[List[str], str]) -> Generator[Metadata, None, None]:
        """Get annotation by name or accession in gene family

        :param value: Names of annotation searched
        :param accession: Accession number of annotation searched

        :return: annotation searched
        """
        assert value is not None and accession is not None
        value = value if isinstance(value, list) else [value]
        accession = accession if isinstance(accession, list) else [accession]

        for annotation in self.metadata:
            if annotation.value in value or annotation.accession in accession:
                yield annotation

    def add_metadata(self, source: str, metadata: Metadata):
        """ Add annotation to gene family

        :param source: Name of database source
        :param metadata: Identifier of the annotation
        """
        source_annot = self.get_source(source)
        same_value = False
        if source_annot is not None:
            index_annot = 0
            insert_bool = False
            while index_annot < len(source_annot):
                current_annot = source_annot[index_annot]
                if current_annot.value == metadata.value:
                    same_value = True
                if current_annot.score is not None and metadata.score is not None:
                    if current_annot.score < metadata.score:
                        if same_value:
                            source_annot[index_annot] = metadata
                        else:
                            source_annot.insert(index_annot, metadata)
                            insert_bool = True
                    elif current_annot.score == metadata.score:
                        if current_annot.e_val is not None and metadata.e_val is not None:
                            if current_annot.e_val > metadata.e_val:
                                if same_value:
                                    source_annot[index_annot] = metadata
                                else:
                                    source_annot.insert(index_annot, metadata)
                                    insert_bool = True
                elif current_annot.e_val is not None and metadata.e_val is not None:
                    if current_annot.e_val > metadata.e_val:
                        if same_value:
                            source_annot[index_annot] = metadata
                        else:
                            source_annot.insert(index_annot, metadata)
                            insert_bool = True
                if not insert_bool and not same_value:
                    index_annot += 1
                else:
                    break
            if not insert_bool and not same_value:
                source_annot.append(metadata)
        else:
            self._metadataGetter[source] = [metadata]
