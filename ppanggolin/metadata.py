#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import List, Union

# installed libraries
from pandas import isna

class Metadata:
    """
    This represents a metadata lionk to gene families organisms or ...

    :param source: source of the metadata
    :param name: name of the annotation/function
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
        assert any(x is not None for x in [score, e_val])
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