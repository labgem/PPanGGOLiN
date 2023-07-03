#! /usr/bin/env python3

import pytest
from random import choices, randint, sample
from typing import Generator, Set

from ppanggolin.metadata import Metadata, MetaFeatures


class TestMetadata:
    def test_create_metadata_with_attributes(self):
        metadata = Metadata('source', attribute='value')
        assert metadata.__getattribute__("attribute") == 'value'
        assert metadata.source == 'source'

    def test_create_metadata_with_no_attributes(self):
        with pytest.raises(Exception):
            Metadata(source='source')

    def test_get_existing_attribute_value(self):
        metadata = Metadata('source', attribute='value')
        assert metadata.get('attribute') == 'value'

    def test_get_non_existing_attribute_value(self):
        metadata = Metadata('source', attribute='value')
        with pytest.raises(AttributeError):
            metadata.get('non_existing_attribute')

    def test_get_all_attributes(self):
        metadata = Metadata('source', attribute='value', another_attribute='another_value')
        assert metadata.fields == ['source', 'attribute', 'another_attribute']

    def test_join_list_attribute(self):
        metadata = Metadata('source', attribute=['value1', 'value2'])
        assert metadata.get("attribute") == 'value1,value2'

    def test_metadata_number_of_attributes(self):
        metadata = Metadata('source', attribute='value', another_attribute='another_value')
        assert metadata.number_of_attribute() == 3


class TestMetaFeatures:
    #  Tests that metadata can be added to MetaFeatures and checking if it was added successfully
    def test_add_metadata(self):
        meta_features = MetaFeatures()
        metadata = Metadata('source1', attribute1='value1')
        meta_features.add_metadata('source1', metadata)
        assert meta_features._metadataGetter['source1'] == [metadata]

    #  Tests that metadata can be gotten from MetaFeatures by source and checking if it returns the correct metadata
    def test_get_metadata_by_source(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='value1')
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert meta_features.get_source('source1') == [metadata1]
        assert meta_features.get_source('source2') == [metadata2]

    #  Tests that metadata can be gotten from MetaFeatures by attribute and checking if it returns the correct metadata
    def test_get_metadata_by_attribute(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='value1')
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert list(meta_features.get_metadata(attribute1='value1')) == [metadata1]

    #  Tests that all metadata can be gotten from MetaFeatures and checking if it returns all metadata
    def test_get_all_metadata(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='value1')
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert list(meta_features.metadata) == [metadata1, metadata2]

    #  Tests that all metadata sources can be gotten from MetaFeatures and checking if it returns all sources
    def test_get_all_metadata_sources(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='value1')
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert meta_features.sources == ['source1', 'source2']

    #  Tests that the source with the maximum number of metadata can be gotten from MetaFeatures and checking if it returns the correct source and number
    def test_get_source_with_maximum_metadata(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='value1')
        metadata2 = Metadata('source1', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source1', metadata2)
        assert meta_features.max_metadata_by_source() == ('source1', 2)

    #  Tests that getting metadata from MetaFeatures with non-existent source returns None
    def test_get_metadata_with_non_existent_source_returns_none(self):
        meta_features = MetaFeatures()
        metadata = Metadata('source1', attribute1='value1')
        meta_features.add_metadata('source1', metadata)
        assert meta_features.get_source('source2') == None

    #  Tests that getting metadata from MetaFeatures with non-existent attribute returns None
    def test_get_metadata_with_non_existent_attribute_returns_none(self):
        meta_features = MetaFeatures()
        metadata = Metadata('source1', attribute1='value1')
        meta_features.add_metadata('source1', metadata)
        assert list(meta_features.get_metadata(attribute2='value2')) == []

    #  Tests that getting metadata from MetaFeatures with empty attribute value returns the correct metadata
    def test_get_metadata_with_empty_attribute_value_returns_correct_metadata(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1='')
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert list(meta_features.get_metadata(attribute1='')) == [metadata1]

    #  Tests that getting metadata from MetaFeatures with list attribute value returns the correct metadata
    def test_get_metadata_with_list_attribute_value_returns_correct_metadata(self):
        meta_features = MetaFeatures()
        metadata1 = Metadata('source1', attribute1=['value1', 'value2'])
        metadata2 = Metadata('source2', attribute2='value2')
        meta_features.add_metadata('source1', metadata1)
        meta_features.add_metadata('source2', metadata2)
        assert list(meta_features.get_metadata(attribute1='value1,value2')) == [metadata1]
