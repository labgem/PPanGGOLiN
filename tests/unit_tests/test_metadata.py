#! /usr/bin/env python3

import pytest
from random import randint
from typing import Generator, Set

from ppanggolin.metadata import Metadata, MetaFeatures


class TestMetadata:
    @pytest.fixture
    def metadata(self) -> Generator[Metadata, None, None]:
        """Create a simple metadata"""
        yield Metadata("source", attribute1="value1", attribute2=["value2", "value3"])

    def test_constructor(self, metadata):
        """Tests that the Metadata object is created successfully with a valid source and attributes"""
        assert metadata.source == "source"
        assert metadata.attribute1 == "value1"
        assert metadata.attribute2 == "value2,value3"

    def test_constructor_with_empty_source_name(self):
        """Tests that a ValueError is raised when creating a Metadata object with an empty source name"""
        with pytest.raises(ValueError):
            Metadata("", attribute="value")

    def test_constructor_with_non_string_source_name(self):
        """Tests that a TypeError is raised when creating a Metadata object with a non-string source name"""
        with pytest.raises(TypeError):
            Metadata(123, attribute="value")

    def test_constructor_with_no_attributes(self):
        """Tests that an Exception is raised when creating a Metadata object with no attributes"""
        with pytest.raises(Exception):
            Metadata("source")

    def test_get_existing_attribute_value(self, metadata):
        """Tests that the value of an existing attribute is returned correctly"""
        assert metadata.attribute1 == "value1"

    def test_get_non_existing_attribute_value(self, metadata):
        """Tests that an AttributeError is raised when getting the value of a non-existing attribute"""
        with pytest.raises(AttributeError):
            _ = metadata.non_existing_attribute

    def test_attribute_fields(self, metadata):
        """Tests that the 'fields' method returns a list of all the attributes in the Metadata object"""
        assert metadata.fields == ["attribute1", "attribute2"]

    def test_length(self, metadata):
        """Tests that the number_of_attribute method returns the correct number of attributes in the Metadata object"""
        assert isinstance(len(metadata), int)
        assert len(metadata) == 2


class TestMetaFeatures:
    @pytest.fixture
    def metadata(self) -> Generator[Set[Metadata], None, None]:
        """Create a random number of metadata

        :return: Set of metadata
        """
        metadata = set()
        for i in range(randint(5, 10)):
            metadata.add(
                Metadata(
                    f"source_{i}", **{f"attr_{j}": j for j in range(randint(1, 5))}
                )
            )
        yield metadata

    @pytest.fixture
    def metafeatures(self, metadata) -> Generator[MetaFeatures, None, None]:
        """Create a simple metafeature object

        :return: metafeature fill with metadata
        """
        metafeatures = MetaFeatures()
        for meta in metadata:
            metafeatures.add_metadata(meta)
        yield metafeatures

    def test_add_metadata(self, metafeatures, metadata):
        """Tests that metadata can be added to the metadata getter"""
        assert all(
            list(metafeatures._metadata_getter[meta.source].values()) == [meta]
            for meta in metadata
        )

    def test_get_metadata_feature_corresponding_to_source(self, metafeatures, metadata):
        """Tests that all the metadata features corresponding to a source can be retrieved"""
        assert all(
            list(metafeatures.get_metadata_by_source(meta.source).values()) == [meta]
            for meta in metadata
        )

    def test_remove_source_from_feature(self, metafeatures):
        """Tests that a source can be removed from the feature"""
        metadata = Metadata("source_del", attribute1="value")
        metafeatures.add_metadata(metadata)
        metafeatures.del_metadata_by_source("source_del")
        assert metafeatures.get_metadata_by_source("source_del") is None

    def test_generate_all_metadata_sources(self, metafeatures, metadata):
        """Tests that all metadata sources can be generated"""
        assert list(metafeatures.sources) == [meta.source for meta in metadata]

    def test_get_metadata_by_attribute_values(self, metafeatures):
        """Tests that metadata can be retrieved based on attribute values"""
        meta = Metadata("source_test", attribute1="value_to_retrieve")
        # meta_list = Metadata("source_list", attribute1=["val_1", "val_2"])
        metafeatures.add_metadata(meta)
        # metafeatures[meta_list.source] = meta_list
        assert list(
            metafeatures.get_metadata_by_attribute(attribute1="value_to_retrieve")
        ) == [meta]
        # assert list(metafeatures.get_metadata(attribute1="val_1")) == [meta_list]

    def test_get_maximum_number_of_metadata_for_one_source(
        self, metafeatures, metadata
    ):
        """Tests that the maximum number of metadata for one source can be retrieved"""
        metadata1 = Metadata("source_max", attribute1="value1")
        metadata2 = Metadata("source_max", attribute2="value2")
        metafeatures.add_metadata(metadata1)
        metafeatures.add_metadata(metadata2)
        assert metafeatures.max_metadata_by_source() == ("source_max", 2)

    def test_metadata_is_not_with_type_metadata(self, metafeatures):
        """Tests that an AssertionError is raised when metadata is not with type Metadata"""
        with pytest.raises(AssertionError):
            metafeatures.add_metadata("not_metadata")
