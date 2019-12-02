#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the aichem module.
"""
import pytest

from aichem import aichem


def test_something():
    assert True


def test_with_error():
    with pytest.raises(ValueError):
        # Do something that raises a ValueError
        raise(ValueError)


# Fixture example
@pytest.fixture
def an_object():
    return {}


def test_aichem(an_object):
    assert an_object == {}
