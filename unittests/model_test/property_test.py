"""Unit tests related to the property handler classes"""
from simu.model.property import PropertyDefinition
from simu.utilities import assert_reproduction, SymbolQuantity
from simu.model.utils import ModelStatus


def test_provide():
    """Provide a single property"""
    props = PropertyDefinition()
    props.provide("area", "m**2")
    assert "area" in props._PropertyDefinition__symbols


def test_define():
    """Provide and define a single property"""
    props = PropertyDefinition()
    props.provide("area", "m**2")
    props.status = ModelStatus.DEFINE
    length = SymbolQuantity("length", "m")
    props["area"] = length * length
    area = props["area"]
    assert_reproduction(area)
    return props


def test_create_handler():
    """Create a property handler object"""
    props = PropertyDefinition()
    props.provide("area", "m**2")
    props.status = ModelStatus.DEFINE
    length = SymbolQuantity("length", "m")
    props["area"] = length * length

    props.status = ModelStatus.READY
    handler = props.create_handler()
    return handler


def test_handler_status():
    handler = test_create_handler()
    sym = handler.local_symbols
    assert_reproduction(sym)
