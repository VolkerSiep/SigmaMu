"""Unit tests related to the parameter handler classes"""
from simu.model.parameter import ParameterDefinition
from simu.utilities import assert_reproduction, SymbolQuantity
from simu.model.utils import ModelStatus


def test_definition():
    param = ParameterDefinition()
    param.define("length")
    assert param.defined("length")


def test_definition_default():
    param = ParameterDefinition()
    param.define("length", 3.14159, "m")
    assert param.defined("length")


def test_definition_symbols():
    param = ParameterDefinition()
    param.define("length", 3.14159, "m")

    param.status = ModelStatus.DEFINE

    assert "length" in param.symbols

    quantity = param["length"]
    assert_reproduction(quantity)


def test_create_handler():
    param = ParameterDefinition()
    param.define("length", 3.14159, "m")
    param.status = ModelStatus.READY
    handler = param.create_handler()
    assert handler.defined("length")
    return handler


def test_update():
    handler = test_create_handler()
    handler.update(length="1 cm")
    handler.status = ModelStatus.FINALISED
    assert handler.free_symbols
    values = handler.values
    assert_reproduction(values)


def test_provide():
    handler = test_create_handler()
    length = SymbolQuantity("new_length", "light_year")
    handler.provide(length=length)
    handler.status = ModelStatus.FINALISED
    assert not handler.values
    assert not handler.free_symbols