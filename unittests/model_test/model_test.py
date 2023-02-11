"""Unit tests related to the model class"""

from typing import Type
from pytest import mark

from simu import Model
from simu.utilities import assert_reproduction, SymbolQuantity


class SquareTestModel(Model):
    """A simple model that calculates the square of a parameter as a surface"""

    def interface(self):
        """Here I can nicely document the inteface of the model"""
        self.parameters.define("length", 10.0, "m")
        self.properties.declare("area", unit="m**2")

    def define(self):
        """Here I can document the internal function of the model.
        I can use this distinction to include this doc-string only for
        detailed documenations."""
        length = self.parameters["length"]
        self.properties["area"] = length * length


class HierarchyTestModel(Model):
    """A simple hierarchical model, where the child module calculates the
    square (surface) of a parameter, and the parent model calculates the
    volume as function of the calculated surface and a depth parameter."""

    def interface(self):
        self.parameters.define("depth", 5.0, "cm")
        self.properties.declare("volume", unit="m**3")

    def define(self):
        with self.hierarchy.add("square", SquareTestModel()) as child:
            pass

        volume = child.properties["area"] * self.parameters["depth"]
        self.properties["volume"] = volume


class HierarchyTestModel2(Model):
    """A simple hierarchical model, where the parent model calculates the
    length of a parameter, and the child model calculates the
    surface as function of the calculated length."""

    def interface(self):
        self.parameters.define("radius", 5.0, "cm")
        self.parameters.define("depth", 10, "cm")
        self.properties.declare("volume", "m**3")

    def define(self):
        length = 2 * self.parameters["radius"]

        # # later shortcut syntax proposal:
        # with self.child("square", SquareTestModel()) as child:
        with self.hierarchy.add("square", SquareTestModel()) as child:
            child.parameters.provide(length=length)

        volume = self.parameters["depth"] * child.properties["area"]
        self.properties["volume"] = volume


class NestedHierarchyTestModel(Model):
    """A new model to test nested hierarchy. Just redefining one of the
    parameters, leaving the other one free."""

    def interface(self):
        self.parameters.define("radius", 10.0, "cm")
        self.properties.declare("volume", "m**3")

    def define(self):
        with self.hierarchy.add("volumator", HierarchyTestModel2()) as child:
            child.parameters.provide(radius=self.parameters["radius"])
        self.properties["volume"] = child.properties["volume"]


class MaterialTestModel(Model):

    def interface(self):
        no_gas = MaterialSpec(species=["NO", "NO2", "O2", "*"])
        self.material.require("inlet", no_gas)

    def define(self):
        inlet = self.material["inlet"]
        intermediate = self.material.create("intermediate", NOxGas())
        intermediate_2 = self.material.create("intermediate_2", inlet.type)

hierarchy_models = [
    HierarchyTestModel, HierarchyTestModel2, NestedHierarchyTestModel
]


def test_square():
    """Test to instantiate the square test model and check symbols"""
    instance, _ = SquareTestModel.as_top_model()

    area = instance.properties["area"]
    length = instance.parameters.free["length"]
    result = {"length": length, "area": area}
    assert_reproduction(result)


def test_two_instances():
    """Check that two instances don't interfer"""
    model = SquareTestModel()
    num = 3
    instances = [model.create_proxy(f"m_{i}") for i in range(num)]
    lengths = [SymbolQuantity(f"l_{i}", "m") for i in range(num)]
    for i in range(num):
        instances[i].parameters.provide(length=lengths[i])
        instances[i].finalise()
    res = [instance.properties["area"].magnitude for instance in instances]
    assert_reproduction(res)


@mark.skip(reason="Didn't reimplement numerics interface yet")
def test_square_numerics():
    """Test evaluating a simple model with just parameters and properties"""
    instance = SquareTestModel().instance()
    instance.parameters.update(length="10 cm")
    instance.finalise()

    numerics = instance.numerics.prepare()
    param = numerics.parameters
    assert_reproduction(param, suffix="param")
    #  should be something like {"length": Q("10 cm")}
    numerics.evaluate()  # evaluate the casadi function
    props = numerics.properties
    #  should be something like {"area": Q("100 cm**2")}
    assert_reproduction(props, suffix="props")


@mark.parametrize("model_cls", hierarchy_models)
def test_hierarchy(model_cls: Type[Model]):
    """Test evaluating a simple hierarchicals model with just parameters and
    properties"""
    name = model_cls.__name__
    func = model_cls.as_top_model()[0].model.function
    assert_reproduction(func.arg_structure, f"{name}_arguments")
    assert_reproduction(func.result_structure, f"{name}_result")


def test_hierarchy_property():
    """Check that calculated property has correct expression"""
    instance, _ = NestedHierarchyTestModel.as_top_model()
    vol = instance.hierarchy["volumator"]
    area = vol.hierarchy["square"].properties["area"].magnitude
    volume = instance.properties["volume"].magnitude
    assert_reproduction([area, volume])
