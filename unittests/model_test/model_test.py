"""Unit tests related to the model class"""
from simu import Model
from simu.utilities.testing import assert_reproduction


class SquareTestModel(Model):
    """A simple model that calculates the square of a parameter as a surface"""

    def interface(self):
        """Here I can nicely document the inteface of the model"""
        self.parameters.define("length", 10.0, "m")
        self.properties.provide("area", unit="m**2")

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
        self.properties.provide("volume", unit="m**3")

    def define(self):
        child = self.hierarchy.add("square", SquareTestModel())

        child = self.hierarchy["square"] = SquareTestModel()
        volume = child.properties["area"] * self.parameters["depth"]
        self.properties["volume"] = volume


def test_square():
    """Test to instantiate the square test model and check symbols"""
    model = SquareTestModel().instance().finalise()

    area = model.properties["area"]
    length = model.parameters.symbols["length"]
    result = {"length": length, "area": area}
    assert_reproduction(result)


def test_square_numerics():
    """Test evaluating a simple model with just parameters and properties"""
    model = SquareTestModel().finalise()
    model.parameters.update(length="10 cm")

    numerics = model.numerics
    # the numerics interface wraps the entire model into a new function
    # and provides numerical structures for parameters, states, residuals, etc.
    param = numerics.parameters
    assert_reproduction(param, suffix="param")
    #  should be something like {"length": Q("10 cm")}
    numerics.evaluate()  # I probably need that to evaluate the casadi function
    props = numerics.properties
    #  should be something like {"area": Q("100 cm**2")}
    assert_reproduction(props, suffix="props")


def test_hierarchy():
    """Test evaluating a simple hierarchical model with just parameters and
    properties"""
    model = HierarchyTestModel().finalise()
    numerics = model.numerics
    param = numerics.parameters
    assert_reproduction(param, suffix="param")

    numerics.evaluate()
    props = numerics.properties
    assert_reproduction(props, suffix="props")
