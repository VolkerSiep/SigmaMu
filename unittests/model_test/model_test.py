from simu import Model
from simu.utilities.testing import assert_reproduction


class SquareTestModel(Model):

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

    def interface(self):
        self.parameters.define("depth", 5.0, "cm")
        self.properties.provide("volume", unit="m**3")

    def define(self):
        child = self.hierarchy["square"] = SquareTestModel()
        volume = child.properties["area"] * self.parameters["depth"]
        self.properties["volume"] = volume


def test_square():
    """Test to instantiate the square test model and check symbols"""
    model = SquareTestModel().finalise()

    area = model.properties["area"]
    length = model.parameters["length"]
    result = {"length": length, "area": area}
    assert_reproduction(result)


def test_square_numerics():
    model = SquareTestModel()
    model.parameters.update(length="10 cm")

    numerics = model.numerics
    # the numerics interface wraps the entire model into a new function
    # and provides numerical structures for parameters, states, residuals, etc.
    param = numerics.parameters
    #  should be something like {"length": Q("10 cm")}
    numerics.evaluate()  # I probably need that to evaluate the casadi function
    properties = numerics.properties
    #  should be something like {"area": Q("100 cm**2")}
    res = [param, properties]
    assert_reproduction(res)


def test_hierarchy():
    model = HierarchyTestModel().finalise()
    numerics = model.numerics
    param = numerics.parameters
    assert_reproduction(param, suffix="param")
