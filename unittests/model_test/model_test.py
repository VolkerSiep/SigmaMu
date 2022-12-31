"""Unit tests related to the model class"""
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

    def define(self):
        length = 2 * self.parameters["radius"]

        # # later shortcut syntax proposal:
        # with self.child("square", SquareTestModel()) as child:
        with self.hierarchy.add("square", SquareTestModel()) as child:
            child.parameters.provide(length=length)


def test_square():
    """Test to instantiate the square test model and check symbols"""
    instance = SquareTestModel().as_top_model()

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


def test_hierarchy():
    """Test evaluating a simple hierarchical model with just parameters and
    properties"""
    numerics = HierarchyTestModel().N

    param = numerics.parameters
    assert_reproduction(param, suffix="param")

    numerics.evaluate()
    props = numerics.properties
    assert_reproduction(props, suffix="props")
