"""Unit tests related to the model class"""

from pytest import raises
import logging

from simu import Model
from simu.utilities import assert_reproduction, SymbolQuantity
from simu.utilities.errors import DataFlowError, DimensionalityError
# from simu.model import Augmentor, MaterialSpec


# class SquareTestModel(Model):
#     """A simple model that calculates the square of a parameter as a surface"""
#
#     def interface(self):
#         """Here I can nicely document the interface of the model"""
#         self.parameters.define("length", 10.0, "m")
#         self.properties.declare("area", unit="m**2")
#
#     def define(self):
#         """Here I can document the internal function of the model.
#         I can use this distinction to include this doc-string only for
#         detailed documentations."""
#         length = self.parameters["length"]
#         self.properties["area"] = length * length
#


#
#
# class NestedHierarchyTestModel(Model):
#     """A new model to test nested hierarchy. Just redefining one of the
#     parameters, leaving the other one free."""
#
#     def interface(self):
#         self.parameters.define("radius", 10.0, "cm")
#         self.properties.declare("volume", "m**3")
#
#     def define(self):
#         with self.hierarchy.add("volumator", HierarchyTestModel2()) as child:
#             child.parameters.provide(radius=self.parameters["radius"])
#         self.properties["volume"] = child.properties["volume"]
#
#
# class MaterialTestModel(Model):
#
#     def interface(self):
#         gas_spec = MaterialSpec(species=["NO", "NO2", "O2", "*"])
#         # gas_spec.require(FancyAugmentor)
#         self.material.require("inlet_1", gas_spec)
#         self.material.require("inlet_2", gas_spec)
#
#     def define(self):
#         inlet_1 = self.material["inlet_1"]
#         inlet_2 = self.material["inlet_2"]
#         # should be defined upfront for project
#         gas: MaterialDefinition = "This is used as a template"
#         intermediate_1 = self.material.create("intermediate_1", gas)
#         intermediate_2 = self.material.create("intermediate_2", inlet_1.definition)


class ParameterTestModel(Model):
    """A simple model to test parameters"""

    def interface(self):
        pdef = self.parameters.define
        pdef("length", 10, "m")
        pdef("width", unit="m")

    def define(self):
        par = self.parameters
        area = par["length"] * par["width"]
        logging.debug(f"area = {area:~}")


def test_parameters_update(caplog):
    caplog.set_level(logging.DEBUG)
    with ParameterTestModel.proxy() as proxy:
        proxy.parameters.update("width", 2, "m")
    assert "area = (length*width) m ** 2" in caplog.text
    assert "width" in proxy.parameters.free
    assert "length" in proxy.parameters.free


def test_parameters_provide():
    width = SymbolQuantity("width", "m")
    with ParameterTestModel.proxy() as proxy:
        proxy.parameters.provide(width=width)
    assert "width" not in proxy.parameters.free
    assert "length" in proxy.parameters.free


def test_parameter_error():
    with raises(KeyError) as err:
        with ParameterTestModel.proxy("Peter") as proxy:
            proxy.parameters.update("Hansi", 2, "m")
    assert "Parameter 'Hansi' not defined in 'Peter'" in str(err)


def test_parameter_missing():
    with raises(DataFlowError) as err:
        with ParameterTestModel.proxy("Peter"):
            pass
    assert "Model 'Peter' has unresolved parameters: 'width'" in str(err.value)


def test_parameter_double():
    width = SymbolQuantity("width", "m")
    with ParameterTestModel.proxy("Peter") as proxy:
        proxy.parameters.provide(width=width)
        with raises(KeyError) as err:
            proxy.parameters.update("width", 2, "m")
    assert "Parameter 'width' already provided in 'Peter'" in str(err.value)


def test_parameter_incompatible():
    temperature = SymbolQuantity("temperature", "K")
    with ParameterTestModel.proxy("Peter") as proxy:
        with raises(DimensionalityError):
            proxy.parameters.provide(width=temperature)
        with raises(DimensionalityError):
            proxy.parameters.update("width", 100, "K")
        proxy.parameters.update("width", 100, "cm")


class PropertyTestModel(Model):
    """A simple model to test properties"""

    def interface(self) -> None:
        self.parameters.define("length", 10, "m")
        self.properties.declare("area", "m^2")

    def define(self) -> None:
        self.properties["area"] = self.parameters["length"] ** 2


def test_parameters_access():
    with PropertyTestModel.proxy() as proxy:
        pass
    area = proxy.properties["area"]
    assert f"{area:~}" == "sq(length) m ** 2"


def test_parameters_access_too_early():
    with PropertyTestModel.proxy() as proxy:
        with raises(DataFlowError):
            area = proxy.properties["area"]


def test_parameters_dont_define():
    model = PropertyTestModel()
    model.define = lambda: None
    with raises(DataFlowError):
        with model.create_proxy():
            pass


def test_parameters_define_other():
    model = PropertyTestModel()

    def mydef():
        model.properties["area"] = model.parameters["length"]

    model.define = mydef
    with raises(DimensionalityError):
        with model.create_proxy():
            pass


class HierarchyTestModel(Model):
    """A simple hierarchical model, where the child module calculates the
    square (surface) of a parameter, and the parent model calculates the
    volume as function of the calculated surface and a depth parameter."""

    def interface(self):
        self.parameters.define("depth", 5.0, "cm")
        self.properties.declare("volume", unit="m**3")

    def define(self):
        with self.hierarchy.add("square", PropertyTestModel) as child:
            pass
        volume = child.properties["area"] * self.parameters["depth"]
        self.properties["volume"] = volume


class HierarchyTestModel2(Model):
    """A simple hierarchical model, where the parent model calculates the
    length of a parameter, and the child model calculates the
    surface as function of the calculated length."""

    def interface(self):
        # make this one available via hierarchy proxy.
        #  declared class can be that one or a subclass!?
        self.hierarchy.declare("square", PropertyTestModel)

        self.parameters.define("radius", 5.0, "cm")
        self.parameters.define("depth", 10, "cm")
        self.properties.declare("volume", "m**3")

    def define(self):
        radius = self.parameters["radius"]
        length = 2 * radius
        # # later shortcut syntax proposal:
        # with self.child("square", SquareTestModel) as child:
        with self.hierarchy.add("square", PropertyTestModel) as child:
            child.parameters.provide(length=length)  # TODO: Why type error?

        volume = self.parameters["depth"] * child.properties["area"]
        self.properties["volume"] = volume


def test_hierarchy():
    """Test evaluating a simple hierarchicals model with just parameters and
    properties"""
    proxy = HierarchyTestModel.top()
    volume = proxy.properties["volume"]
    assert f"{volume:~}" == "(sq(length)*depth) cm * m ** 2"


def test_hierarchy2():
    """Test evaluating a simple hierarchicals model with just parameters and
    properties"""
    proxy = HierarchyTestModel2.top()
    volume = proxy.properties["volume"]
    assert f"{volume:~}" == "(depth*sq((2*radius))) cm ** 3"



# def test_square():
#     """Test to instantiate the square test model and check symbols"""
#     instance, _ = SquareTestModel.as_top_model()
#
#     area = instance.properties["area"]
#     length = instance.parameters.free["length"]
#     result = {"length": length, "area": area}
#     assert_reproduction(result)
#
#
# def test_two_instances():
#     """Check that two instances don't interfere"""
#     model = SquareTestModel()
#     num = 3
#     instances = [model.create_proxy(f"m_{i}") for i in range(num)]
#     lengths = [SymbolQuantity(f"l_{i}", "m") for i in range(num)]
#     for i in range(num):
#         with instances[i] as unit:
#             unit.parameters.provide(length=lengths[i])
#     res = [instance.properties["area"].magnitude for instance in instances]
#     assert_reproduction(res)
#
#
# @mark.skip(reason="Didn't reimplement numerics interface yet")
# def test_square_numerics():
#     """Test evaluating a simple model with just parameters and properties"""
#     instance = SquareTestModel().instance()
#     instance.parameters.update(length="10 cm")
#     instance.finalise()
#
#     numerics = instance.numerics.prepare()
#     param = numerics.parameters
#     assert_reproduction(param, suffix="param")
#     #  should be something like {"length": Q("10 cm")}
#     numerics.evaluate()  # evaluate the casadi function
#     props = numerics.properties
#     #  should be something like {"area": Q("100 cm**2")}
#     assert_reproduction(props, suffix="props")
#
#

#
#
# def test_hierarchy_property():
#     """Check that calculated property has correct expression"""
#     instance, _ = NestedHierarchyTestModel.as_top_model()
#     vol = instance.hierarchy["volumator"]
#     area = vol.hierarchy["square"].properties["area"].magnitude
#     volume = instance.properties["volume"].magnitude
#     assert_reproduction([area, volume])
