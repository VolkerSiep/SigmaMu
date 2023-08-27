"""Test models for unit testing"""

import logging
from simu import Model


class ParameterTestModel(Model):
    """A simple model to test parameters"""

    def interface(self):
        with self.parameters as p:
            p.define("length", 10, "m")
            p.define("width", unit="m")

    def define(self):
        par = self.parameters
        area = par["length"] * par["width"]
        logging.debug(f"area = {area:~}")


class PropertyTestModel(Model):
    """A simple model to test properties"""

    def interface(self) -> None:
        self.parameters.define("length", 10, "m")
        self.properties.declare("area", "m^2")

    def define(self) -> None:
        self.properties["area"] = self.parameters["length"] ** 2


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
        self.hierarchy.declare("square", PropertyTestModel)
        with self.parameters as p:
            p.define("radius", 5.0, "cm")
            p.define("depth", 10, "cm")
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
