"""Test models for unit testing"""

import logging
from simu import Model
from simu.core.thermo import InitialState
from simu.core.thermo.material import MaterialSpec, MaterialDefinition
from simu.core.thermo.species import SpeciesDefinition
from simu.app.thermo.factories import ExampleThermoFactory
from simu.core.thermo import ThermoParameterStore
from simu.core.utilities import SymbolQuantity, Quantity


RK_LIQ = "Boston-Mathias-Redlich-Kwong-Liquid"


class SimpleParameterTestModel(Model):
    """A simple model to test parameters"""

    def interface(self):
        with self.parameters as p:
            p.define("length", 10, "m")

    def define(self):
        pass


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


class MaterialTestModel(Model):
    def interface(self):
        spec = MaterialSpec(["H2O", "*"])
        self.materials.define_port("inlet", spec)

    def define(self):
        test_material = define_a_material(["H2O", "NO2"])
        _ = self.materials["inlet"]
        _ = self.materials.create_flow("local", test_material)


class MaterialTestModel2(Model):
    def interface(self):
        spec = MaterialSpec(["H2O", "*"])
        self.materials.define_port("inlet", spec)

    def define(self):
        inlet = self.materials["inlet"]
        _ = self.materials.create_state("local", inlet.definition)


class MaterialTestModel3(Model):
    def interface(self):
        pass

    def define(self):
        test_material = define_a_material(["H2O", "NO2"])
        self.materials.create_flow("local", test_material)


class ResidualTestModel(Model):
    """A simple model to test residuals"""

    def interface(self):
        self.parameters.define("length", 10, "m")
        self.parameters.define("area_spec", 50, "m**2")

    def define(self):
        res = self.parameters["length"] ** 2 - self.parameters["area_spec"]
        self.residuals.add("area", res, "m**2")


class ResidualTestModel2(Model):
    def interface(self):
        pass

    def define(self):
        res = SymbolQuantity("Hubert", "K")
        self.residuals.add("Hubert", res, "degC")


class SquareTestModel(Model):
    def __init__(self):
        super().__init__()
        self.no2sol = define_a_material(["CH3-CH2-CH3", "CH3-(CH2)2-CH3"])

    def interface(self):
        with self.parameters as p:
            p.define("T", 10, "degC")
            p.define("p", 10, "bar")
            p.define("N", 1, "mol/s")
            p.define("x_c3", 10, "%")

    def define(self):
        flow = self.materials.create_flow("local", self.no2sol)  # 4 DOF
        flow["N"] = flow["n"].sum()
        flow["x"] = flow["n"] / flow["N"]

        param, radd = self.parameters, self.residuals.add
        radd("N", flow["N"] - param["N"], "mol/s")
        res = flow["N"] * param["x_c3"] - flow["n"]["CH3-CH2-CH3"]
        radd("x", res, "mol/s")
        radd("T", flow["T"] - param["T"], "K")
        radd("p", flow["p"] - param["p"], "bar")


class MaterialParentTestModel(Model):

    class Child(Model):
        def __init__(self, mat_def):
            super().__init__()
            self.mat_def = mat_def

        def interface(self):
            self.materials.define_port("port")

        def define(self):
            flow = self.materials.create_flow("m_child", self.mat_def)

    def define(self):
        mat_def = define_a_material(["H2O", "NO2"])
        flow = self.materials.create_flow("m_parent", mat_def)
        with self.hierarchy.add("child", self.Child, mat_def) as child:
            child.materials.connect("port", flow)


def define_a_material(species) -> MaterialDefinition:
    """Defines a material to use. Normally, this would be a singelton somewhere
    in the project."""

    factory = ExampleThermoFactory()
    speciesdb = {s: SpeciesDefinition(s) for s in species}
    frame = factory.create_frame(speciesdb, RK_LIQ)
    store = ThermoParameterStore()
    initial_state = InitialState.from_std(len(species))
    initial_state.temperature = Quantity("10 degC")
    initial_state.pressure = Quantity("10 bar")
    return MaterialDefinition(frame, initial_state, store)
