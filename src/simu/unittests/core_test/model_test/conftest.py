from pytest import fixture, mark

from simu import (
    Model, NumericHandler, StringDictThermoSource, SpeciesDefinition,
    ThermoParameterStore, InitialState, MaterialDefinition, Material,
    MaterialSpec)

from simu.app.thermo.factories import ExampleThermoFactory


@fixture(scope="session")
def material_model_function(material_test_model3):
    """Make a function out of a model defining materials"""
    proxy = material_test_model3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results


@fixture(scope="session")
def material_h2o_rk_liq() -> Material:
    rk_liq = "Boston-Mathias-Redlich-Kwong-Liquid"
    factory = ExampleThermoFactory()
    species = {"H2O": SpeciesDefinition("H2O")}
    frame = factory.create_frame(species, rk_liq)
    store = ThermoParameterStore()
    initial_state = InitialState.from_std(1)
    material_def = MaterialDefinition(frame, initial_state, store)
    material = material_def.create_flow()
    return material


@fixture(scope="session")
def material_def_h2o_with_param() -> MaterialDefinition:
    factory = ExampleThermoFactory()
    species_db = {"H2O": SpeciesDefinition("H2O")}
    frame = factory.create_frame(species_db, "Ideal-Solid")
    store = ThermoParameterStore()

    source = StringDictThermoSource({
        'H0S0ReferenceState': {
            's_0': {'H2O': '0 J / K / mol '},
            'dh_form': {'H2O': '0 J / mol '},
            'T_ref': '298.15 K', 'p_ref': '1 bar'},
        'LinearHeatCapacity': {'cp_a': {'H2O': '50.0 J / K / mol '},
                               'cp_b': {'H2O': '0 J / K ** 2 / mol '}},
        'ConstantGibbsVolume': {'v_n': {'H2O': '18 ml / mol '}}})
    store.add_source("default", source)
    initial_state = InitialState.from_std(1)
    return MaterialDefinition(frame, initial_state, store)


@fixture(scope="session")
def thermo_param():
    data = {
        'H0S0ReferenceState': {
            's_0': {'CH3-(CH2)2-CH3': '0 J/(mol*K)',
                    'CH3-CH2-CH3': '0 J/(mol*K)'},
            'dh_form': {'CH3-(CH2)2-CH3': '0 kJ/mol',
                        'CH3-CH2-CH3': '0 kJ/mol'},
            'T_ref': '25 degC',
            'p_ref': '1 atm'},
        'LinearHeatCapacity': {
            'cp_a': {'CH3-(CH2)2-CH3': '98 J/(mol*K)',
                     'CH3-CH2-CH3': '75 J/(mol*K)'},
            'cp_b': {'CH3-(CH2)2-CH3': '0 J/(mol*K*K)',
                     'CH3-CH2-CH3': '0 J/(mol*K*K)'}},
        'CriticalParameters': {
            'T_c': {'CH3-(CH2)2-CH3': '425 K', 'CH3-CH2-CH3': '370 K'},
            'p_c': {'CH3-(CH2)2-CH3': '38 bar', 'CH3-CH2-CH3': '42.5 bar'},
            'omega': {'CH3-CH2-CH3': 0.199, 'CH3-(CH2)2-CH3': 0.153}},
        'MixingRule_A': {'T_ref': '25 degC'},
        'VolumeShift': {
            'c_i': {'CH3-(CH2)2-CH3': '0 m ** 3 / mol',
                    'CH3-CH2-CH3': '0 m ** 3 / mol'}},
        'BostonMathiasAlphaFunction': {
            'eta': {'CH3-CH2-CH3': 0, 'CH3-(CH2)2-CH3': 0}}
    }
    return StringDictThermoSource(data)


@fixture(scope="session")
def material_test_model_4(material_def_h2o_with_param):
    """Need to make a fixture for this model, as it needs another fixture
    for the material definition"""
    class  MaterialTestModel4(Model):
        def interface(self):
            pass

        def define(self):
            self.materials.create_flow("flow1", material_def_h2o_with_param)
            self.materials.create_flow("flow2", material_def_h2o_with_param)

    return MaterialTestModel4


@fixture(scope="session")
def material_parent_test_model():
    material_definition = simple_material_definition_function(["H2O", "NO2"])

    class MaterialParentTestModel(Model):
        class Child(Model):
            def __init__(self, mat_def):
                super().__init__()
                self.mat_def = mat_def

            def interface(self):
                self.materials.define_port("port")

            def define(self):
                self.materials.create_flow("m_child", self.mat_def)

        def define(self):
            mat_def = material_definition
            flow = self.materials.create_flow("m_parent", mat_def)
            with self.hierarchy.add("child", self.Child, mat_def) as child:
                child.materials.connect("port", flow)

    return MaterialParentTestModel


@fixture()
def square_test_model():
    material_definition = simple_material_definition_function(
        ["CH3-CH2-CH3", "CH3-(CH2)2-CH3"])

    class SquareTestModel(Model):
        def __init__(self):
            super().__init__()
            self.no2sol = material_definition

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

            d_n = 2 * flow["n"]["CH3-(CH2)2-CH3"] - flow["n"]["CH3-CH2-CH3"]
            self.bounds.add("delta n", d_n)

    return SquareTestModel


@fixture(scope="session")
def material_test_model3():
    material_definition = simple_material_definition_function(["H2O", "NO2"])


    class MaterialTestModel3(Model):
        def interface(self):
            pass

        def define(self):
            self.materials.create_flow("local", material_definition)

    return MaterialTestModel3


@fixture(scope="session")
def material_test_model():
    material_definition = simple_material_definition_function(["H2O", "NO2"])

    class MaterialTestModel(Model):
        def interface(self):
            spec = MaterialSpec(["H2O", "*"])
            self.materials.define_port("inlet", spec)

        def define(self):
            _ = self.materials["inlet"]
            _ = self.materials.create_flow("local", material_definition)

    return MaterialTestModel


@fixture(scope="session")
def simple_material_definition(request) -> MaterialDefinition:
    """Defines a material to use. Normally, this would be a singleton somewhere
    in the project."""
    return simple_material_definition_function(request.param)


def simple_material_definition_function(species) -> MaterialDefinition:
    """Defines a material to use. Normally, this would be a singleton somewhere
    in the project."""
    rk_liq = "Boston-Mathias-Redlich-Kwong-Liquid"
    factory = ExampleThermoFactory()
    species_db = {s: SpeciesDefinition(s) for s in species}
    frame = factory.create_frame(species_db, rk_liq)
    store = ThermoParameterStore()
    initial_state = InitialState.from_cbar(10.0, 10.0, [1.0] * len(species))
    return MaterialDefinition(frame, initial_state, store)

@fixture(scope="session")
def model_with_residual():
    def material_definition():
        ideal_liq = "Ideal-Liquid"
        factory = ExampleThermoFactory()
        species_db = {s: SpeciesDefinition(s)
                      for s in "H2O Na:1+ Cl:1-".split()}
        frame = factory.create_frame(species_db, ideal_liq)
        store = ThermoParameterStore()
        initial_state = InitialState.from_cbar(10.0, 1.0, [1.0, 0.1, 0.1])
        return MaterialDefinition(frame, initial_state, store)

    class ModelWithThermoResidual(Model):
        def interface(self):
            pass

        def define(self):
            definition = material_definition()
            self.materials.create_flow("liq", definition)

    return ModelWithThermoResidual