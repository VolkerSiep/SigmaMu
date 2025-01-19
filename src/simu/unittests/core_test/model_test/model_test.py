"""Unit tests related to the model class"""

from pytest import raises, mark

from simu.core.utilities.errors import DataFlowError, DimensionalityError
from .models import *


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


def test_parameters_access():
    with PropertyTestModel.proxy() as proxy:
        pass
    area = proxy.properties["area"]
    assert f"{area:~}" == "sq(length) m ** 2"


def test_parameters_access_too_early():
    with PropertyTestModel.proxy() as proxy:
        with raises(DataFlowError):
            _ = proxy.properties["area"]


def test_parameters_dont_define():
    model = PropertyTestModel()
    model.define = lambda: None
    with raises(DataFlowError):
        with model.create_proxy():
            pass


def test_parameters_define_other():
    model = PropertyTestModel()

    def my_def():
        """replace define"""
        model.properties["area"] = model.parameters["length"]

    model.define = my_def
    with raises(DimensionalityError):
        with model.create_proxy():
            pass


def test_hierarchy():
    """Test evaluating a simple hierarchical model with just parameters and
    properties"""
    proxy = HierarchyTestModel.top()
    volume = proxy.properties["volume"]
    assert f"{volume:~}" == "(sq(length)*depth) cm * m ** 2"


def test_hierarchy2():
    """Test evaluating a simple hierarchical model with just parameters and
    properties"""
    proxy = HierarchyTestModel2.top()
    volume = proxy.properties["volume"]
    assert f"{volume:~}" == "(depth*sq((2*radius))) cm ** 3"

@mark.parametrize("simple_material_definition", [["H2O", "NO2"]], indirect=True)
def test_material(simple_material_definition, material_test_model):
    material = simple_material_definition.create_flow()

    with material_test_model().create_proxy() as model:
        assert "inlet" in model.materials
        model.materials.connect("inlet", material)


@mark.parametrize("simple_material_definition", [["KMnO4"]], indirect=True)
def test_wrong_material(simple_material_definition, material_test_model):
    material = simple_material_definition.create_flow()
    with raises(ValueError) as err:
        with material_test_model().create_proxy() as model:
            model.materials.connect("inlet", material)
    assert "incompatible" in str(err)


@mark.parametrize("simple_material_definition", [["H2O", "NO2"]], indirect=True)
def test_material_reuse_def(simple_material_definition):
    material = simple_material_definition.create_flow()

    with MaterialTestModel2().create_proxy() as model:
        assert "inlet" in model.materials
        model.materials.connect("inlet", material)


def test_residual_def():
    _ = ResidualTestModel().top()


def test_residual():
    with ResidualTestModel2().create_proxy() as model:
        pass
    res = model.residuals["Hubert"]
    assert f"{res.value:~}" == "Hubert K"
