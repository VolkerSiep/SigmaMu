from pytest import raises as pt_raises
from casadi import SX

from simu.core.utilities import (Quantity, SymbolQuantity, assert_reproduction,
                                 jacobian, qsum, log, exp, qsqrt, qpow, conditional,
                                 base_unit, QFunction, flatten_dictionary,
                                 unflatten_dictionary, extract_units_dictionary)


def test_symbol_quantity():
    """Test propagation of symbolic quantities"""
    x_1 = SymbolQuantity("x1", "m", ["A", "B"])
    x_2 = SymbolQuantity("x2", "kg", 2)
    a = SymbolQuantity("a", "1/s")
    y_1 = a * x_1
    y_2 = a * x_2
    assert_reproduction([y_1, y_2])


def test_quantity():
    """Test instantiation of quantities"""
    x = [
        Quantity(1, "cm"),
        Quantity("1 cm"),
        Quantity(1, "J/mol").to_base_units()
    ]
    x.append(Quantity(x[0]))
    assert_reproduction(x)


def test_jacobian():
    """Test Jacobian function for quantities"""
    x = SymbolQuantity("x1", "m", ["A", "B"])
    a = SymbolQuantity("a", "1/s")
    y = a * x
    z = jacobian(y, x)
    assert_reproduction(z)


def test_sum1():
    """Test sum function for quantities"""
    x = SymbolQuantity("x1", "m", "ABCDEFU")
    y = qsum(x)
    assert_reproduction(y)


def test_log():
    """Test log function for quantities"""
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        log(x1)

    z = log(x1 / x2)
    assert_reproduction(z)


def test_exp():
    """Test exp function for quantities"""
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        exp(x1)

    z = exp(x1 / x2)
    assert_reproduction(z)


def test_sqrt():
    """Test sqrt function for quantities"""
    x = SymbolQuantity("x", "m^2", "AB")
    z = qsqrt(x)
    assert_reproduction(z)


def test_pow():
    """Test pow function for quantities"""
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        qpow(x1, x2)

    z = qpow(x1 / x2, x2 / x1)
    assert_reproduction(z)


def test_conditional():
    """Test conditional function for branching"""
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    # try invalid condition
    cond: SX = x1 > x2  # this is not a boolean, but an SX node.
    z = conditional(cond, x1, x2)
    assert_reproduction(z)


def test_base_unit():
    """Test converting units into base-units"""
    candidates = [
        "cm", "hour", "kmol", "t/hr", "barg", "W/mol", "m**2/s", "C", "V", "pi"
    ]
    res = {c: base_unit(c) for c in candidates}
    assert_reproduction(res)


def test_qfunction():
    """test QFunction with flat parameters and results"""
    x = SymbolQuantity("x1", "m", ["A", "B"])
    a = SymbolQuantity("a", "1/s")
    y = a * x
    f = QFunction({"x": x, "a": a}, {"y": y})

    # and now with quantities
    x = Quantity([1, 2], "cm")
    a = Quantity("0.1 kHz")
    y = f({"x": x, "a": a})["y"]

    assert_reproduction(y)


def create_flat():
    orig = {"C": {"A": 1, "B": 2}, "A": 3}
    flat = flatten_dictionary(orig)
    return orig, flat


def test_simple_flatten():
    """Check that a simple dict is flattened as expected"""
    orig, flat = create_flat()
    assert_reproduction(flat)


def test_unflatten():
    """Check that original dict is reproduced by un-flattening"""
    orig, flat = create_flat()
    rep = unflatten_dictionary(flat)
    assert rep == orig


def test_flatten_escaped():
    """Check that (un-)flattening works with escaped separators"""
    orig = {"C": {"A/D": 1, "B": 2}, "A": 3}
    flat = flatten_dictionary(orig)
    assert_reproduction(flat)
    rep = unflatten_dictionary(flat)
    assert rep == orig


def test_q_function_nested():
    """test QFunction with nested parameters and results"""
    x = SymbolQuantity("x1", "m", ["A", "B"])
    a = SymbolQuantity("a", "1/s")
    y = a * x
    f = QFunction({"x": x, "b": {"a": a}}, {"z": {"y": y}})

    # and now with quantities
    x = Quantity([1, 2], "cm")
    a = Quantity("0.1 kHz")
    y = f({"x": x, "b": {"a": a}})["z"]["y"]
    assert_reproduction(y)


def test_extract_units_dictionary():
    """Test extraction of units from a nested dictionary"""
    x = Quantity([1, 2], "cm")
    a = Quantity("0.1 kHz")
    struct = {"x": x, "b": {"a": a}}
    units = extract_units_dictionary(struct)
    assert_reproduction(units)
