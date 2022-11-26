from pytest import raises as pt_raises

from simu.utilities import (Quantity, SymbolQuantity, assert_reproduction,
                            jacobian, sum1, log, exp, sqrt, qpow, conditional,
                            base_unit, QFunction)


def qstr(q):
    return f"{q:~}"


def test_symbol_quantity():
    x_1 = SymbolQuantity("x1", "m", ["A", "B"])
    x_2 = SymbolQuantity("x2", "kg", 2)
    a = SymbolQuantity("a", "1/s")
    y_1 = a * x_1
    y_2 = a * x_2
    assert_reproduction([qstr(y_1), qstr(y_2)])


def test_quantity():
    x = [
        Quantity(1, "cm"),
        Quantity("1 cm"),
        Quantity(1, "J/mol").to_base_units()
    ]
    x.append(Quantity(x[0]))
    assert_reproduction(list(map(qstr, x)))


def test_jacobian():
    x = SymbolQuantity("x1", "m", ["A", "B"])
    a = SymbolQuantity("a", "1/s")
    y = a * x
    z = jacobian(y, x)
    assert_reproduction(qstr(z))


def test_sum1():
    x = SymbolQuantity("x1", "m", "ABCDEFU")
    y = sum1(x)
    assert_reproduction(qstr(y))


def test_log():
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        z = log(x1)

    z = log(x1 / x2)
    assert_reproduction(qstr(z))


def test_exp():
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        z = exp(x1)

    z = exp(x1 / x2)
    assert_reproduction(qstr(z))


def test_sqrt():
    x = SymbolQuantity("x", "m^2", "AB")
    z = sqrt(x)
    assert_reproduction(qstr(z))


def test_exp():
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    with pt_raises(TypeError):
        z = qpow(x1)

    z = qpow(x1 / x2, x2 / x1)
    assert_reproduction(qstr(z))


def test_conditional():
    x1 = SymbolQuantity("x1", "m", "AB")
    x2 = SymbolQuantity("x2", "cm", "AB")

    # try invalid condition
    cond = x1 > x2
    z = conditional(cond, x1, x2)
    assert_reproduction(qstr(z))


def test_base_unit():
    candidates = [
        "cm", "hour", "kmol", "t/hr", "barg", "W/mol", "m**2/s", "C", "V", "pi"
    ]
    res = {c: base_unit(c) for c in candidates}
    assert_reproduction(res)


def test_qfunction():
    x = SymbolQuantity("x1", "m", ["A", "B"])
    a = SymbolQuantity("a", "1/s")
    y = a * x
    f = QFunction({"x": x, "a": a}, {"y": y})

    # and now with quantities
    x = Quantity([1, 2], "cm")
    a = Quantity("0.1 kHz")
    y = f(x=x, a=a)["y"]

    assert_reproduction(qstr(y))
