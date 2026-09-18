from typing import Sequence

from numpy import squeeze, array, atleast_1d, argmin, isfinite, argmax, abs
from numpy.typing import NDArray

from simu import NumericHandler, NHKeys, Quantity
from simu.core.utilities.types import NestedMap, NestedMutMap
from simu.core.utilities.errors import NonSquareSystem


class ModelContext:
    def __init__(self, model: NumericHandler):
        self._parameters = model.function.arg_structure.get(
            NHKeys.MODEL_PARAMS, {}
        )
        self._properties = model.function.result_structure.get(
            NHKeys.MODEL_PROPS, {}
        )

    def parameter_unit(self, path: Sequence[str]) -> str:
        return self._extract(path, self._parameters)

    def property_unit(self, path: Sequence[str]) -> str:
        return self._extract(path, self._properties)

    @staticmethod
    def _extract(path: Sequence[str], structure: NestedMap[str]) -> str:
        result = structure
        try:
            for p in path:
                result = result[p]
        except (KeyError, TypeError) as e:
            raise KeyError(f"Invalid path: '{'.'.join(path)}'") from e
        if not isinstance(result, str):
            raise KeyError(f"Invalid path: '{'.'.join(path)}'")
        return result


def relax(b: NDArray, a: NDArray,
          bound_names: Sequence[str],
          gamma: float) -> tuple[float, str]:
    a, b = [squeeze(array(x)) for x in (a, b)]
    # are there bounds violated?
    invalid = [n for n, m_i in zip(bound_names, atleast_1d(b <= 0)) if m_i]

    if invalid:
        msg = f"Bound violation of: {', '.join(invalid)}"
        raise ValueError(msg)

    mask = (a > 0)
    a = a[mask]
    alpha, min_alpha_name = 1.0, ""
    if not len(a):
        return alpha, min_alpha_name

    min_a_idx = int(argmin(a))
    if a[min_a_idx] * gamma < 1:
        alpha = a[min_a_idx] * gamma
        bn = [b for b, m in zip(bound_names, mask) if m]
        min_alpha_name = bn[min_a_idx]
    return alpha, min_alpha_name


def not_finite(vector: NDArray, names: Sequence[str]) -> Sequence[str]:
    finite = isfinite(vector)
    if False in isfinite(finite):
        return [names[i] for i, f in enumerate(finite) if not f]
    return []


def assess_residuals(vector: NDArray,
                     names: Sequence[str]) -> tuple[float, str]:
    if len(vector):
        # assess error
        max_err_idx = int(argmax(abs(vector)))
        max_name = names[max_err_idx]
        max_err = float(abs(vector[max_err_idx]))
        return max_err, max_name
    else:  # trivial model, nothing to solve
        return 0, ""


def check_model_square(model: NumericHandler):
    args = model.arguments
    # store size of state
    state_size = args[NHKeys.VECTORS][NHKeys.STATES].m.size()[0]
    res_size = len(model.vector_res_names(NHKeys.RESIDUALS))

    if state_size != res_size:
        raise NonSquareSystem(state_size, res_size)


class DataRowConverter:
    def __init__(self, from_uom: Sequence[str], to_uom: Sequence[str]):
        self._from = from_uom
        self._to = to_uom

    def __call__(self, row: Sequence[float]) -> Sequence[float]:
        return [
            Quantity(r, f).to(t).magnitude
            for r, f, t in zip(row, self._from, self._to)
    ]


def replace_qty(
        arguments: NestedMutMap[Quantity],
        item: Quantity, path: Sequence[str]
):
    """Replace an item"""
    prev = None
    for p in path:
        if not p in arguments:
            return  # Thermo-parameter is not in model, skip
        prev, arguments = arguments, arguments[p]
    prev[path[-1]] = item


def extract_qty(results: NestedMap[Quantity], path: Sequence[str]) -> Quantity:
    for p in path:
        results = results[p]
    return results
