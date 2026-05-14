from typing import Sequence

from numpy import squeeze, array, atleast_1d, argmin, isfinite, argmax, abs
from numpy.typing import NDArray

from simu import NumericHandler, NHKeys
from simu.core.utilities.errors import NonSquareSystem


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