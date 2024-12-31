from time import time
from pathlib import Path
from numpy import logspace, zeros, log10, ravel
from numpy.random import random, seed
from scipy.sparse import random as srandom, csc_array
from scipy.sparse.linalg import spsolve
from pypardiso import spsolve as pspsolve
from casadi import DM, solve, Sparsity
from matplotlib import pyplot

N = 25000
NUM = 100
TIME_LIMIT = 10
FIG_PATH = Path(__file__).parents[3] / "doc" / "source" / "figures"

SCIPY, CASADI, CAS_SCI, PYPAR = range(4)

def create_objects(size, density):
    a_s = srandom(size, size, density, format="csc", dtype=float)
    a_s.setdiag(random(size) + size)
    b = random(size)
    sparsity = Sparsity(size, size, a_s.indptr, a_s.indices, True)
    a_c = DM(sparsity, a_s.data)
    return a_s, a_c, b

def conv_solve(a, b):
    # Convert from Casadi, solve in Scipy
    a = csc_array(a)
    spsolve(a, b)


def main():
    def measure(a, b, which):
        nonlocal functions
        func = functions[which]
        if func is None:
            return float("nan")
        start = time()
        func(a, b)
        result = time() - start
        if result > TIME_LIMIT:
            functions[which] = None
        return result

    functions = {
        SCIPY: spsolve,
        CASADI: solve,
        CAS_SCI: conv_solve,
        PYPAR: pspsolve}

    seed(0)  # reproducible random numbers
    sizes = logspace(1, log10(N), num=NUM)
    densities = [0.01, 0.02]
    times = zeros((NUM, 4, 2))
    for k, size in enumerate(sizes):
        size = int(size)
        active = [i for i, f in functions.items() if f is not None]
        if not active:
            times[k:, :, :] = float("nan")
            break
        for d, density in enumerate(densities):
            a_s, a_c, b = create_objects(size, density)
            times[k][SCIPY][d] = measure(a_s, b, SCIPY)
            times[k][CASADI][d] = measure(a_c, b, CASADI)
            times[k][CAS_SCI][d] = measure(a_c, b, CAS_SCI)
            times[k][PYPAR][d] = measure(a_s, b, PYPAR)
        print(k, size, ravel(times[k, :, 1]))


    pyplot.loglog(sizes, times[:, SCIPY, 0], "k-", label="Scipy, $\\varrho=0.01$")
    pyplot.loglog(sizes, times[:, CASADI, 0], "r-", label="Casadi")
    pyplot.loglog(sizes, times[:, CAS_SCI, 0], "g-", label="Casadi $\\to$ Scipy")
    pyplot.loglog(sizes, times[:, PYPAR, 0], "b-", label="pypardiso")
    pyplot.loglog(sizes, times[:, SCIPY, 1], "k--", label="$\\varrho=0.02$")
    pyplot.loglog(sizes, times[:, CASADI, 1], "r--")
    pyplot.loglog(sizes, times[:, CAS_SCI, 1], "g--")
    pyplot.loglog(sizes, times[:, PYPAR, 1], "b--")
    pyplot.grid()
    pyplot.legend(loc="best")
    pyplot.xlabel("System size")
    pyplot.ylabel("Computation time [s]")
    pyplot.savefig(FIG_PATH / "solver_comparison.png", bbox_inches="tight")
    pyplot.show()



if __name__ == '__main__':
    main()