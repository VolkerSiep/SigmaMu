from functools import cached_property

from numpy import asarray, hstack, squeeze
import casadi as cas
from pint.facets.plain import PlainUnit


from .structures import flatten_dictionary, unflatten_dictionary
from .quantity import Quantity, simplify_quantity
from .types import Map, NestedMap, MutMap, NestedMutMap


class QFunction:
    """Wrapper around ``casadi.Function`` to consider units of measurements.
    This derived function object is defined by a dictionary of arguments and
    results, both of which are ``Quantity`` instances with casadi symbols as
    magnitudes. These symbols are independent symbols for ``args`` and derived
    symbols for ``results``.

    The function object is then called as well with a dictionary of
    ``Quantity`` values. The units of measurement must be consistent, and a
    conversion is done to the initially defined units for the individual
    arguments. The result is given back as a dictionary of ``Quantity`` values
    in the same units as initially defined.
    """

    def __init__(self, args: NestedMap[Quantity], results: NestedMap[Quantity],
                 func_name: str = "f", simplify_units: bool = True):
        args_flat = flatten_dictionary(args)
        arg_sym = cas.vertcat(*[v.magnitude for v in args_flat.values()])
        results_flat = flatten_dictionary(results)
        if simplify_units:
            results_flat = {k: simplify_quantity(v)
                            for k, v in results_flat.items()}

        self.__res_shapes = {}
        res_sym = []
        for key, value in results_flat.items():
            mag = value.magnitude
            try:
                shape = mag.shape
            except AttributeError:
                shape = (1, 1)
            self.__res_shapes[key] = shape
            res_sym.append(cas.reshape(mag, shape[0] * shape[1], 1))
        res_sym = cas.vertcat(*res_sym)

        self.arg_units = {k: v.units for k, v in args_flat.items()}
        self.res_units = {k: v.units for k, v in results_flat.items()}
        self.func = cas.Function(func_name, [arg_sym], [res_sym], ["x"], ["y"])

    def __call__(self, args: NestedMap[Quantity]) -> NestedMutMap[Quantity]:
        """Call operator for the function object, as described above.

        :param args: The (partially) symbolic arguments of the function
        """
        args_flat = cas.vertcat(*self.prepare_arguments(args).values())
        result = self.func(args_flat)  # calling CasADi function
        result = self.__extract_entities(result)
        result = {k: Quantity(v, self.res_units[k])
                  for k, v in result.items()}
        return unflatten_dictionary(result)

    def evaluate(self, args: NestedMap[Quantity]) -> NestedMutMap[Quantity]:
        """Call operator for the function object, as described above.

        :param args: The numeric arguments of the function
        """
        return self.evaluate_prepared(self.prepare_arguments(args))

    def evaluate_prepared(self, args: MutMap[float]) -> NestedMutMap[Quantity]:
        args_flat = hstack(list(args.values()))
        result = self.func(args_flat)  # calling CasADi function
        result = self.__extract_entities(result)
        result = {k: squeeze(asarray(v)) for k, v in result.items()}
        result = {k: Quantity(v, self.res_units[k])
                  for k, v in result.items()}
        return unflatten_dictionary(result)

    def prepare_arguments(self, args: NestedMap[Quantity]) -> MutMap[float]:
        """Flatten the arguments and convert the quantities into the units
        that are expected by the inner function definition. Using this method
        when reusing large parts of the argument can speed up the evaluation
        process.
        """
        return {
            key: value.to(self.arg_units[key]).magnitude
            for key, value in flatten_dictionary(args).items()
        }

    def __extract_entities(self, raw_result: cas.DM) -> Map[cas.DM]:
        result: MutMap[cas.DM] = {}
        idx = 0
        for key, shape in self.__res_shapes.items():
            size = shape[0] * shape[1]
            result[key] = raw_result[idx:idx + size].reshape(shape)
            idx += size
        return result

    @cached_property
    def result_structure(self) -> NestedMutMap[str]:
        """Return the result structure as a nested dictionary, only including
        the units of measurements as values of end nodes"""
        simplify = self.__simplify_unit
        units = {k: f"{simplify(v):~}" for k, v in self.res_units.items()}
        return unflatten_dictionary(units)

    @cached_property
    def arg_structure(self) -> NestedMutMap[str]:
        """Return the argument structure as a nested dictionary, only including
        the units of measurements as values of end nodes"""
        simplify = self.__simplify_unit
        units = {k: f"{simplify(v):~}" for k, v in self.arg_units.items()}
        return unflatten_dictionary(units)

    @staticmethod
    def __simplify_unit(unit: PlainUnit) -> PlainUnit:
        return simplify_quantity(Quantity(1, unit)).units
