from ..utilities import flatten_dictionary


class NumericHandler:

    def __init__(self, parent):
        self.__parent = parent
        self.__result = None

    @property
    def parameters(self):
        params = flatten_dictionary({
            name: child.numerics.parameters
            for name, child in self.__parent.hierarchy.items()
        })
        params.update(self.__parent.parameters.values)
        return params

    def evaluate(self):
        # TODO: need to use function defined under prepare.
        #  maybe create public QFunction object
        args = {"parameters": self.parameters}
        self.__result = self.__parent.function(args)

    def prepare(self):
        # define function for the entire model
        # define parameters as flat array
        #   maybe even states
        #  ... and properties of course
        #  ... what about thermodynamic properties?
        pass

    @property
    def properties(self):
        # todo: dive into hierarchy to get child results
        return self.__result["properties"]