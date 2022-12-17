class NumericHandler:

    def __init__(self, parent):
        self.__parent = parent
        self.__result = None

    @property
    def parameters(self):
        # todo: dive into hierarchy to get child parameters
        return self.__parent.parameters.values

    def evaluate(self):
        args = {"parameters": self.parameters}
        self.__result = self.__parent.function(args)

    @property
    def properties(self):
        return self.__result["properties"]