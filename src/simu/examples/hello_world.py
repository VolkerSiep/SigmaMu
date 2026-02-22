from simu import Model, NumericHandler, Quantity, NHKeys


class Square(Model):
    """A model of a square"""

    def interface(self):
        self.parameters.define("length", 10, "m")
        self.properties.declare("area", "m^2")

    def define(self):
        self.properties["area"] = self.parameters["length"] ** 2


def main():
    numeric = NumericHandler(Square.top())
    func = numeric.function
    print(func.arg_structure)
    print(func.result_structure)

    args = numeric.arguments
    print(args)

    args[NHKeys.MODEL_PARAMS]["length"] = Quantity(20, "cm")
    result = func(args)

    print(f"{result[NHKeys.MODEL_PROPS]['area']:.3fP~}")


if __name__ == '__main__':
    main()

