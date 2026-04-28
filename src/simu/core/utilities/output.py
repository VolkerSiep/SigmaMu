# stdlib
from sys import stdout
from typing import Tuple

# internal
from .types import Map, OutputIOStream


class ProgressTableOutput:
    """Based on defined column headings and printing format strings for the
    table elements, print the provided data successively to the given stream.

    The header is not printed before the first row of data is provided. The
    space taken by the row elements then determines the individual column
    widths. Example:

    >>> from sys import stdout
    >>> from simu import  Quantity

    >>> cols = {"m": ("Magnitude", "{:9.2g}"), "u": ("Unit", "{:>15s}")}
    >>> table = ProgressTableOutput(cols, stdout)
    >>> table.row(Quantity(20, "m/s"))
    Magnitude           Unit
    --------- --------------
           20 meter / second
    >>> table.row(Quantity(10, "degC"))
           10 degree_Celsius

    """
    def __init__(self, columns: Map[Tuple[str, str]], output: OutputIOStream):
        """
        The constructor configures the table based on the following parameters:

        :param columns:  The keys of the dictionary are the names of the
          attributes to tabulate. THe value tuple represents the column header
          (first element) and the formatting string (second element).
        :param output:  The used output stream. If ``None``, no output will be
          generated.
        """
        self.__cols = columns
        self._first = True
        self._write = (lambda t: None) if (output is None) else output.write

    def row(self, data: object):
        """Print a data row, whereas the attributes are extracted from the given
        ``object``. If this is the first row of the table, the headers will be
        printed with it, and the column width thereby determined.

        :param data: The object, which must have the attributes as defined as
          keys in the ``columns`` mapping given to the constructor.
        """
        write = self._write

        elem = [c[1].format(getattr(data, k)) for k, c in self.__cols.items()]

        if self._first:
            headings = [f"{c[0]:>{len(e)}s}"
                        for e, c in zip(elem, self.__cols.values())]
            write(" ".join(headings) + "\n")
            write(" ".join("-" * len(h) for h in headings) + "\n")
            self._first = False

        write(" ".join(elem) + "\n")
