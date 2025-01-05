from sys import stdout
from typing import Tuple

from .types import Map


class ProgressTableOutput:
    """Based on defined column headings and printing format strings for the
    table elements, print the provided data successively to the given stream.

    The header is not printed before the first row of data is provided. The
    space taken by the row elements then determines the individual column
    widths. Example:

    >>> from simu import  Quantity
    >>> from io import StringIO  # this is only required to make doctest happy.
    >>> stream = StringIO()  # else, we could just use the default stream here.

    >>> cols = {"m": ("Magnitude", "{:9.2g}"), "u": ("Unit", "{:>15s}")}
    >>> table = ProgressTableOutput(cols, stream=stream)
    >>> table.row(Quantity(20, "m/s"))
    >>> table.row(Quantity(10, "degC"))
    >>> print(stream.getvalue())
    Row Magnitude           Unit
    --- --------- --------------
      1        20 meter / second
      2        10 degree_Celsius

    """
    def __init__(self, columns: Map[Tuple[str, str]], row_dig : int|None = 3,
                 row_head = "Row", stream=stdout):
        """
        The constructor configures the table based on the following parameters:

        :param columns:  The keys of the dictionary are the names of the
          attributes to tabulate. THe value tuple represents the column header
          (first element) and the formatting string (second element).
        :param row_dig: The number of digits to be reserved for the row number.
          If ``None``, do not print the row number. Otherwise, the row will be
          the first column of the table.
        :param row_head: The heading for the row column.
        :param stream:  The used output stream. If ``None``, no output will be
          generated.
        """
        self.__cols = columns
        self._first = True
        self._row_fmt = None if row_dig is None else f"{{:{row_dig}d}}"
        self._row_head = row_head
        self._write = (lambda t: None) if (stream is None) else stream.write
        self._row = 1

    def row(self, data: object, row: int = None):
        """Print a data row, whereas the attributes are extracted from the given
        ``object``. If this is the first row of the table, the headers will be
        printed with it, and the column width thereby determined.

        :param data: The object, which must have the attributes as defined as
          keys in the ``columns`` mapping given to the constructor.
        :param row:  The row number to be printed. If omitted, the row number
          is incremented from the previous row, starting with one.
        """
        write = self._write
        row = self._row if row is None else row

        elem = [c[1].format(getattr(data, k)) for k, c in self.__cols.items()]
        row_txt = "" if self._row_fmt is None else self._row_fmt.format(row)

        if self._first:
            headings = [f"{c[0]:>{len(e)}s}"
                        for e, c in zip(elem, self.__cols.values())]
            if self._row_fmt is not None:
                row_head = f"{self._row_head:{len(row_txt)}s}"
                headings = [row_head] + headings

            write(" ".join(headings) + "\n")
            write(" ".join("-" * len(h) for h in headings) + "\n")
            self._first = False

        if self._row_fmt is not None:
            elem = [row_txt] + elem
        write(" ".join(elem) + "\n")
        self._row = row + 1
