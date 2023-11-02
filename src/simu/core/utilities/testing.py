# -*- coding: utf-8 -*-
"""auxiliary routines for unit testing"""

# stdlib modules
from json import dumps, load, loads, JSONEncoder
from inspect import currentframe, getouterframes
from pathlib import Path
from difflib import Differ

from casadi import SX

SX.__json__ = SX.__str__
FILENAME = "refdata.json"


class CustomEncoder(JSONEncoder):
    """Custom encoder for Objects"""

    def default(self, o):
        """Try to json-encode objects by calling their ``__json__`` method"""
        try:
            return o.__json__()
        except AttributeError:
            return super().default(o)


def user_agree(message):
    """Ask for confirmation, but assume "no", if not run interactively"""
    try:
        ans = input(f"{message} y/[n]? ")
        return ans.lower() == "y"
    except (OSError, EOFError):  # run automatically with std streams caught
        return False


def assert_reproduction(result, suffix=None):
    """Assert the json-dump of the data to be the same as before.
    This method will (if run interactively) ask the user to accept the
    new or changed data, if no old reference data exists or if it differs.

    The ``sorted=True`` flag is set when dumping the data, allowing
    (nested) dictionaries to be compared correctly.
    """

    def load_file():
        """try to open refdata file. If it doesn't exist, dump and return
        an empty dictionary"""
        try:
            with open(filename, "r", encoding="utf-8") as file:
                data = load(file)
        except FileNotFoundError:
            data = {}
            with open(filename, "w", encoding="utf-8") as file:
                file.write(dumps(data))
        return data

    frame = getouterframes(currentframe())[1]
    caller_file = Path(frame.filename)
    filename = caller_file.absolute().parent / FILENAME
    ref_data_all = load_file()

    # to align and assure compatibility
    result = loads(dumps(result, cls=CustomEncoder))
    func_name = frame.function  # get name of calling function
    func_name = f"{caller_file.name}::{func_name}"
    if suffix:
        func_name = f"{func_name}_{suffix}"
    ref_data = ref_data_all.get(func_name, None)

    def save_data(data):
        """Save the reference data to the file"""
        ref_data_all[func_name] = data
        with open(filename, "w", encoding="utf-8") as file:
            file.write(dumps(ref_data_all, sort_keys=True, indent=2))

    if ref_data is None:
        msg = (f"No reference data exists for {func_name}. " +
               "The following data is generated now:\n\n" +
               dumps(result, indent=2, sort_keys=True) +
               "\n\nDo you accept this data")
        if user_agree(msg):
            save_data(result)
        else:
            raise AssertionError("Reference data rejected by user")
    else:
        ref = dumps(ref_data, indent=2, sort_keys=True)
        val = dumps(result, indent=2, sort_keys=True)
        try:
            assert ref == val
        except AssertionError:
            differ = Differ()
            diff = differ.compare(ref.splitlines(), val.splitlines())
            msg = (
                f"Deviation from reference data detected ({func_name}):\n"
                + "\n".join(diff) + "\n\nDo you accept the new data")
            if user_agree(msg):
                save_data(result)
            else:
                raise
