# -*- coding: utf-8 -*-

# stdlib modules
from json import dumps, load, loads
from sys import _getframe, exc_info
from types import TracebackType
from pathlib import Path
from difflib import Differ

def user_agree(message):
    try:
        ans = input(f"{message} y/[n]? ")
        return ans.lower() == "y"
    except OSError:  # run automatically with std streams caught
        return False

def assert_reproduction(data, suffix=None):
    """Assert the json-dump of the data to be the same as before.
    This method will (if run interactively) ask the user to accept the
    new or changed data, if no old reference data exists or if it differs.

    The ``sorted=True`` flag is set when dumping the data, allowing
    (nested) dictionaries to be compared correctly.
    """

    def user_agree(message):
        try:
            ans = input(f"{message} y/[n]? ")
            return ans.lower() == "y"
        except OSError:  # run automatically with std streams caught
            return False

    def load_file():
        # try to open file
        # if not exists, dump and return empty dictionary
        try:
            with open(filename, "r") as file:
                data = load(file)
        except FileNotFoundError:
            data = {}
            with open(filename, "w") as file:
                file.write(dumps(data))
        return data

    def save_data(data):
        ref_data_all[meth_name] = data
        with open(filename, "w") as file:
            file.write(dumps(ref_data_all, sort_keys=True, indent=2))

    caller_file = Path(_getframe(1).f_code.co_filename)
    filename = caller_file.absolute().parent / "refdata.json"
    ref_data_all = load_file()

    data = loads(dumps(data))  # to align and assure compatibility
    meth_name = _getframe(1).f_code.co_name  # get name of calling method
    meth_name = f"{caller_file.name}::{meth_name}"
    if suffix:
        meth_name = f"{meth_name}_{suffix}"
    ref_data = ref_data_all.get(meth_name, None)

    try:
        if not ref_data:
            msg = (f"No reference data exists for {meth_name}. " +
                   "The following data is generated now:\n\n" +
                   dumps(data, indent=2, sort_keys=True) +
                   "\n\nDo you accept this data")
            if user_agree(msg):
                save_data(data)
            else:
                raise AssertionError("Reference data rejected by user")
        else:
            ref = dumps(ref_data, indent=2, sort_keys=True)
            val = dumps(data, indent=2, sort_keys=True)
            try:
                assert ref == val
            except AssertionError:
                differ = Differ()
                diff = differ.compare(ref.splitlines(), val.splitlines())
                msg = (f"Deviation from reference data detected ({meth_name}):" +
                       "\n".join(diff) + "\n\nDo you accept the new data")
                if user_agree(msg):
                    save_data(data)
                else:
                    raise
    except AssertionError as err:
        traceback = exc_info()[2]
        back_frame = traceback.tb_frame.f_back
        back_tb = TracebackType(tb_next=None,
                                tb_frame=back_frame,
                                tb_lasti=back_frame.f_lasti,
                                tb_lineno=back_frame.f_lineno)
        raise err.with_traceback(back_tb)
