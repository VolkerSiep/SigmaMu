"""Simple module to define the version - such that it can be read from other
places as well."""

MAJOR : int = 0  #: Major version
MINOR : int = 2  #: Minor version
FLAG : None | str = None  #: Version flag (None if release, else alpha or beta)
BUILD : int = 1  #: Build flag

VERSION = f"{MAJOR}.{MINOR}"
if FLAG is not None:
    VERSION = f"{VERSION}-{FLAG[0]}{BUILD}"
