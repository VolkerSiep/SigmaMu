"""Simple module to define the version - such that it can be read from other
places as well."""

MAJOR = 1  #: Major version
MINOR = 0  #: Minor version
FLAG = "alpha"  #: Version flag (None if release, else alpha or beta)
BUILD = 4  #: Build flag

VERSION = f"{MAJOR}.{MINOR}"
if FLAG is not None:
    VERSION = f"{VERSION}-{FLAG[0]}{BUILD}"
