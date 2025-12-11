"""Simple module to define the version - such that it can be read from other
places as well."""

MAJOR = 1  #: Major version
MINOR = 0  #: Minor version
FLAG = "beta"  #: Version flag (None if release, else alpha or beta)
BUILD = 1  #: Build flag

version = f"{MAJOR}.{MINOR}"
if FLAG is not None:
    version = f"{version}-{FLAG[0]}{BUILD}"
