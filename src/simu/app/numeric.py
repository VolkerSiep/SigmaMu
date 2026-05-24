from collections.abc import Collection

from simu import PropertyFilter

class ExclusionFilter(PropertyFilter):
    """This class is a simple exclusion filter with two purposes:

    - to demonstrate how to subclass  :class:`~simu.PropertyFilter` with ease
    - to provide a concrete filter class that serves 80 % of all likely uses

    This filter only filters on property names, but not specifically on the
    ``sub_key`` attributes given. As such. one can for instance filter our all
    standard state chemical potentials ``mu_std``, but not just the entries for
    selected species.
    """
    def __init__(self, excluded_keys: Collection[str]):
        """Create a filter object based on excluded keys.

        :param excluded_keys: A collection of property identifiers to exclude
        """
        self._excluded_keys = excluded_keys

    def keep_property(self, name: str, sub_key: str = None) -> bool:
        return not name in self._excluded_keys
