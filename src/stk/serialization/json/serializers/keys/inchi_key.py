"""
InChIKey
========

"""

from .molecule_keys import MoleculeKey
from .utilities import get_inchi_key


class InchiKey(MoleculeKey):
    """
    """

    def __init__(self):
        """
        Initialize an :class:`InchiKey` instance.

        """

        super().__init__('InChIKey', get_inchi_key)
