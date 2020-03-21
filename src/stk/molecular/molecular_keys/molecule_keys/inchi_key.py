"""
InChIKey
========

"""

from .molecule_key import MoleculeKey
from .utilities import get_inchi_key


class InchiKey(MoleculeKey):
    """
    Used to get the InChIKey of molecules.

    Examples
    --------

    """

    def __init__(self):
        """
        Initialize an :class:`InchiKey` instance.

        """

        super().__init__('InChIKey', get_inchi_key)
