"""
InChI
=====

"""

from .molecule_key import MoleculeKey
from .utilities import get_inchi


class Inchi(MoleculeKey):
    """
    Used to get the InChI of molecules.

    """

    def __init__(self):
        """
        Initialize a :class:`Inchi` instance.

        """

        super().__init__('InChI', get_inchi)
