"""
InChIKey
========

"""

from .molecule import MoleculeKeyMaker
from .utilities import get_inchi_key


class InchiKey(MoleculeKeyMaker):
    """
    Used to get the InChIKey of molecules.

    Examples
    --------
    *Adding InChIKey to a Molecule's JSON*

    You want to use the InChIKey as part of a JSON representation of a
    molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            key_makers=(stk.InchiKey(), ),
        )
        # Get the JSON representation, including an InChIKey.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self):
        """
        Initialize an :class:`InchiKey` instance.

        """

        super().__init__('InChIKey', get_inchi_key)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'InchiKey()'
