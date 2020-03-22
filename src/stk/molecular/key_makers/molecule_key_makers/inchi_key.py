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
    You want to use the InChIKey as part of a JSON representation of a
    molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            molecule_keys=(stk.InchiKey(), ),
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
