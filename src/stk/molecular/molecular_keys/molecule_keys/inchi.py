"""
InChI
=====

"""

from .molecule_key import MoleculeKey
from .utilities import get_inchi


class Inchi(MoleculeKey):
    """
    Used to get the InChI of molecules.

    Examples
    --------
    You want to use the InChI as part of a JSON representation of a
    molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            molecule_keys=(stk.Inchi(), ),
        )
        # Get the JSON representation, including an InChI.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self):
        """
        Initialize a :class:`Inchi` instance.

        """

        super().__init__('InChI', get_inchi)
