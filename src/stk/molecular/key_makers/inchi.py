"""
InChI
=====

"""

from .molecule import MoleculeKeyMaker
from .utilities import get_inchi


class Inchi(MoleculeKeyMaker):
    """
    Used to get the InChI of molecules.

    Examples
    --------
    *Adding InChI to a Molecule's JSON*

    You want to use the InChI as part of a JSON representation of a
    molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            key_makers=(stk.Inchi(), ),
        )
        # Get the JSON representation, including an InChI.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self):
        """
        Initialize a :class:`Inchi` instance.

        """

        super().__init__('InChI', get_inchi)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'Inchi()'
