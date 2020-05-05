"""
SMILES
======

"""

from .molecule import MoleculeKeyMaker
from .utilities import get_smiles


class Smiles(MoleculeKeyMaker):
    """
    Used to get the SMILES of molecules.

    Examples
    --------
    *Adding SMILES to a Molecule's JSON*

    You want to use the isomeric, canonical SMILES from RDKit as part
    of a JSON representation of a molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            key_makers=(stk.Smiles(), ),
        )
        # Get the JSON representation, including an SMILES.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self):
        """
        Initialize a :class:`.Smiles` instance.

        """

        super().__init__('SMILES', get_smiles)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'Smiles()'
