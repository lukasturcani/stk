"""
SMILES
======

"""

from . import molecule as _molecule
from . import utilities as _utilities


__all__ = (
    'Smiles',
)


class Smiles(_molecule.MoleculeKeyMaker):
    """
    Used to get the SMILES of molecules.

    Examples:

        *Adding SMILES to a Molecule's JSON*

        You want to use the isomeric, canonical SMILES from RDKit as
        part of a JSON representation of a molecule

        .. testcode:: adding-smiles-to-a-molecules-json

            import stk

            jsonizer = stk.MoleculeJsonizer(
                key_makers=(stk.Smiles(), ),
            )
            # Get the JSON representation, including an SMILES.
            json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

        .. testcode:: adding-smiles-to-a-molecules-json
            :hide:

            assert json['molecule']['SMILES'] == 'NCCN'
            assert json['matrix']['SMILES'] == 'NCCN'

    """

    def __init__(self) -> None:
        """
        Initialize a :class:`.Smiles` instance.

        """

        _molecule.MoleculeKeyMaker.__init__(
            self=self,
            key_name='SMILES',
            get_key=_utilities.get_smiles,
        )

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return 'Smiles()'
