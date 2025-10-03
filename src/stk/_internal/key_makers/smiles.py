from stk._internal.molecule import Molecule

from .molecule import MoleculeKeyMaker
from .utilities import get_smiles


class Smiles(MoleculeKeyMaker):
    """
    Used to get the SMILES of molecules.

    .. warning::

        SMILES strings generated with ``stk`` version ``v2025.07.17.0`` or sooner
        (using ``rdkit`` version ``2024.9.1`` or sooner) will be different than newer
        versions due to a change in handling the valence of organic atoms bound to
        metals. Details can be found in the ``rdkit`` release notes `2025_03_1`_.
        No changes occur in ``stk`` construction. An example change in SMILES:
        ``CCCO->[Fe+2]`` becomes ``CCC[OH]->[Fe+2]``.

    .. _`2025_03_1`: https://github.com/rdkit/rdkit/releases/tag/Release_2025_03_1

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
        return

    def get_key_name(self) -> str:
        return "SMILES"

    def get_key(self, molecule: Molecule) -> str:
        return get_smiles(molecule)

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "Smiles()"
