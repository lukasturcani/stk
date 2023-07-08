from stk._internal.key_makers.molecule import MoleculeKeyMaker
from stk._internal.key_makers.utilities import get_inchi
from stk._internal.molecule import Molecule


class Inchi(MoleculeKeyMaker):
    """
    Used to get the InChI of molecules.

    Examples:

        *Adding InChI to a Molecule's JSON*

        You want to use the InChI as part of a JSON representation of a
        molecule

        .. testcode:: adding-inchi-to-a-molecules-json

            import stk

            jsonizer = stk.MoleculeJsonizer(
                key_makers=(stk.Inchi(), ),
            )
            building_block = stk.BuildingBlock('NCCN')
            # Get the JSON representation, including an InChI.
            json = jsonizer.to_json(building_block)

        .. testcode:: adding-inchi-to-a-molecules-json
            :hide:

            _inchi = stk.Inchi()
            _expected_inchi = 'InChI=1S/C2H8N2/c3-1-2-4/h1-4H2'
            _molecule_json = json['molecule']
            assert (
                _molecule_json[_inchi.get_key_name()]
                == _expected_inchi
            )

            _matrix_json = json['matrix']
            assert (
                _matrix_json[_inchi.get_key_name()]
                == _expected_inchi
            )

    """

    def __init__(self) -> None:
        return

    def get_key_name(self) -> str:
        return "InChI"

    def get_key(self, molecule: Molecule) -> str:
        return get_inchi(molecule)

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "Inchi()"
