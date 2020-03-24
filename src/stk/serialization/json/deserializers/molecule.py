"""
Molecule DeJSONizer
===================

"""

from stk.molecular import Molecule

from .utilities import to_atom, to_bond


class MoleculeDejsonizer:
    """

    """

    def from_json(self, json, position_matrix):
        """

        """

        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(json['atoms'])
        )
        return Molecule(
            atoms=atoms,
            bonds=tuple(
                to_bond(atoms, bond_json)
                for bond_json in json['bonds']
            ),
            position_matrix=position_matrix,
        )
