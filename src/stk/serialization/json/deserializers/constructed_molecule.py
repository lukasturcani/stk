"""
Constructed Molecule DeJSONizer
===============================

"""

from stk.molecular import ConstructedMolecule

from .molecule import MoleculeDejsonizer
from .utilities import to_atom, to_bond, to_atom_info, to_bond_info


class ConstructedMoleculeDejsonizer:
    """

    """

    def from_json(
        self,
        molecule_json,
        building_block_jsons,
        constructed_molecule_json,
        building_block_position_matrices,
        position_matrix,
    ):

        dejsonizer = MoleculeDejsonizer()
        building_blocks = tuple(
            dejsonizer.from_json(json, position_matrix)
            for json, position_matrix
            in zip(
                building_block_jsons,
                building_block_position_matrices,
            )
        )
        num_building_blocks = (
            (building_block, num)
            for building_block, num in zip(
                building_blocks,
                constructed_molecule_json['num_building_blocks'],
            )
        )
        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json
            in enumerate(constructed_molecule_json['atoms'])
        )
        bonds = tuple(
            to_bond(atoms, bond_json)
            for bond_json in constructed_molecule_json['bonds']
        )
        return ConstructedMolecule.init(
            atoms=atoms,
            bonds=bonds,
            position_matrix=position_matrix,
            atom_infos=tuple(
                to_atom_info(
                    building_blocks=building_blocks,
                    atom=atoms[atom_id],
                    json=atom_info_json,
                )
                for atom_id, atom_info_json
                in enumerate(constructed_molecule_json['atom_infos'])
            ),
            bond_infos=tuple(
                to_bond_info(
                    building_blocks=building_blocks,
                    bond=bonds[bond_id],
                    json=bond_info_json,
                )
                for bond_id, bond_info_json
                in enumerate(constructed_molecule_json['bond_infos'])
            ),
            num_building_blocks=num_building_blocks,
        )
