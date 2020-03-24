"""
Constructed Molecule DeJSONizer
===============================

"""

from stk.molecular import ConstructedMolecule

from .utilities import to_atom, to_bond, to_atom_info, to_bond_info


class ConstructedMoleculeDejsonizer:
    """
    Creates :class:`.ConstructedMolecule` instances from JSONs.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.ConstructedMoleculeDejsonizer` instance.

        """

        return

    def from_json(
        self,
        molecule_json,
        constructed_molecule_json,
        position_matrix,
        building_blocks,
    ):
        """
        Get a :class:`.ConstructedMolecule` from a JSON.

        Parameters
        ----------
        molecule_json : :class:`dict`
            A JSON of the molecular information of the constructed
            molecule.

        constructed_molecule_json : :class:`dict`
            A JSON of the constructed molecule information of the
            constructed molecule.

        position_matrix : :class:`numpy.ndarray`
            The position matrix of the constructed molecule.

        building_blocks : :class:`tuple` of :class:`.Molecule`
            The building blocks of the constructed molecule.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The constructed molecule.

        """

        num_building_blocks = (
            (building_block, num)
            for building_block, num in zip(
                building_blocks,
                constructed_molecule_json['nBB'],
            )
        )
        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(molecule_json['a'])
        )
        bonds = tuple(
            to_bond(atoms, bond_json)
            for bond_json in molecule_json['b']
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
                in enumerate(constructed_molecule_json['aI'])
            ),
            bond_infos=tuple(
                to_bond_info(
                    building_blocks=building_blocks,
                    bond=bonds[bond_id],
                    json=bond_info_json,
                )
                for bond_id, bond_info_json
                in enumerate(constructed_molecule_json['bI'])
            ),
            num_building_blocks=num_building_blocks,
        )