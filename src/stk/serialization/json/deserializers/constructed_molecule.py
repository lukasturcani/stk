"""
Constructed Molecule DeJSONizer
===============================

"""

import numpy as np

from stk.molecular import ConstructedMolecule

from .molecule import MoleculeDejsonizer
from .utilities import to_atom, to_atom_info, to_bond, to_bond_info


class ConstructedMoleculeDejsonizer:
    """
    Abstract base class for creating constructed molecules from JSONs.

    See Also
    --------
    :class:`.MoleculeDejsonizer`

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. These are just default implementations,
    which can be safely ignored or overridden, when implementing
    subclasses. However, the default implementation can be used
    directly, if it suits your needs.

    """

    def __init__(self):
        """
        Initialize a :class:`.ConstructedMoleculeDejsonizer` instance.

        """

        self._dejsonizer = MoleculeDejsonizer()

    def from_json(self, json):
        """
        Get a :class:`.ConstructedMolecule` from a JSON.

        Parameters
        ----------
        json : :class:`dict`
            A JSON of the constructed molecule.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The constructed molecule.

        """

        building_blocks = tuple(map(
            self._dejsonizer.from_json,
            json['buildingBlocks'],
        ))

        num_building_blocks = (
            (building_block, num)
            for building_block, num in zip(
                building_blocks,
                json['constructedMolecule']['nBB'],
            )
        )
        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(json['molecule']['a'])
        )
        bonds = tuple(
            to_bond(atoms, bond_json)
            for bond_json in json['molecule']['b']
        )
        return ConstructedMolecule.init(
            atoms=atoms,
            bonds=bonds,
            position_matrix=np.array(json['matrix']['m']),
            atom_infos=tuple(
                to_atom_info(
                    building_blocks=building_blocks,
                    atom=atoms[atom_id],
                    json=atom_info_json,
                )
                for atom_id, atom_info_json
                in enumerate(json['constructedMolecule']['aI'])
            ),
            bond_infos=tuple(
                to_bond_info(
                    building_blocks=building_blocks,
                    bond=bonds[bond_id],
                    json=bond_info_json,
                )
                for bond_id, bond_info_json
                in enumerate(json['constructedMolecule']['bI'])
            ),
            num_building_blocks=num_building_blocks,
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}()'
