"""
MongoDB Constructed Molecule Cache
==================================

"""

import numpy as np
from stk.molecular import InchiKey
from stk.serialization import (
    MoleculeJsonizer,
    ConstructedMoleculeJsonizer,
    ConstructedMoleculeDejsonizer,
)
from .constructed_molecule import ConstructedMoleculeCache


class MongoDbConstructedMoleculeCache(ConstructedMoleculeCache):
    """

    """

    def __init__(
        self,
        mongo_client,
        database='stk',
        molecule_collection='molecules',
        constructed_molecule_collection='constructed_molecules',
        position_matrix_collection='position_matrices',
        key_makers=(InchiKey(), ),
        molecule_jsonizer=MoleculeJsonizer(),
        constructed_molecule_jsonizer=ConstructedMoleculeJsonizer(),
        molecule_dejsonizer=ConstructedMoleculeDejsonizer(),
    ):
        """

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._constructed_molecules = database[
            constructed_molecule_collection
        ]
        self._position_matrices = database[position_matrix_collection]
        self._key_makers = key_makers
        self._molecule_jsonizer = molecule_jsonizer
        self._constructed_molecule_jsonizer = (
            constructed_molecule_jsonizer
        )
        self._molecule_dejsonizer = molecule_dejsonizer

    def put(self, molecule):
        molecule = molecule.with_canonical_atom_ordering()

        position_matrix_json = {
            'position_matrix':
                molecule.get_position_matrix().tolist()
        }
        for key_maker in self._key_makers:
            position_matrix_json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        self._position_matrices.insert_one(position_matrix_json)

        self._molecules.insert_one(
            document=self._molecule_jsonizer.to_json(molecule),
        )
        self._constructed_molecules.insert_one(
            document=self._constructed_molecule_jsonizer.to_json(
                molecule=molecule,
            )
        )

    def get(self, key, building_blocks, default=None):
        molecule_json = self._molecules.find_one(key)
        if molecule_json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif molecule_json is None:
            return default

        constructed_molecule_json = (
            self._constructed_molecules.find_one(key)
        )
        if constructed_molecule_json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif constructed_molecule_json is None:
            return default

        position_matrix = np.array(
            self._position_matrices.find_one(key)['position_matrix'],
            dtype=np.float64,
        )
        return self._molecule_dejsonizer.from_json(
                molecule_json=molecule_json,
                constructed_molecule_json=constructed_molecule_json,
                position_matrix=position_matrix,
                building_blocks=building_blocks,
            )
