"""
MongoDB Molecular Cache
=======================

"""

import numpy as np
from stk.molecular import ConstructedMolecule, InchiKey
from stk.serialization import (
    MoleculeJsonizer,
    ConstructedMoleculeJsonizer,
    MoleculeDejsonizer,
    ConstructedMoleculeDejsonizer,
)
from .molecular_cache import MolecularCache


class MongoDbMolecularCache(MolecularCache):
    """
    Uses MongoDB to store and retrieve molecules.

    Examples
    --------

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
        molecule_dejsonizer=MoleculeDejsonizer(),
        constructed_molecule_dejsonizer=(
            ConstructedMoleculeDejsonizer()
        ),
        building_block_cache=None,
    ):
        """
        Initialize a :class:`.MongoDbMolecularCache` instance.

        Parameters
        ----------
        mongo_client : :class:`pymongo.MongoClient`
            The database client.

        database : :class:`str`
            The name of the database to use.

        molecule_collection : :class:`str`
            The name of the collection which stores molecular
            information.

        constructed_molecule_collection : :class:`str`
            The name of the collection which stores additional
            constructed molecule information.

        molecule_jsonizer : :class:`.MoleculeJsonizer`
            Used to create the JSON representations of molecules
            stored in the database.

        constructed_molecule_jsonizer : \
                :class:`.ConstructedMoleculeJsonizer`
            Used to create the JSON representations of constructed
            molecule data stored in the database.

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._constructed_molecules = database[
            constructed_molecule_collection
        ]
        self._position_matrices = database[position_matrix_collection]
        self._molecule_jsonizer = molecule_jsonizer
        self._constructed_molecule_jsonizer = (
            constructed_molecule_jsonizer
        )
        self._molecule_dejsonizer = molecule_dejsonizer
        self._constructed_molecule_dejsonizer = (
            constructed_molecule_dejsonizer
        )
        if building_block_cache is None:
            building_block_cache = self
        self._building_block_cache = building_block_cache

    def put(self, molecule):
        molecule = molecule.with_canonical_atom_ordering()

        position_matrix_json = {
            'position_matrix':
                molecule.get_position_matrix().tolist(),
        }
        for key_maker in self._key_makers:
            position_matrix_json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        self._position_matrices.insert_one(position_matrix_json)

        json = self._molecule_jsonizer.to_json(molecule)
        self._molecules.insert_one(json)

        if not isinstance(molecule, ConstructedMolecule):
            return

        json = self._constructed_molecule_jsonizer.to_json(molecule)
        self._constructed_molecules.insert_one(json)

    def get(self, key, default=None):
        molecule_json = self._molecules.find_one(key)
        if molecule_json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif molecule_json is None:
            return default

        position_matrix = np.array(
            self._position_matrices.find_one(key),
            dtype=np.float64,
        )

        constructed_molecule_json = (
            self._constructed_molecules.find_one(key)
        )
        if constructed_molecule_json is None:
            return self._molecule_dejsonizer.from_json(
                json=molecule_json,
                position_matrix=position_matrix,
            )

        building_block_jsons = self._molecules.find(
            {'$or': constructed_molecule_json['BB']}
        )
        building_block_matrices = self._position_matrices.find(

        )

        return self._constructed_molecule_dejsonizer.from_json(
            molecule_json=molecule_json,
            building_block_jsons=building_block_jsons,
            constructed_molecule_json=constructed_molecule_json
        )
