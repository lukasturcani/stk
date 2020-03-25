"""
MongoDB Molecule Cache
======================

"""

import numpy as np
from stk.molecular import InchiKey
from stk.serialization import (
    MoleculeJsonizer,
    MoleculeDejsonizer,
)
from . molecule import MoleculeCache


class MongoDbMoleculeCache(MoleculeCache):
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
        position_matrix_collection='position_matrices',
        key_makers=(InchiKey(), ),
        molecule_jsonizer=MoleculeJsonizer(),
        molecule_dejsonizer=MoleculeDejsonizer(),
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

        position_matrix_collection : :class:`str`
            The name of the collection which stores the position
            matrices of the molecules put into and retrieved from
            the cache.

        key_makers : :class:`tuple` of :class:`.MoleculeKeyMaker`
            Used to make the keys for molecules, which are used to
            reference the position matrices in the
            `position_matrix_collection`.

        molecule_jsonizer : :class:`.MoleculeJsonizer`
            Used to create the JSON representations of molecules
            stored in the database.

        molecule_dejsonizer : :class:`.MoleculeDejsonizer`
            Used to create :class:`.Molecule` instances from their
            JSON representations.

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._position_matrices = database[position_matrix_collection]
        self._key_makers = key_makers
        self._molecule_jsonizer = molecule_jsonizer
        self._molecule_dejsonizer = molecule_dejsonizer

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
            self._position_matrices.find_one(key)['position_matrix'],
            dtype=np.float64,
        )
        return self._molecule_dejsonizer.from_json(
                json=molecule_json,
                position_matrix=position_matrix,
            )
