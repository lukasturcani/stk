"""
MongoDB Molecular Cache
=======================

"""

from stk.molecular import ConstructedMolecule
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
        molecule_jsonizer=MoleculeJsonizer(),
        constructed_molecule_jsonizer=ConstructedMoleculeJsonizer(),
        molecule_dejsonizer=MoleculeDejsonizer(),
        constructed_molecule_dejsonizer=(
            ConstructedMoleculeDejsonizer()
        ),
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
        self._molecule_collection = database[molecule_collection]
        self._constructed_molecule_collection = database[
            constructed_molecule_collection
        ]
        self._molecule_jsonizer = molecule_jsonizer
        self._constructed_molecule_jsonizer = (
            constructed_molecule_jsonizer
        )
        self._molecule_dejsonizer = molecule_dejsonizer
        self._constructed_molecule_dejsonizer = (
            constructed_molecule_dejsonizer
        )

    def put(self, molecule):
        json = self._molecule_jsonizer.to_json(molecule)
        self._molecule_collection.insert_one(json)

        if not isinstance(molecule, ConstructedMolecule):
            return

        json = self._constructed_molecule_jsonizer.to_json(molecule)
        self._constructed_molecule_collection.insert_one(json)

    def get(self, key, default=None):
        molecule_json = self._molecule_collection.find_one(key)
        if molecule_json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif molecule_json is None:
            return default

        constructed_molecule_json = (
            self._constructed_molecule_collection.find_one(key)
        )
        if constructed_molecule_json is None:
            return self._molecule_dejsonizer.from_json(molecule_json)

        building_blocks = constructed_molecule_json['building_blocks']
        building_block_jsons = self._molecule_collection.find(
            {'$or': building_blocks}
        )
        return self._constructed_molecule_dejsonizer.from_json(
            molecule_json=molecule_json,
            building_block_jsons=building_block_jsons,
            constructed_molecule_json=constructed_molecule_json
        )
