"""
Constructed Molecule MongoDB
============================

"""

from functools import lru_cache

from stk.serialization import (
    ConstructedMoleculeJsonizer,
    ConstructedMoleculeDejsonizer,
)
from ..constructed_molecule import ConstructedMoleculeDatabase
from .utilities import HashableDict


class ConstructedMoleculeMongoDb(ConstructedMoleculeDatabase):
    """
    Uses MongoDB to store and retrieve constructed molecules.

    See Also
    --------
    :class:`.MoleculeMongoDb`
        If you need to store and retrieve molecules, which are not
        :class:`.ConstructedMolecule` instances, use a \
        :class:`.MoleculeMongoDb`.

    Examples
    --------
    Before using this class, make sure you have :mod:`pymongo` and
    that its working properly. I recommend reading at least the
    introductory and installation
    documentation of :mod:`pymongo` before using this class. Those
    docs can be found here__.

    __ https://api.mongodb.com/python/current/

    *Usage*

    You want to store and retrieve a :class:`.ConstructedMolecule`
    from the database

    .. code-block:: python

        import stk
        # pymongo does not come with stk, you have to install it
        # explicitly with "pip install pymongo".
        import pymongo

        # Connect to a MongoDB. This example connects to a local
        # MongoDB, but you can connect to a remote DB too with
        # MongoClient() - read the documentation for pymongo to see how
        # to do that.
        client = pymongo.MongoClient()
        db = stk.ConstructedMoleculeMongoDb(client)

        # Create a molecule.
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )

        # Place it into the database.
        db.put(polymer)

        # Retrieve it from the database.
        key_maker = stk.InchiKey()
        retrieved = db.get({
            key_maker.get_key_name(): key_maker.get_key(polymer),
        })

    Note that the molecule retrieved from that database can have
    a different atom ordering than the one put into it. So while the
    molecule will have the same structure, the order of the atoms
    may be different to the molecule placed into the database. This
    is because the database gives the molecules a canonical atom
    ordering, which allows position matrices to be used across
    different atom id orderings.

    By default, the only molecular key the database stores, is the
    InChIKey. However, additional keys can be added to the JSON stored
    in the database by using a different
    :class:`.ConstructedMoleculeJsonizer`

    .. code-block:: python

        db = stk.ConstructedMoleculeMongoDb(
            mongo_client=client,
            # Store the InChI and the InChIKey of molecules in
            # the JSON representation.
            jsonizer=stk.ConstructedMoleculeJsonizer(
                key_makers=(stk.Inchi(), stk.InchiKey()),
            ),
        )
        # Places the JSON of the molecule into the database. In this
        # case, the JSON includes both the InChI and the InChIKey.
        db.put(polymer)

        # You can now use the InChI or the InChIKey to retrieve the
        # molecule from the database.
        key_maker = stk.Inchi()
        retrieved = db.get({
            key_maker.get_key_name(): key_maker.get_key(polymer),
        })

    Obviously, most of the time, you won't have the molecule you are
    trying to retrieve from the database. Maybe you only have the
    SMILES of the molecule. You can still retrieve it.

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit

        retrieved = db.get(
            'InChI': rdkit.MolToInchi(rdkit.MolFromSmiles('BrCCCCBr'))
        )

    As long as you have the name of the key, and the expected value
    of the key, you can retrieve your molecule from the database.

    Note that you can create your own keys and add them to the database

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit

        # Create your own key. This one is called "SMILES" and the
        # value is the SMILES of the molecule.
        class Smiles(stk.MoleculeKeyMaker):
            def __init__(self):
                return

            def get_key_name(self):
                return 'SMILES'

            def get_key(self, molecule):
                return rdkit.MolToSmiles(molecule.to_rdkit_mol())

        db = stk.ConstructedMoleculeMongoDb(
            mongo_client=client,
            jsonizer=stk.ConstructedMoleculeJsonizer(
                # Include your own custom key maker in the JSON
                # representation.
                key_makers = (stk.Inchi(), stk.InchiKey(), Smiles()),
            ),
        )

        # Place the JSON of your molecule into the database. In this
        # case the JSON will include a key called "SMILES" and
        # the value will be the SMILES of the molecule.
        db.put(polymer)

        # You can now find your molecule by using SMILES as the key.
        retrieved = db.get({'SMILES': 'BrCCCCBr'})

    Often, it is unnecessary to create a whole subclass for a your
    custom key

    .. code-block:: python

        smiles = stk.MoleculeKeyMaker(
            key_name='SMILES',
            get_key=lambda molecule:
                rdkit.MolToSmiles(molecule.to_rdkit_mol()),
        )
        db = stk.ConstructedMoleculeMongoDb(
            mongo_client=client,
            jsonizer=stk.ConstructedMoleculeJsonizer(
            key_makers=(stk.InchiKey(), smiles),
        )

    Note that the key you use to get the molecule back from the
    database should be unique. In other words, there should always just
    be one molecule which has that key in the database. Using a key
    that is matched by multiple molecules will likely cause your data
    to be jumbled. For example, you might return the atoms of one
    of the molecules matched by the key but holding the position matrix
    of the second molecule matched by the key.

    """

    def __init__(
        self,
        mongo_client,
        database='stk',
        molecule_collection='molecules',
        constructed_molecule_collection='constructed_molecules',
        position_matrix_collection='position_matrices',
        jsonizer=ConstructedMoleculeJsonizer(),
        dejsonizer=ConstructedMoleculeDejsonizer(),
        lru_cache_size=128,
    ):
        """
        Initialize a :class:`.ConstructedMoleculeMongoDb`.

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
            The name of the collection which stored constructed
            molecule information, that does not belong in the
            `molecule_collection`.

        position_matrix_collection : :class:`str`
            The name of the collection which stores the position
            matrices of the molecules put into and retrieved from
            the database.

        jsonizer : :class:`.ConstructedMoleculeJsonizer`
            Used to create the JSON representations of molecules
            stored in the database.

        dejsonizer : :class:`.ConstructedMoleculeDejsonizer`
            Used to create :class:`.Molecule` instances from their
            JSON representations.

        lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            reading and writing to the database repeatedly. This sets
            the number of molecules which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._constructed_molecules = database[
            constructed_molecule_collection
        ]
        self._position_matrices = database[position_matrix_collection]
        self._jsonizer = jsonizer
        self._dejsonizer = dejsonizer

        self._get = lru_cache(maxsize=lru_cache_size)(self._get)
        self._put = lru_cache(maxsize=lru_cache_size)(self._put)

    def put(self, molecule):
        molecule = molecule.with_canonical_atom_ordering()
        json = self._jsonizer.to_json(molecule)
        # lru_cache requires that the parameters to the cached function
        # are hashable objects.
        json['matrix']['m'] = tuple(
            tuple(row) for row in json['matrix']['m']
        )
        json['matrix'] = HashableDict(json['matrix'])
        json['molecule'] = HashableDict(json['molecule'])
        json['constructedMolecule'] = HashableDict(
            json['constructedMolecule']
        )
        json['constructedMolecule']['BB'] = tuple(map(
            HashableDict,
            json['constructedMolecule']['BB'],
        ))

        def make_hashable(json):
            json['matrix']['m'] = tuple(
                tuple(row) for row in json['matrix']['m']
            )
            json['matrix'] = HashableDict(json['matrix'])
            json['molecule'] = HashableDict(json['molecule'])
            return HashableDict(json)

        json['buildingBlocks'] = tuple(map(
            make_hashable,
            json['buildingBlocks'],
        ))
        return self._put(HashableDict(json))

    def _put(self, json):
        self._position_matrices.insert_one(json['matrix'])
        self._molecules.insert_one(json['molecule'])
        self._constructed_molecules.insert_one(
            document=json['constructedMolecule'],
        )
        for building_block_json in json['buildingBlocks']:
            self._molecules.insert_one(building_block_json['molecule'])
            self._position_matrices.insert_one(
                document=building_block_json['matrix'],
            )

    def get(self, key):
        # lru_cache requires that the parameters to the cached function
        # are hashable objects.
        return self._get(HashableDict(key))

    def _get(self, key):
        """
        Get the molecule with `key` from the cache.

        Parameters
        ----------
        key : :class:`.HashableDict`
            The key of a molecule, which is to be returned from the
            database.

        Returns
        -------
        :class:`.Molecule`
            The molecule held in the database under `key`.

        """

        molecule_json = self._molecules.find_one(key)
        if molecule_json is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )

        constructed_molecule_json = (
            self._constructed_molecules.find_one(key)
        )
        if constructed_molecule_json is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )

        position_matrix = self._position_matrices.find_one(key)
        if position_matrix is None:
            raise KeyError(
                'No position matrix found in the database with a key '
                f'of: {key}'
            )

        return self._dejsonizer.from_json(
            json={
                'molecule': molecule_json,
                'constructedMolecule': constructed_molecule_json,
                'matrix': position_matrix,
                'buildingBlocks': tuple(map(
                    self._get_building_block,
                    constructed_molecule_json['BB'],
                ))
            },
        )

    def _get_building_block(self, key):
        return {
            'molecule': self._molecules.find_one(key),
            'matrix': self._position_matrices.find_one(key),
        }
