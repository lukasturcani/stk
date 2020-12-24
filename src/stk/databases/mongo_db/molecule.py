"""
Molecule MongoDB
================

"""

from functools import lru_cache
import warnings

from stk.serialization import (
    MoleculeJsonizer,
    MoleculeDejsonizer,
)
from ..molecule import MoleculeDatabase
from .utilities import HashableDict


class MoleculeMongoDb(MoleculeDatabase):
    """
    Uses MongoDB to store and retrieve molecules.

    See Also
    --------
    :class:`.ConstructedMoleculeMongoDb`
        You can store :class:`.ConstructedMolecule` instances with \
        :class:`.MoleculeMongoDb`, however you can only retrieve \
        them as :class:`.Molecule` instances. If you want to store \
        and retrieve :class:`.ConstructedMolecule` instances, you \
        should use :class:`.ConstructedMoleculeMongoDb`.

    Examples
    --------
    *Storing and Retrieving Constructed Molecules*

    You want to store and retrieve a molecule from the database

    .. code-block:: python

        import stk
        import pymongo

        # Connect to a MongoDB. This example connects to a local
        # MongoDB, but you can connect to a remote DB too with
        # MongoClient() - read the documentation for pymongo to see how
        # to do that.
        client = pymongo.MongoClient()
        db = stk.MoleculeMongoDb(client)

        # Create a molecule.
        molecule = stk.BuildingBlock('NCCN')

        # Place it into the database.
        db.put(molecule)

        # Retrieve it from the database.
        key_maker = stk.InchiKey()
        retrieved = db.get({
            key_maker.get_key_name(): key_maker.get_key(molecule)
        })

    Note that the molecule retrieved from that database can have
    a different atom ordering than the one put into it. So while the
    molecule will have the same structure, the order of the atoms
    may be different to the molecule placed into the database. This
    is because the database gives the molecules a canonical atom
    ordering, which allows position matrices to be used across
    different atom id orderings.

    *Iterating over All Entries in the Database*

    All entries in a database can be iterated over very simply

    .. code-block:: python

        for entry in db.get_all():
            # Do something to entry.

    *Using Alternative Keys for Retrieving Molecules*

    By default, the only molecular key the database stores, is the
    InChIKey. However, additional keys can be added to the JSON stored
    in the database by using a different :class:`.MoleculeJsonizer`

    .. code-block:: python

        db = stk.MoleculeMongoDb(
            mongo_client=client,
            # Store the InChI and the InChIKey of molecules in
            # the JSON representation.
            jsonizer=stk.MoleculeJsonizer(
                key_makers=(stk.Inchi(), stk.InchiKey()),
            )
        )
        # Places the JSON of the molecule into the database. In this
        # case, the JSON includes both the InChI and the InChIKey.
        db.put(molecule)

        # You can now use the InChI or the InChIKey to retrieve the
        # molecule from the database.
        key_maker = stk.Inchi()
        retrieved = db.get({
            key_maker.get_key_name(): key_maker.get_key(molecule),
        })

    Obviously, most of the time, you won't have the molecule you are
    trying to retrieve from the database. Maybe you only have the
    SMILES of the molecule. You can still retrieve it.

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit

        retrieved = db.get(
            'InChI': rdkit.MolToInchi(rdkit.MolFromSmiles('NCCN'))
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

        db = stk.MoleculeMongoDb(
            mongo_client=client,
            jsonizer=stk.MoleculeJsonizer(
                # Include your own custom key maker in the JSON
                # representation.
                key_makers = (stk.Inchi(), stk.InchiKey(), Smiles()),
            ),
        )

        molecule = stk.BuildingBlock('BrBr')

        # Place the JSON of your molecule into the database. In this
        # case the JSON will include a key called "SMILES" and
        # the value will be the SMILES of the molecule.
        db.put(molecule)

        # You can now find your molecule by using SMILES as the key.
        retrieved = db.get({'SMILES': 'BrBr'})

    Often, it is unnecessary to create a whole subclass for a your
    custom key

    .. code-block:: python

        smiles = stk.MoleculeKeyMaker(
            key_name='SMILES',
            get_key=lambda molecule:
                rdkit.MolToSmiles(molecule.to_rdkit_mol()),
        )
        db = stk.MoleculeMongoDb(
            mongo_client=client,
            jsonizer=stk.MoleculeJsonizer(
                key_makers=(stk.InchiKey(), smiles),
            ),
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
        position_matrix_collection=None,
        jsonizer=MoleculeJsonizer(),
        dejsonizer=MoleculeDejsonizer(),
        lru_cache_size='',
        put_lru_cache_size=128,
        get_lru_cache_size=128,
        indices=('InChIKey', ),
    ):
        """
        Initialize a :class:`.MoleculeMongoDb` instance.

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
            the cache. When ``None``, defaults to
            'building_block_position_matrices`.

        jsonizer : :class:`.MoleculeJsonizer`
            Used to create the JSON representations of molecules
            stored in the database.

        dejsonizer : :class:`.MoleculeDejsonizer`
            Used to create :class:`.Molecule` instances from their
            JSON representations.

        lru_cache_size : :class:`int`, optional
            This argument is deprecated and will be removed in any
            version of :mod:`stk` released on, or after, 01/01/21.
            Use the `put_lru_cache_size` and `get_lru_cache_size`
            arguments instead.

            A RAM-based least recently used cache is used to avoid
            reading and writing to the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        put_lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            writing to the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        get_lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            reading from the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        indices : :class:`tuple` of :class:`str`, optional
            The names of molecule keys, on which an index should be
            created, in order to minimize lookup time.

        """

        if position_matrix_collection is None:
            position_matrix_collection = (
                'building_block_position_matrices'
            )
            warnings.warn(
                "The default value of the position_matrix_collection "
                "parameter in MoleculeMongoDb has changed from "
                "'position_matrices' "
                "to 'building_block_position_matrices'! To remove "
                "this warning, set the position_matrix_collection "
                "parameter explicitly. You can still use the default "
                "value, for example, position_matrix_collection='"
                "building_block_position_matrices' will remove the "
                "warning but use the default value. "
                "This warning will be removed "
                "in any version of stk released on, or after, "
                "1/03/21.",
                UserWarning,
            )

        if lru_cache_size != '':
            warnings.warn(
                'The lru_cache_size argument is deprecated and will '
                'be removed in any version of stk released on, or '
                'after, 01/01/21. Use the put_lru_cache_size and '
                'get_lru_cache_size arguments instead.',
                FutureWarning,
            )
            put_lru_cache_size = lru_cache_size
            get_lru_cache_size = lru_cache_size

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._position_matrices = database[position_matrix_collection]
        self._jsonizer = jsonizer
        self._dejsonizer = dejsonizer

        self._get = lru_cache(maxsize=get_lru_cache_size)(self._get)
        self._put = lru_cache(maxsize=put_lru_cache_size)(self._put)

        for index in indices:
            # Do not create the same index twice.
            if f'{index}_1' not in self._molecules.index_information():
                self._molecules.create_index(index)
            if (
                f'{index}_1'
                not in self._position_matrices.index_information()
            ):
                self._position_matrices.create_index(index)

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
        return self._put(HashableDict(json))

    def _put(self, json):
        keys = dict(json['matrix'])
        keys.pop('m')

        query = {'$or': []}
        for key, value in keys.items():
            query['$or'].append({key: value})

        self._molecules.update_many(
            filter=query,
            update={
                '$set': json['molecule'],
            },
            upsert=True,
        )
        self._position_matrices.update_many(
            filter=query,
            update={
                '$set': json['matrix'],
            },
            upsert=True,
        )

    def get(self, key):
        # lru_cache requires that the parameters to the cached function
        # are hashable objects.
        return self._get(HashableDict(key))

    def _get(self, key):
        """
        Get the molecule with `key` from the database.

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

        json = self._molecules.find_one(key)
        if json is None:
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

        return self._dejsonizer.from_json({
            'molecule': json,
            'matrix': position_matrix,
        })

    def _get_molecule_keys(self, entry):

        # Ignore keys reserved by position matrix collections.
        reserved_keys = ('m', '_id')

        for key, value in entry.items():
            if key not in reserved_keys:
                yield key, value

    def get_all(self):
        for entry in self._position_matrices.find():
            # Do 'or' query over all key value pairs.
            query = {'$or': [
                {key: value}
                for key, value in self._get_molecule_keys(entry)
            ]}

            json = self._molecules.find_one(query)
            if json is None:
                raise KeyError(
                    'No molecule found in the database associated '
                    f'with a position matrix with query: {query}. '
                    'This suggests your database is corrupted.'
                )

            yield self._dejsonizer.from_json({
                'molecule': json,
                'matrix': entry,
            })
