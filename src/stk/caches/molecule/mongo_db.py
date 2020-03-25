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
    Before using this class, make sure you have :mod:`pymongo` and
    that its working properly. I recommend reading at least the
    introductory and installation
    documentation of :mod:`pymongo` before using this class. Those
    docs can be found here__.

    __ https://api.mongodb.com/python/current/

    *Usage*

    You want to store and retrieve a molecule from the database

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
        db = stk.MongoDbMoleculeCache(client)

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

    By default, the database stores the InChIKey of molecules only.
    However, additional keys can be added to the JSON stored in the
    database by using a different :class:`.MoleculeJsonizer`

    .. code-block:: python

        db = stk.MongoDbMoleculeCache(
            mongo_client=client,
            jsonizer=stk.MoleculeJsonizer(
                # Store the InChI and the InChI key of molecules in
                # their JSON representation.
                key_makers=(stk.Inchi(), stk.InchiKey()),
            )
        )
        # Places the JSON of the molecule into the database. In this
        # case, the JSON includes both the InChI and the InChIKey.
        db.put(molecule)

        # You can now use the InChI or the InChI key to retrieve the
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

        key_maker = stk.Inchi()
        retrieved = db.get(
            key_maker.get_key_name():
                rdkit.MolToInchi(rdkit.MolFromSmiles('NCCN'))
        )

    As long as you have the name of the key, and the expected value
    of the key you can retrieve your molecule from the database.

    Note that you can create your own keys and add them to the database

    .. code-block:: python

        # Create your own key. This one is called NumAtoms and its
        # value is the number of atoms in the molecule.
        class NumAtoms(stk.MoleculeKeyMaker):
            def __init__(self):
                return

            def get_key_name(self):
                return 'NumAtoms'

            def get_key(self, molecule):
                return molecule.get_num_atoms()

        db = stk.MongoDbMoleculeCache(
            mongo_client=client,
            jsonizer=stk.MoleculeJsonizer(
                # Include your own custom key maker in the JSON
                # representation.
                key_makers=(stk.Inchi(), stk.InchiKey(), NumAtoms())
            ),
        )

        molecule = stk.BuildingBlock('BrBr')

        # Place the JSON of your molecule into the database. In this
        # case the JSON will include a key called "NumAtoms" and
        # the value will be the number of atoms in the molecule.
        db.put(molecule)

        # You can now find your molecule by putting in the number of
        # atoms.
        retrieved = db.get({'NumAtoms': 2})

    Often, it is unnecessary to create a whole subclass for a your
    custom key

    .. code-block:: python

        num_bonds = stk.MoleculeKey(
            key_name='NumBonds',
            get_key=lambda molecule: molecule.get_num_bonds(),
        )

        db = stk.MongoDbMoleculeCache(
            mongo_client=client,
            jsonizer=stk.MoleculeJsonizer(
                # Include a "NumAtoms" key and a "NumBonds" key.
                key_makers=(stk.InchiKey(), NumAtoms(), num_bonds),
            ),
        )

        # Place a JSON of your molecule into the database. In this
        # example, the JSON will include the keys "NumAtoms" and
        # "NumBond", and their respective values for your molecule.
        db.put(molecule)

        # You can now find your molecule by putting in the number
        # of bonds
        retrieved = db.get({'NumBonds': 1})

        # Or the number of bonds and atoms.
        retrieved = db.get({'NumAtoms': 2, 'NumBonds': 1})

    """

    def __init__(
        self,
        mongo_client,
        database='stk',
        molecule_collection='molecules',
        position_matrix_collection='position_matrices',
        key_makers=(InchiKey(), ),
        jsonizer=MoleculeJsonizer(),
        dejsonizer=MoleculeDejsonizer(),
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

        jsonizer : :class:`.MoleculeJsonizer`
            Used to create the JSON representations of molecules
            stored in the database.

        dejsonizer : :class:`.MoleculeDejsonizer`
            Used to create :class:`.Molecule` instances from their
            JSON representations.

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._position_matrices = database[position_matrix_collection]
        self._key_makers = key_makers
        self._jsonizer = jsonizer
        self._dejsonizer = dejsonizer

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

        json = self._jsonizer.to_json(molecule)
        self._molecules.insert_one(json)

    def get(self, key, default=None):
        json = self._molecules.find_one(key)
        if json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif json is None:
            return default

        position_matrix = np.array(
            self._position_matrices.find_one(key)['position_matrix'],
            dtype=np.float64,
        )
        return self._dejsonizer.from_json(json, position_matrix)
