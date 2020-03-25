"""
MongoDB Molecule Cache
======================

"""

from stk.serialization import (
    MoleculeJsonizer,
    MoleculeDejsonizer,
)
from . molecule import MoleculeCache


class MongoDbMoleculeCache(MoleculeCache):
    """
    Uses MongoDB to store and retrieve molecules.

    See Also
    --------
    :class:`.MongoDbConstructedMoleculeCache`
        You can store :class:`.ConstructedMolecule` instances with \
        :class:`.MongoDbMoleculeCache`, however you can only retrieve \
        them as :class:`.Molecule` instances. If you want to store \
        and retrieve :class:`.ConstructedMolecule` instances, you \
        should use :class:`.MongoDbConstructedMoleculeCache`.

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

    By default, the only molecular key the database stores, is the
    InChIKey. However, additional keys can be added to the JSON stored
    in the database by using a different :class:`.MoleculeJsonizer`

    .. code-block:: python

        db = stk.MongoDbMoleculeCache(
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

        key_maker = stk.Inchi()
        retrieved = db.get(
            key_maker.get_key_name():
                rdkit.MolToInchi(rdkit.MolFromSmiles('NCCN'))
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

        db = stk.MongoDbMoleculeCache(
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
        db = stk.MongoDbMoleculeCache(
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
        position_matrix_collection='position_matrices',
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
        self._jsonizer = jsonizer
        self._dejsonizer = dejsonizer

    def put(self, molecule):
        molecule = molecule.with_canonical_atom_ordering()
        json = self._jsonizer.to_json(molecule)
        self._molecules.insert_one(json['molecule'])
        self._position_matrices.insert_one(json['matrix'])

    def get(self, key, default=None):
        json = self._molecules.find_one(key)
        if json is None and default is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        elif json is None:
            return default

        return self._dejsonizer.from_json({
            'molecule': json,
            'matrix': self._position_matrices.find_one(key),
        })
