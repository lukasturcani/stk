"""
MongoDB Constructed Molecule Cache
==================================

"""

from stk.serialization import (
    ConstructedMoleculeJsonizer,
    ConstructedMoleculeDejsonizer,
)
from .constructed_molecule import ConstructedMoleculeCache


class MongoDbConstructedMoleculeCache(ConstructedMoleculeCache):
    """
    Uses MongoDB to store and retrieve constructed molecules.

    The database for constructed molecules makes a separate collection
    for the molecular information, constructed molecule information
    and position matrices. These pieces of data reference each other
    across these collections through the keys provided by
    :class:`.MoleculeKeyMaker` and
    :class:`.ConstructedMoleculeKeyMaker`. The
    :class:`.MoleculeKeyMaker` is used to create keys, which reference
    the molecular data and :class:`.ConstructedMoleculeKeyMaker` is
    used to create keys which reference constructed molecule data.
    Constructed molecule data is the data held by a
    :class:`.ConstructedMolecule`, which is not relevant to a plain
    :class:`.Molecule` instance.

    See Also
    --------
    :class:`.MongoDbMoleculeCache`
        If you need to store and retrieve molecules, which are not
        :class:`.ConstructedMolecule` instances, use a \
        :class:`.MongoDbMoleculeCache`.

    Examples
    --------
    Before using this class, make sure you have :mod:`pymongo` and
    that its working properly. I recommend reading at least the
    introductory and installation
    documentation of :mod:`pymongo` before using this class. Those
    docs can be found here__.

    __ https://api.mongodb.com/python/current/

    This class has a lot in common with
    :class:`.MongoDbMoleculeCache`, so make sure you read that
    documentation first.

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
        db = stk.MongoDbConstructedMoleculeCache(client)

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
        retrieved = db.get(
            key={
                key_maker.get_key_name(): key_maker.get_key(polymer),
            },
            building_blocks=tuple(polymer.get_building_blocks()),
        )

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
    and :class:`.ConstructedMoleculeJsonizer`

    .. code-block:: python


        # Store the InChI and the InChIKey of molecules in
        # the JSON representation.
        key_makers = (stk.Inchi(), stk.InchiKey())
        db = stk.MongoDbConstructedMoleculeCache(
            mongo_client=client,
            molecule_jsonizer=stk.MoleculeJsonizer(key_makers),
            constructed_molecule_jsonizer=(
                stk.ConstructedMoleculeJsonizer(key_makers),
            ),
            # Make sure the position matrices stored in the database
            # can be referenced by the same keys as the molecules.
            key_makers=key_makers,
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
                rdkit.MolToInchi(rdkit.MolFromSmiles('BrCCCCBr'))
        )

    As long as you have the name of the key, and the expected value
    of the key, you can retrieve your molecule from the database.

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

        # Include your own custom key maker in the JSON representation.
        key_makers = (stk.Inchi(), stk.InchiKey(), NumAtoms())
        db = stk.MongoDbConstructedMoleculeCache(
            mongo_client=client,
            molecule_jsonizer=stk.MoleculeJsonizer(key_makers),
            constructed_molecule_jsonizer=(
                stk.ConstructedMoleculeJsonizer(key_makers)
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

        num_bonds = stk.MoleculeKeyMaker(
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

    In fact, you can use any valid MongoDB query to get your molecule.

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
    ):
        """

        """

        database = mongo_client[database]
        self._molecules = database[molecule_collection]
        self._constructed_molecules = database[
            constructed_molecule_collection
        ]
        self._position_matrices = database[position_matrix_collection]
        self._jsonizer = jsonizer
        self._dejsonizer = dejsonizer

    def put(self, molecule):
        molecule = molecule.with_canonical_atom_ordering()
        json = self._jsonizer.to_json(molecule)
        self._position_matrices.insert_one(json['matrix'])
        self._molecules.insert_one(json['molecule'])
        self._constructed_molecules.insert_one(
            document=json['constructedMolecule'],
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

        return self._molecule_dejsonizer.from_json(
            json={
                'molecule': molecule_json,
                'constructedMolecule': constructed_molecule_json,
                'matrix': self._position_matrices.find_one(key),
            },
            building_blocks=building_blocks,
        )
