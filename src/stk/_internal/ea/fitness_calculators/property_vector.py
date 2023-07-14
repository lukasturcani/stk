import typing
from collections.abc import Callable, Iterable
from typing import Any

from stk._internal.databases.value import ValueDatabase
from stk._internal.ea.molecule_record import MoleculeRecord

from .fitness_calculator import FitnessCalculator

T = typing.TypeVar("T", bound=MoleculeRecord)


class PropertyVector(FitnessCalculator[T]):
    """
    Uses multiple molecular properties as a fitness value.

    Examples:

        *Calculating Fitness Values*

        .. testcode:: calculating-fitness-values

            import stk

            # First, create the functions which calculate the properties
            # of molecules.
            def get_num_atoms(record):
                return record.get_molecule().get_num_atoms()

            def get_num_bonds(record):
                return record.get_molecule().get_num_bonds()

            def get_diameter(record):
                return record.get_molecule().get_maximum_diameter()

            # Next, create the fitness calculator.
            fitness_calculator = stk.PropertyVector(
                property_functions=(
                    get_num_atoms,
                    get_num_bonds,
                    get_diameter,
                ),
            )

            # Calculate the fitness value of a molecule.
            # "value" is a tuple, holding the number of atoms, number of
            # bonds and the diameter of the molecule.
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(stk.BuildingBlock('BrCCBr'), ),
                    repeating_unit='A',
                    num_repeating_units=1,
                ),
            )
            value = fitness_calculator.get_fitness_value(record)

        .. testcode:: calculating-fitness-values
            :hide:

            _bb = stk.BuildingBlock('BrCCBr')
            assert value == (
                _bb.get_num_atoms(),
                _bb.get_num_bonds(),
                _bb.get_maximum_diameter(),
            )

        *Storing Fitness Values in a Database*

        Sometimes you want to store fitness values in a database, you
        can do this by providing the `output_database` parameter.

        .. testsetup:: storing-fitness-values-in-a-database

            import stk

            # Change the database used, so that when a developer
            # runs the doctests locally, their "stk" database is not
            # contaminated.
            _test_database = '_stk_doctest_database'
            _old_init = stk.ValueMongoDb
            stk.ValueMongoDb = lambda mongo_client, collection: (
                _old_init(
                    mongo_client=mongo_client,
                    database=_test_database,
                    collection=collection,
                )
            )

            # Change the database MongoClient will connect to.

            import os
            import pymongo

            _mongo_client = pymongo.MongoClient
            _mongodb_uri = os.environ.get(
                'MONGODB_URI',
                'mongodb://localhost:27017/'
            )
            pymongo.MongoClient = lambda: _mongo_client(_mongodb_uri)

        .. testcode:: storing-fitness-values-in-a-database

            import stk
            import pymongo

            # Create a database which stores the fitness value of each
            # molecule.
            fitness_db = stk.ValueMongoDb(
                # This connects to a local database - so make sure you have
                # local MongoDB server running. You can also connect to
                # a remote MongoDB with MongoClient(), read to pymongo
                # docs to see how to do that.
                mongo_client=pymongo.MongoClient(),
                collection='fitness_values',
            )

            # Define the functions which calculate molecular properties.
            def get_num_atoms(record):
                return record.get_molecule().get_num_atoms()

            def get_num_bonds(record):
                return record.get_molecule().get_num_bonds()

            def get_diameter(record):
                return record.get_molecule().get_maximum_diameter()

            # Create the fitness calculator.
            fitness_calculator = stk.PropertyVector(
                property_functions=(
                    get_num_atoms,
                    get_num_bonds,
                    get_diameter,
                ),
                output_database=fitness_db,
            )

            # Calculate fitness values.
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(stk.BuildingBlock('BrCCBr'), ),
                    repeating_unit='A',
                    num_repeating_units=1,
                ),
            )
            value1 = fitness_calculator.get_fitness_value(record)

            # You can retrieve the fitness values from the database.
            value2 = fitness_db.get(record.get_molecule())

        .. testcode:: storing-fitness-values-in-a-database
            :hide:

            assert value1 == tuple(value2)

        .. testcleanup:: storing-fitness-values-in-a-database

            stk.ValueMongoDb = _old_init
            pymongo.MongoClient().drop_database(_test_database)
            pymongo.MongoClient = _mongo_client

        *Caching Fitness Values*

        Usually, if you calculate the fitness value of a molecule, you
        do not want to re-calculate it, because this may be expensive,
        and the fitness value is going to be the
        same anyway. By using the `input_database` parameter, together
        with the `output_database` parameter, you can make sure you store
        and retrieve calculated fitness values instead of repeating the
        same calculation multiple times.

        The `input_database` is checked before a calculation happens, to
        see if the value already exists, while the `output_database` has
        the calculated fitness value deposited into it.

        .. testsetup:: caching-fitness-values

            import stk

            # Change the database used, so that when a developer
            # runs the doctests locally, their "stk" database is not
            # contaminated.
            _test_database = '_stk_doctest_database'
            _old_init = stk.ValueMongoDb
            stk.ValueMongoDb = lambda mongo_client, collection: (
                _old_init(
                    mongo_client=mongo_client,
                    database=_test_database,
                    collection=collection,
                )
            )

            # Change the database MongoClient will connect to.

            import os
            import pymongo

            _mongo_client = pymongo.MongoClient
            _mongodb_uri = os.environ.get(
                'MONGODB_URI',
                'mongodb://localhost:27017/'
            )
            pymongo.MongoClient = lambda: _mongo_client(_mongodb_uri)

        .. testcode:: caching-fitness-values

            import stk
            import pymongo

            # You can use the same database for both the input_database
            # and output_database parameters.
            fitness_db = stk.ValueMongoDb(
                # This connects to a local database - so make sure you have
                # local MongoDB server running. You can also connect to
                # a remote MongoDB with MongoClient(), read to pymongo
                # docs to see how to do that.
                mongo_client=pymongo.MongoClient(),
                collection='fitness_values',
            )

            # Define the functions which calculate molecular properties.
            def get_num_atoms(record):
                return record.get_molecule().get_num_atoms()

            def get_num_bonds(record):
                return record.get_molecule().get_num_bonds()

            def get_diameter(record):
                return record.get_molecule().get_maximum_diameter()

            # Create the fitness calculator.
            fitness_calculator = stk.PropertyVector(
                property_functions=(
                    get_num_atoms,
                    get_num_bonds,
                    get_diameter,
                ),
                input_database=fitness_db,
                output_database=fitness_db,
            )

            # Assuming that a fitness value for this molecule was not
            # deposited into the database in a previous session, this
            # will calculate the fitness value.
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(stk.BuildingBlock('BrCCBr'), ),
                    repeating_unit='A',
                    num_repeating_units=1,
                ),
            )
            value1 = fitness_calculator.get_fitness_value(record)
            # This will not re-calculate the fitness value, instead,
            # value1 will be retrieved from the database.
            value2 = fitness_calculator.get_fitness_value(record)

        .. testcode:: caching-fitness-values
            :hide:

            value3 = fitness_calculator.get_fitness_value(record)
            assert value2 is value3

        .. testcleanup:: caching-fitness-values

            stk.ValueMongoDb = _old_init
            pymongo.MongoClient().drop_database(_test_database)
            pymongo.MongoClient = _mongo_client

    """

    def __init__(
        self,
        property_functions: Iterable[Callable[[T], Any]],
        input_database: ValueDatabase | None = None,
        output_database: ValueDatabase | None = None,
    ) -> None:
        """
        Parameters:
            property_functions \
(list[~collections.abc.Callable[[T], typing.Any]]):
                A group of functions, each of which is used to
                calculate a single property of the molecule.

            input_database:
                A database to check before calling `fitness_function`. If a
                fitness value exists for a molecule in the database, the
                stored value is returned, instead of calling
                `fitness_function`.

            output_database:
                A database into which the calculate fitness value is
                placed.
        """
        self._property_functions = tuple(property_functions)
        self._input_database = input_database
        self._output_database = output_database

    def get_fitness_value(self, record: T) -> typing.Any:
        if self._input_database is not None:
            try:
                fitness_value = self._input_database.get(
                    molecule=record.get_molecule(),
                )
            except KeyError:
                fitness_value = tuple(
                    property_function(record)
                    for property_function in self._property_functions
                )
        else:
            fitness_value = tuple(
                property_function(record)
                for property_function in self._property_functions
            )

        if self._output_database is not None:
            self._output_database.put(record.get_molecule(), fitness_value)

        return fitness_value
