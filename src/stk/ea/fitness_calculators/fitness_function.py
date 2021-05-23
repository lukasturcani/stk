"""
Fitness Calculator
==================

"""

from .fitness_calculator import FitnessCalculator


class FitnessFunction(FitnessCalculator):
    """
    Takes a function and uses it as a fitness calculator.

    Examples
    --------
    *Calculating Fitness Values*

    .. testcode:: calculating-fitness-values

        import stk

        # First create a function to use as the fitness function. For
        # example, assume that the fitness value is given by the number
        # of atoms.
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        # Then use it to create a fitness calculator.
        fitness_calculator = stk.FitnessFunction(
            fitness_function=get_num_atoms,
        )

        # Use the fitness calculator to get fitness values.
        value1 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

    .. testcode:: calculating-fitness-values
        :hide:

        assert value1 == stk.BuildingBlock('BrCCBr').get_num_atoms()

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

        # Define a fitness function.
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        # Create a fitness calculator.
        fitness_calculator = stk.FitnessFunction(
            fitness_function=get_num_atoms,
            output_database=fitness_db,
        )

        # Calculate fitness values.
        value = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

        # You can retrieve the fitness values from the database.
        value1 = fitness_db.get(stk.BuildingBlock('BrCCBr'))

    .. testcode:: storing-fitness-values-in-a-database
        :hide:

        assert value == value1

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

        # Define a fitness function.
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        # Create a fitness calculator.
        fitness_calculator = stk.FitnessFunction(
            fitness_function=get_num_atoms,
            input_database=fitness_db,
            output_database=fitness_db,
        )

        # Assuming that a fitness value for this molecule was not
        # deposited into the database in a previous session, this
        # will calculate the fitness value.
        value1 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )
        # This will not re-calculate the fitness value, instead,
        # value1 will be retrieved from the database.
        value2 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

    .. testcode:: caching-fitness-values
        :hide:

        assert value1 == value2

    .. testcleanup:: caching-fitness-values

        stk.ValueMongoDb = _old_init
        pymongo.MongoClient().drop_database(_test_database)
        pymongo.MongoClient = _mongo_client

    """

    def __init__(
        self,
        fitness_function,
        input_database=None,
        output_database=None,
    ):
        """
        Initialize a :class:`.FitnessFunction` instance.

        Parameters
        ----------
        fitness_function : :class:`callable`
            Takes a single parameter, the
            :class:`.ConstructedMolecule` whose
            fitness needs to be calculated, and returns its
            fitness value.

        input_database : :class:`.ValueDatabase`, optional
            A database to check before calling `fitness_function`. If a
            fitness value exists for a molecule in the database, the
            stored value is returned, instead of calling
            `fitness_function`.

        output_database : :class:`.ValueDatabase`, optional
            A database into which the calculate fitness value is
            placed.

        """

        self._fitness_function = fitness_function
        self._input_database = input_database
        self._output_database = output_database

    def get_fitness_value(self, molecule):
        if self._input_database is not None:
            try:
                fitness_value = self._input_database.get(molecule)
            except KeyError:
                fitness_value = self._fitness_function(molecule)
        else:
            fitness_value = self._fitness_function(molecule)

        if self._output_database is not None:
            self._output_database.put(molecule, fitness_value)

        return fitness_value
