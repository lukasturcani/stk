"""
Property Vector
===============

"""

from .fitness_calculator import FitnessCalculator


class PropertyVector(FitnessCalculator):
    """
    Uses multiple molecular properties as a fitness value.

    Examples
    --------
    *Calculating Fitness Values*

    .. code-block:: python

        import stk

        # First, create the functions which calculate the properties
        # of molecules.
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        def get_num_bonds(molecule):
            return molecule.get_num_bonds()

        def get_diameter(molecule):
            return molecule.get_maximum_diameter()

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
        value = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

    *Storing Fitness Values in a Database*

    Sometimes you want to store fitness values in a database, you
    can do this by providing the `output_database` parameter.

    .. code-block:: python

        import stk
        # pymongo does not come with stk, you have to install it
        # separately with "pip install pymongo"
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
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        def get_num_bonds(molecule):
            return molecule.get_num_bonds()

        def get_diameter(molecule):
            return molecule.get_maximum_diameter()

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
        value1 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

        # You can retrieve the fitness values from the database.
        value = fitness_db.get(stk.BuildingBlock('BrCCBr'))

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

    .. code-block:: python

        import stk

        # pymongo does not come with stk, you have to install it
        # separately with "pip install pymongo"
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
        def get_num_atoms(molecule):
            return molecule.get_num_atoms()

        def get_num_bonds(molecule):
            return molecule.get_num_bonds()

        def get_diameter(molecule):
            return molecule.get_maximum_diameter()

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
        value1 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )
        # This will not re-calculate the fitness value, instead,
        # value1 will be retrieved from the database.
        value2 = fitness_calculator.get_fitness_value(
            molecule=stk.BuildingBlock('BrCCBr'),
        )

    """

    def __init__(
        self,
        property_functions,
        input_database=None,
        output_database=None,
    ):

        """
        Initialize a :class:`.PropertyVector` instance.

        Parameters
        ----------
        property_fns : :class:`tuple` of :class:`callable`
            A group of :class:`function`, each of which is used to
            calculate a single property of the molecule. Each function
            must take one parameter, `mol`, which accepts
            a :class:`.Molecule` object. This is the molecule used to
            calculate the property.

        input_database : :class:`.ValueDatabase`, optional
            A database to check before calling `fitness_function`. If a
            fitness value exists for a molecule in the database, the
            stored value is returned, instead of calling
            `fitness_function`.

        output_database : :class:`.ValueDatabase`, optional
            A database into which the calculate fitness value is
            placed.

        """

        self._property_functions = property_functions
        self._input_database = input_database
        self._output_database = output_database

    def get_fitness_value(self, molecule):
        return tuple(
            property_function(molecule)
            for property_function in self._property_functions
        )
