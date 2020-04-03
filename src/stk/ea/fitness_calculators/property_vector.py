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

        # Create the molecules which are to have fitness values
        # evaluated.
        mol1 = stk.BuildingBlock('NCCN', ['amine'])
        mol2 = stk.BuildingBlock('NC[Si]CCCN', ['amine'])
        mol3 = stk.BuildingBlock('O=CCC(C=O)CCC=O', ['aldehyde'])

        # Create the functions which calculate the molecule properties.
        def atom_number(mol):
            return len(mol.atoms)

        def diameter(mol):
            return mol.maximum_diamater()

        def energy(mol):
            energy_calculator = stk.MMFFEnergy()
            return energy_calculator.get_energy(mol)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(
            atom_number,
            diameter,
            energy
        )

        # Calculate the fitness vector of mol1. It will be a list
        # holding the number of atoms, diameter and energy of mol1,
        # respectively.
        mol1_fitness = fitness_calculator.get_fitness(mol1)

        # Calculate the fitness vector of mol2. It will be a list
        # holding the number of atoms, diameter and energy of mol2,
        # respectively.
        mol2_fitness = fitness_calculator.get_fitness(mol2)

        # Calculate the fitness vector of mol3. It will be a list
        # holding the number of atoms, diameter and energy of mol3,
        # respectively.
        mol3_fitness = fitness_calculator.get_fitness(mol3)


    Use on :class:`.ConstructedMolecule` objects

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # calculate.
        bb1 = stk.BuildingBlock('[Br]CC[Br]', ['bromine'])
        polymer1 = stk.ConstructedMolecule(
            building_blocks=[bb1],
            topology_graph=stk.polymer.Linear('A', [0], n=5)
        )

        bb2 = stk.BuildingBlock('[Br]CCNNCC[Br]', ['bromine'])
        polymer2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        # Create the functions which calculate the molecule properties.
        def atom_number(mol):
            return len(mol.atoms)

        def diameter(mol):
            return mol.maximum_diameter()

        def monomer_number(mol):
            return sum(mol.building_block_counter.values())

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(
            atom_number,
            diameter,
            monomer_number
        )

        # Calculate the fitness vector of polymer1. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer1, respectively.
        polymer1_fitness = fitness_calculator.get_fitness(polymer1)

        # Calculate the fitness vector of polymer2. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer2, respectively.
        polymer2_fitness = fitness_calculator.get_fitness(polymer2)

    """

    def __init__(
        self,
        property_functions,
        input_database=None,
        output_database=None,
    ):

        """
        Initialize a :class:`CageFitness` instance.

        Parameters
        ----------
        *property_fns : :class:`tuple` of :class:`function`
            A group of :class:`function`, each of which is used to
            calculate a single property of the molecule. Each function
            must take one parameter, `mol`, which accepts
            a :class:`.Molecule` object. This is the molecule used to
            calculate the property.

        use_cache : :class:`bool`, optional
            If ``True`` then fitness values for molecules already held
            in the cache are not re-calculated but the value already
            stored is used.

        """

        self._property_fns = property_fns
        super().__init__(use_cache=use_cache)

    def _get_fitness(self, mol):
        """
        Get the fitness of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose property vector should be calculated.

        Returns
        -------
        :class:`list`
            A :class:`list` of properties of the `mol`. Will be
            ``None`` if any of the properties fail.

        """

        property_vector = []
        for property_fn in self._property_fns:
            logger.info(
                f'Using {property_fn.__name__} on "{mol}".'
            )
            property_vector.append(property_fn(mol))
        return property_vector
