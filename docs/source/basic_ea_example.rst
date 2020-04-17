================
Basic EA Example
================

Introduction
============

This tutorial will introduce you the basic components needed for
building an EA with :mod:`stk`. You can see a basic outline of the
EA in :class:`.EvolutionaryAlgorithm`.

You can get all the code associated with tutorial by running::

    $ git clone https://github.com/lukasturcani/basic_ea

This tutorial relies on a Python library, which does not come
with :mod:`stk`, called :mod:`vabene`. It used to make a
library of building blocks, which the EA uses, in order to find the
optimal molecules for our task. You can install it with::

    $ pip install vabene


EA Components
=============

There are a couple of components of the EA, which we need to define,
so that it can run. All of these provide an avenue for you to
customize the algorithm. We will need to use a

#. :class:`.Selector` - There are three different selectors we will
   need to define. One for selecting which molecules get mutated,
   one for selecting which molecules undergo crossover and one for
   selecting which molecules are part of the next generation of the
   EA.
#. :class:`.Mutator` - Used to carry out mutation operations on
   molecules.
#. :class:`.Crosser` - Used to carry out crossover operations on
   molecules.
#. :class:`.FitnessCalculator` - Used to calculate the fitness value
   of molecules.
#. :class:`.FitnessNormalizer` - Used to normalize the fitness values
   of molecules. Both the :class:`.FitnessCalculator` and the
   :class:`.FitnessNormalizer` assign fitness values to molecules,
   however the difference between them is that a
   :class:`.FitnessCalculator` can only return a fitness value based
   on the molecule itself, while a :class:`.FitnessNormalizer`
   can take into account all the other fitness values in the
   population. For example, it could divide the fitness values by the
   population average. A fitness value is always first assigned by a
   :class:`.FitnessCalculator`, and can then be optionally followed
   by a :class:`.FitnessNormalizer`


Basic EA Loop
=============

To build a really simple EA, we would begin by writing a really simple
Python script, such as

.. code-block:: python

    import stk


    ea = stk.EvolutionaryAlgorithm(
        # We will define these components later.
        initial_population=...,
        fitness_calculator=...,
        mutator=...,
        crosser=...,
        generation_selector=...,
        mutation_selector=...,
        crossover_selector=...,
        fitness_normalizer=...,
    )

    # Go through 15 generations of the EA.
    for i, generation in enumerate(ea.get_generations(15)):

        # The generation object gives you access to the molecules
        # found in the generation.

You can place whatever code you like into the loop, for example,
you can write each molecule in each generation to a file

.. code-block:: python

    # Go through 15 generations of the EA.
    for i, generation in enumerate(ea.get_generations(15)):
        # Go through the molecules in the generation, and write them
        # to a file.
        for molecule_id, molecule_record in enumerate(
            generation.get_molecule_records()
        ):
            molecule_record.get_molecule().write(
                path=f'generation_{i}_molecule_{molecule_id}.mol',
            )

While this is a perfectly valid EA loop, we can make it a lot better.


Adding a Database
-----------------

One of the main things that will significantly improve our quality of
life, is replacing our file writing, with a molecular database.
This means using a subclass of :class:`.ConstructedMoleculeDatabase`,
because the molecules produced by the EA are always constructed
molecules.

We won't define which :class:`.ConstructedMoleculeDatabase` we want to
use just yet, for now, all we need to know is that a
:class:`.ConstructedMoleculeDatabase` guarantees the methods
:meth:`~.ConstructedMoleculeDatabase.put` and
:meth:`~.ConstructedMoleculeDatabase.get`.
When using
:meth:`~.ConstructedMoleculeDatabase.put`,
the molecules are immediately deposited into the database, there is no
staging area.

So let's first assume we have defined some kind of
:class:`.ConstructedMoleculeDatabase`

.. code-block:: python

    # This will be a ConstructedMoleculeDatabase instance, which we
    # will define later.
    db = ...


Now we can modify the EA loop to use the database instead of
writing a bunch of files


.. code-block:: python

    # Go through 15 generations of the EA.
    for generation in ea.get_generations(15):
        for record in generation.get_molecule_records():
            db.put(record.get_molecule())


Already our EA loop is much nicer.


Plotting the EA Progress
========================

Usually, when we run an EA, we want to be able evaluate its
performance somehow. A very simple way to do this, is to plot how
the fitness of the population changes with generations. You
can use a :class:`.ProgressPlotter` to do this.

The :class:`.ProgressPlotter` needs to know what generations it
should plot, so we have to modify our loop so that it stores the
previous generations

.. code-block:: python

    generations = []
    for generation in ea.get_generations(15):
        for record in generation.get_molecule_records():
            db.put(record.get_molecule())
        generations.append(generation)

Now that we have the generations, we can use a
:class:`.ProgressPlotter` to plot them

.. code-block:: python

    fitness_progress = stk.ProgressPlotter(
        generations=generations,
        get_property=lambda record: record.get_fitness_value(),
        y_label='Fitness Value',
    )
    fitness_progress.write('fitness_progress.png')


Review
======

Ok, we now have a half-decent EA loop, so let's review it.

.. code-block:: python

    import stk

    db = ...
    ea = stk.EvolutionaryAlgorithm(
        initial_population=...,
        fitness_calculator=...,
        mutator=...,
        crosser=...,
        generation_selector=...,
        mutation_selector=...,
        crossover_selector=...,
        fitness_normalizer=...,
    )

    # Go through 15 generations of the EA.
    generations = []
    for generation in ea.get_generations(15):
        for record in generation.get_molecule_records():
            db.put(record.get_molecule())
        generations.append(generation)

    fitness_progress = stk.ProgressPlotter(
        generations=generations,
        get_property=lambda record: record.get_fitness_value(),
        y_label='Fitness Value',
    )
    fitness_progress.write('fitness_progress.png')

The only thing thats left to do, is define the components of the EA
that we want to use. There are a lot of options, so for the sake of
example, I will just use a couple of straight-forward ones.


Defining EA Components
======================

When defining EA components, there are two major questions that the
user must answer

* What molecular properties do I want to optimize?
* What kinds of molecular structures do I want to consider?

The user answers the first question by defining a
:class:`.FitnessCalculator`. The :class:`.FitnessCalculator` returns
a fitness value, and this is the value that the EA will optimize.
The simplest way to define a :class:`.FitnessCalculator` is to
first define a simple Python function, which takes a
:class:`.ConstructedMolecule` instance, and returns the fitness
value of that molecule.

For example, in many applications it is desirable to have rigid
molecules. One way to measure how rigid a molecule is, is to
calculate the number of rotatable bonds it has. The more rotatable
bonds, the less rigid the molecule. Therefore, if we want the EA to
produce rigid molecules, our fitness function should give a high
fitness to molecules with *few* rotatable bonds. We can therefore
define a function which returns the inverse of the number of rotatable
bonds in a molecule

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit

    def get_rigidity(molecule):
        num_rotatable_bonds = rdkit.CalcNumRotatableBonds(
            mol=rdkit_molecule,
        )
        # Add 1 to the denominator to prevent division by 0.
        return 1 / (num_rotatable_bonds + 1)

In addition to minimizing the number of rotatable bonds, we also
want to minimize the molecular complexity, so that they look
at least somewhat reasonable. :mod:`rdkit` provides a simple
way to calculate complexity

.. code-block:: python

    from rdkit.Chem.GraphDescriptors import BertzCT

    def get_complexity(molecule):
        return BertzCT(molecule)

Now that we can combine the two values into a single fitness value,
that the EA can optimize. There are multiple way to do this, but a
easy way is take the ratio of the value we want to maximize,
the rigidity, with the value with want to minimize, the complexity.

.. code-block:: python

    def get_fitness_value(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return (
            get_rigidity(rdkit_molecule)
            / get_complexity(rdkit_molecule)
        )

Now that we have our function, we can turn it into a
:class:`.FitnessCalculator` by using :class:`.FitnessFunction`

.. code-block:: python

    fitness_calculator = stk.FitnessFunction(get_fitness_value)

Now we only have to answer the second question,
*What kinds of molecular structures do I want to consider?*

The user answers this question by defining an initial population of
molecules the EA should use, as well as the mutation and crossover
operations. These operations will determine which molecules the EA
can construct.

Lets begin by defining an initial population. The first thing we will
need is a set of building blocks, with which we can build our
molecules. We can use the :mod:`vabene` library to create random
molecular graphs, which can be used to make our building blocks

.. code-block:: python

    import vabene as vb
    import numpy as np
    import itertools as it

    def get_building_block(
        generator,
        atomic_number,
        functional_group_factory,
    ):
        # The number of atoms, excluding hydrogen, in our building
        # block.
        num_atoms = generator.randint(7, 15)
        # The distance between the bromine or fluorine atoms in our
        # building block.
        fg_separation = generator.randint(1, num_atoms-3)

        atom_factory = vb.RandomAtomFactory(
            # Our building blocks will be composed from a mix of
            # carbon atoms with maximum valance of 4, carbon atoms
            # with a maximum valence of 3 and oxygen atoms with a
            # maximum valence of 2.
            atoms=(
                vb.Atom(6, 0, 4),
                vb.Atom(6, 0, 3),
                vb.Atom(8, 0, 2),
            ),
            # All of our building blocks will have 2 halogen atoms,
            # separated by a random number of carbon atoms.
            required_atoms=(
                (vb.Atom(atomic_number, 0, 1), )
                +
                (vb.Atom(6, 0, 4), ) * fg_separation
                +
                (vb.Atom(atomic_number, 0, 1), )
            ),
            num_atoms=num_atoms,
            random_seed=generator.randint(0, 1000),
        )
        atoms = tuple(atom_factory.get_atoms())
        bond_factory = vb.RandomBondFactory(
            required_bonds=tuple(
                vb.Bond(i, i+1, 1) for i in range(fg_separation+1)
            ),
            random_seed=generator.randint(0, 1000),
        )
        bonds = bond_factory.get_bonds(atoms)

        building_block = stk.BuildingBlock.init_from_rdkit_mol(
            molecule=vabene_to_rdkit(vb.Molecule(atoms, bonds)),
            functional_groups=[functional_group_factory],
        )
        # We can give random coordinates to the building block,
        # because it's fast and doesn't matter in this case.
        return building_block.with_position_matrix(
            position_matrix=generator.uniform(
                low=-100,
                high=100,
                size=(building_block.get_num_atoms(), 3),
            ),
        )

Note that we need to define a function which converts
:mod:`vabene` molecules into :mod:`rdkit` molecules

.. code-block:: python

    def vabene_to_rdkit(molecule):
        editable = rdkit.EditableMol(rdkit.Mol())
        for atom in molecule.get_atoms():
            rdkit_atom = rdkit.Atom(atom.get_atomic_number())
            rdkit_atom.SetFormalCharge(atom.get_charge())
            editable.AddAtom(rdkit_atom)

        for bond in molecule.get_bonds():
            editable.AddBond(
                beginAtomIdx=bond.get_atom1_id(),
                endAtomIdx=bond.get_atom2_id(),
                order=rdkit.BondType(bond.get_order()),
            )

        rdkit_molecule = editable.GetMol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit_molecule = rdkit.AddHs(rdkit_molecule)
        rdkit.Kekulize(rdkit_molecule)
        rdkit_molecule.AddConformer(
            conf=rdkit.Conformer(rdkit_molecule.GetNumAtoms()),
        )
        return rdkit_molecule

Once these functions are defined, we can use :func:`get_building_block`
to generate out building blocks

.. code-block:: python

    # Use a random seed to get reproducible results.
    random_seed = 4
    generator = np.random.RandomState(random_seed)

    # Make 1000 fluoro building bocks.
    fluoros = tuple(
        get_building_block(generator, 9, stk.FluoroFactory())
        for i in range(1000)
    )
    # Make 1000 bromo building blocks.
    bromos = tuple(
        get_building_block(generator, 35, stk.BromoFactory())
        for i in range(1000)
    )


In this example, the EA will create ``AB`` dimers, using the
:class:`.Linear`  topology graph. The initial population of 25 such
dimers can be made by taking the first 5 ``bromo`` and ``fluoro``
building blocks

.. code-block:: python

    def get_initial_population(fluoros, bromos):
        for fluoro, bromo in it.product(fluoros[:5], bromos[:5]):
            yield stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(fluoro, bromo),
                    repeating_unit='AB',
                    num_repeating_units=1,
                ),
            )

    initial_population = tuple(get_initial_population(fluoros, bromos)

Next, we can define our mutation operations. There are a multiple
options, as you can see in the sidebar. One thing that you might
notice immediately, is that there are multiple :class:`.Mutator`
types you would like to use during the EA, but the
:class:`.EvolutionaryAlgorithm` only takes a single
:class:`.Mutator`. To get around this, we can use a compound
:class:`.Mutator`, such as the :class:`.RandomMutator`. When you create
a :class:`.RandomMutator`, you define it in terms of other mutators
you want to use, for example

.. code-block:: python

    def get_functional_group_type(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__

    def is_fluoro(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__ is stk.Fluoro

    def is_bromo(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__ is stk.Bromo

    mutator = stk.RandomMutator(
        mutators=(
            stk.RandomBuildingBlock(
                building_blocks=fluoros,
                is_replaceable=is_fluoro,
                random_seed=generator.randint(0, 1000),
            ),
            stk.SimilarBuildingBlock(
                building_blocks=fluoros,
                is_replaceable=is_fluoro,
                random_seed=generator.randint(0, 1000),
            ),
            stk.RandomBuildingBlock(
                building_blocks=bromos,
                is_replaceable=is_bromo,
                random_seed=generator.randint(0, 1000),
            ),
            stk.SimilarBuildingBlock(
                building_blocks=bromos,
                is_replaceable=is_bromo,
                random_seed=generator.randint(0, 1000),
            ),
        ),
        random_seed=generator.randint(0, 1000),
    )

When :meth:`~.Mutator.mutate` is called on a :class:`.RandomMutator`,
it randomly selects one of the mutators you gave it during
initialization, and asks it to perform the mutation operation on its
behalf. In this way, all of the mutators you provided it will get used
during the EA.

Now we can put all of these components together, and fill in the
remaining ones too

.. code-block:: python

    ea = stk.EvolutionaryAlgorithm(
        initial_population=tuple(
            get_initial_population(fluoros, bromos)
        ),
        fitness_calculator=stk.FitnessFunction(get_fitness_value),
        mutator=stk.RandomMutator(
            mutators=(
                stk.RandomBuildingBlock(
                    building_blocks=fluoros,
                    is_replaceable=is_fluoro,
                    random_seed=generator.randint(0, 1000),
                ),
                stk.SimilarBuildingBlock(
                    building_blocks=fluoros,
                    is_replaceable=is_fluoro,
                    random_seed=generator.randint(0, 1000),
                ),
                stk.RandomBuildingBlock(
                    building_blocks=bromos,
                    is_replaceable=is_bromo,
                    random_seed=generator.randint(0, 1000),
                ),
                stk.SimilarBuildingBlock(
                    building_blocks=bromos,
                    is_replaceable=is_bromo,
                    random_seed=generator.randint(0, 1000),
                ),
            ),
            random_seed=generator.randint(0, 1000),
        ),
        crosser=stk.GeneticRecombination(
            get_gene=get_functional_group_type,
        ),
        generation_selector=stk.Best(
            num_batches=25,
            duplicate_molecules=False,
        ),
        mutation_selector=stk.Roulette(
            num_batches=5,
            random_seed=generator.randint(0, 1000),
        ),
        crossover_selector=stk.Roulette(
            num_batches=3,
            batch_size=2,
            random_seed=generator.randint(0, 1000),
        ),
        # We don't need to do a normalization in this example.
        fitness_normalizer=stk.NullFitnessNormalizer(),
    )

Defining a Database
-------------------

The last thing we need to do is define the database. The default
database of :mod:`stk` is MongoDB, which can be used with
:class:`.ConstructedMoleculeMongoDb`. Before using this class, make
sure you have :mod:`pymongo` and that its working properly. I recommend
reading at least the introductory and installation
documentation of :mod:`pymongo` before using this class. Those
docs can be found here__.

__ https://api.mongodb.com/python/current/

Note that this is easy to do, and well worth the minimal effort it
requires to setup. Obviously, if you really don't want to
use the database, you do not have to create it, and you can remove
references to it in your EA loop.

Assuming everything is setup, we can create our database instance

.. code-block:: python

    # pymongo does not come with stk, you have to install it
    # explicitly with "pip install pymongo".
    import pymongo

    # Connect to a MongoDB. This example connects to a local
    # MongoDB, but you can connect to a remote DB too with
    # MongoClient() - read the documentation for pymongo to see how
    # to do that.
    client = pymongo.MongoClient()
    db = stk.ConstructedMoleculeMongoDb(client)


Final Version
=============

The final version of our code is

.. code-block:: python

    import stk
    import rdkit.Chem.AllChem as rdkit
    from rdkit.Chem.GraphDescriptors import BertzCT
    import pymongo
    import vabene as vb
    import numpy as np
    import itertools as it
    import logging


    logger = logging.getLogger(__name__)


    def vabene_to_rdkit(molecule):
        editable = rdkit.EditableMol(rdkit.Mol())
        for atom in molecule.get_atoms():
            rdkit_atom = rdkit.Atom(atom.get_atomic_number())
            rdkit_atom.SetFormalCharge(atom.get_charge())
            editable.AddAtom(rdkit_atom)

        for bond in molecule.get_bonds():
            editable.AddBond(
                beginAtomIdx=bond.get_atom1_id(),
                endAtomIdx=bond.get_atom2_id(),
                order=rdkit.BondType(bond.get_order()),
            )

        rdkit_molecule = editable.GetMol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit_molecule = rdkit.AddHs(rdkit_molecule)
        rdkit.Kekulize(rdkit_molecule)
        rdkit_molecule.AddConformer(
            conf=rdkit.Conformer(rdkit_molecule.GetNumAtoms()),
        )
        return rdkit_molecule


    def get_building_block(
        generator,
        atomic_number,
        functional_group_factory,
    ):
        # The number of atoms, excluding hydrogen, in our building
        # block.
        num_atoms = generator.randint(7, 15)
        # The distance between the bromine or fluorine atoms in our
        # building block.
        fg_separation = generator.randint(1, num_atoms-3)

        atom_factory = vb.RandomAtomFactory(
            # Our building blocks will be composed from a mix of
            # carbon atoms with maximum valance of 4, carbon atoms
            # with a maximum valence of 3 and oxygen atoms with a
            # maximum valence of 2.
            atoms=(
                vb.Atom(6, 0, 4),
                vb.Atom(6, 0, 3),
                vb.Atom(8, 0, 2),
            ),
            # All of our building blocks will have 2 halogen atoms,
            # separated by a random number of carbon atoms.
            required_atoms=(
                (vb.Atom(atomic_number, 0, 1), )
                +
                (vb.Atom(6, 0, 4), ) * fg_separation
                +
                (vb.Atom(atomic_number, 0, 1), )
            ),
            num_atoms=num_atoms,
            random_seed=generator.randint(0, 1000),
        )
        atoms = tuple(atom_factory.get_atoms())
        bond_factory = vb.RandomBondFactory(
            required_bonds=tuple(
                vb.Bond(i, i+1, 1) for i in range(fg_separation+1)
            ),
            random_seed=generator.randint(0, 1000),
        )
        bonds = bond_factory.get_bonds(atoms)

        building_block = stk.BuildingBlock.init_from_rdkit_mol(
            molecule=vabene_to_rdkit(vb.Molecule(atoms, bonds)),
            functional_groups=[functional_group_factory],
        )
        # We can give random coordinates to the building block,
        # because it's fast and doesn't matter in this case.
        return building_block.with_position_matrix(
            position_matrix=generator.uniform(
                low=-100,
                high=100,
                size=(building_block.get_num_atoms(), 3),
            ),
        )


    def get_initial_population(fluoros, bromos):
        for fluoro, bromo in it.product(fluoros[:5], bromos[:5]):
            yield stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(fluoro, bromo),
                    repeating_unit='AB',
                    num_repeating_units=1,
                ),
            )


    def get_rigidity(molecule):
        num_rotatable_bonds = rdkit.CalcNumRotatableBonds(
            mol=rdkit_molecule,
        )
        # Add 1 to the denominator to prevent division by 0.
        return 1 / (num_rotatable_bonds + 1)

    def get_complexity(molecule):
        return BertzCT(molecule)


    def get_fitness_value(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return (
            get_rigidity(rdkit_molecule)
            / get_complexity(rdkit_molecule)
        )


    def get_functional_group_type(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__


    def is_fluoro(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__ is stk.Fluoro


    def is_bromo(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__ is stk.Bromo


    def get_num_rotatable_bonds(record):
        molecule = record.get_molecule().to_rdkit_mol()
        rdkit.SanitizeMol(molecule)
        return rdkit.CalcNumRotatableBonds(molecule)


    def main():
        logging.basicConfig(level=logging.INFO)

        # Use a random seed to get reproducible results.
        random_seed = 4
        generator = np.random.RandomState(random_seed)

        logger.info('Making building blocks.')

        # Make 1000 fluoro building bocks.
        fluoros = tuple(
            get_building_block(generator, 9, stk.FluoroFactory())
            for i in range(1000)
        )
        # Make 1000 bromo building blocks.
        bromos = tuple(
            get_building_block(generator, 35, stk.BromoFactory())
            for i in range(1000)
        )

        db = stk.ConstructedMoleculeMongoDb(pymongo.MongoClient())
        ea = stk.EvolutionaryAlgorithm(
            initial_population=tuple(
                get_initial_population(fluoros, bromos)
            ),
            fitness_calculator=stk.FitnessFunction(get_fitness_value),
            mutator=stk.RandomMutator(
                mutators=(
                    stk.RandomBuildingBlock(
                        building_blocks=fluoros,
                        is_replaceable=is_fluoro,
                        random_seed=generator.randint(0, 1000),
                    ),
                    stk.SimilarBuildingBlock(
                        building_blocks=fluoros,
                        is_replaceable=is_fluoro,
                        random_seed=generator.randint(0, 1000),
                    ),
                    stk.RandomBuildingBlock(
                        building_blocks=bromos,
                        is_replaceable=is_bromo,
                        random_seed=generator.randint(0, 1000),
                    ),
                    stk.SimilarBuildingBlock(
                        building_blocks=bromos,
                        is_replaceable=is_bromo,
                        random_seed=generator.randint(0, 1000),
                    ),
                ),
                random_seed=generator.randint(0, 1000),
            ),
            crosser=stk.GeneticRecombination(
                get_gene=get_functional_group_type,
            ),
            generation_selector=stk.Best(
                num_batches=25,
                duplicate_molecules=False,
            ),
            mutation_selector=stk.Roulette(
                num_batches=5,
                random_seed=generator.randint(0, 1000),
            ),
            crossover_selector=stk.Roulette(
                num_batches=3,
                batch_size=2,
                random_seed=generator.randint(0, 1000),
            ),
            # We don't need to do a normalization in this example.
            fitness_normalizer=stk.NullFitnessNormalizer(),
        )

        logger.info('Starting EA.')

        generations = []
        for generation in ea.get_generations(15):
            for record in generation.get_molecule_records():
                db.put(record.get_molecule())
            generations.append(generation)

        logger.info('Making fitness plot.')

        fitness_progress = stk.ProgressPlotter(
            generations=generations,
            get_property=lambda record: record.get_fitness_value(),
            y_label='Fitness Value',
        )
        fitness_progress.write('fitness_progress.png')

        logger.info('Making rotatable bonds plot.')

        rotatable_bonds_progress = stk.ProgressPlotter(
            generations=generations,
            get_property=get_num_rotatable_bonds,
            y_label='Number of Rotatable Bonds',
        )
        rotatable_bonds_progress.write('rotatable_bonds_progress.png')

    if __name__ == '__main__':
        main()


The plot of fitness we produced looks like this:

.. image:: https://i.imgur.com/eNc0f4B.png

which addmittedly doesn't look that great. However, since the
fitness value is just a ratio of two numbers, it doesn't really tell
us much. A better thing to look at is the plot for the number of
rotatable bonds

.. image:: https://i.imgur.com/7retPh8.png


Clearly, our EA was able to quickly minimize the number of rotatable
bonds to the lowest possible value across all members of the
population.

We can also compare the molecules in the initial population

.. image:: https://i.imgur.com/qEvCMk0.png

to those in the final population

.. image:: https://i.imgur.com/alBW1OS.png

where the hydrogen atoms have been left out for clarity. Again, the
EA obviously managed to reduce the complexity of the molecules as well.

Next, you can read the intermediate tutorial, which will show you
some additional customization options for the EA.
