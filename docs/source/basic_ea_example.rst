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
with :mod:`stk`, called :mod:``. It used to make the library of
building blocks, which the EA searches, in order to find the optimal
molecules for our task. You can install it with::

    $ pip install


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

    # Go through 50 generations of the EA.
    for i, generation in enumerate(ea.get_generations(50)):

        # The generation object gives you access to the molecules
        # found in the generation.

You can place whatever code you like into the loop, for example,
you can write each molecule in each generation to a file

.. code-block:: python

    # Go through 50 generations of the EA.
    for i, generation in enumerate(ea.get_generations(50)):
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
:class:`.ConstructedMoleculeDatabase` guarantees the methods,
:meth:`~.ConstructedMoleculeDatabase.put`,
:meth:`~.ConstructedMoleculeDatabase.get` and
:meth:`~.ConstructedMoleculeDatabase.put_many`,
When using
:meth:`~.ConstructedMoleculeDatabase.put`, or
:meth:`~.ConstructedMoleculeDatabase.put_many`,
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

    # Go through 50 generations of the EA.
    for i, generation in enumerate(ea.get_generations(50)):
        molecules = (
            record.get_molecule()
            for record in generation.get_molecule_records()
        )
        db.put_many(molecules)


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
    for i, generation in enumerate(ea.get_generations(50)):
        molecules = (
            record.get_molecule()
            for record in generation.get_molecule_records()
        )
        db.put_many(molecules)
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

    # Go through 50 generations of the EA.
    generations = []
    for i, generation in enumerate(ea.get_generations(50)):
        molecules = (
            record.get_molecule()
            for record in generation.get_molecule_records()
        )
        db.put_many(molecules)
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
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        num_rotatable_bonds = rdkit.CalcNumRotatableBonds(
            mol=rdkit_molecule,
        )
        # Add 1 to the denominator to prevent division by 0.
        return 1 / (num_rotatable_bonds + 1)


Now that we have our function, we can turn it into a
:class:`.FitnessCalculator` by using :class:`.FitnessFunction`

.. code-block:: python

    fitness_calculator = stk.FitnessFunction(get_rigidity)

Now we only have to answer the second question,
*What kinds of molecular structures do I want to consider?*

.. code-block:: python

    import numpy as np

    def get_initial_population(fluoros, bromos, random_seed):
        generator = np.random.RandomState(random_seed)
        intial_fluoros = generator.random.choice(
            a=fluoros,
            size=5,
            replace=False,
        )
        initial bromos = generator.random.choice(
            a=bromos,
            size=5,
            replace=False,
        )

        for fluoro, bromo in it.product(
            intial_fluoros,
            initial_bromos,
        ):
            yield stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles=fluoro,
                            functional_groups=[stk.FluoroFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=bromo,
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=1,
                ),
            )

.. code-block:: python

    def get_functional_group_type(building_block):
        functional_group, = building_block.get_functional_groups(0)
        return functional_group.__class__

    # It's nice to get reproducible results.
    random_seed = 3
    ea = stk.EvolutionaryAlgorithm(
        initial_population=tuple(get_initial_population()),
        fitness_calculator=stk.FitnessFunction(get_rigidity),
        mutator=stk.RandomBuildingBlock(
            building_blocks=...,
            # All building blocks are replaceable.
            is_replaceable=lambda building_block: True,
            random_seed=random_seed,
        ),
        crosser=stk.GeneticRecombination(
            get_gene=get_functional_group_type,
        ),
        generation_selector=stk.Best(
            num_batches=25,
            duplicate_molecules=False,
        ),
        mutation_selector=stk.Roulette(
            num_batches=5
            random_seed=random_seed,
        ),
        crossover_selector=stk.Roulette(
            num_batches=3,
            batch_size=2,
            random_seed=random_seed,
        )
        # We don't need to do a normalization in this example.
        fitness_normalizer=stk.NullFitnessNormalizer(),

    )



Defining a Fitness Calculator
-----------------------------

One remaining EA component we need to define is a
:class:`.FitnessCalculator`.


Defining an Initial Population
------------------------------


Defining a Database
-------------------



Final Version
=============


This is a complete, basic EA. However, it has some obvious limitations:

* Only a single :class:`.Mutator` is defined. This means that only the
  space of bromo building blocks was explored.
* The :class:`.RandomBuildingBlock` mutator throws away an entire
  building block each time it performs a mutation. This means all the
  chemical information in it is lost. We would like a mutation which
  only *modifies* the building block, so that the new building block
  shares a some chemical features with the old one.
* The :class:`.FitnessCalculator` re-calculated the fitness value on
  molecules for which it had already calculated fitness values. This
  is OK in this example, but often a fitness calculation can be
  expensive and repeating it would seriously degrade the performance
  of our EA.

Next, you can read the intermediate tutorial, which will address all
of these limitations, and show you additional
customizations you can make to the EA, which allow it to be more
powerful and more efficient.
