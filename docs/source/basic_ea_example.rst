Basic EA Example
================

Introduction
============

This tutorial will introduce you the basic components needed for
building an EA with :mod:`stk`. You can see a basic outline of the
EA in :class:`.EvolutionaryAlgorithm`.

You can get all the code associated with tutorial by running::

    $ git clone https://github.com/lukasturcani/basic_ea


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

For most most of the EA components, it doesn't really matter what
molecules you are trying to design, any option will do


.. code-block:: python

    ea = stk.EvolutionaryAlgorithm(
        initial_population=...,
        fitness_calculator=...,
        mutator=...,
        crosser=stk.GeneticRecombination(
        ),
        generation_selector=stk.AboveAverage(
            num_batches=20,
            duplicate_molecules=False,
        ),
        mutation_selector=stk.Roulette(
            num_batches=5
            # It's nice to get reproducible results.
            random_seed=5,
        ),
        crossover_selector=stk.Roulette(
            num_batches=3,
            batch_size=2,
            random_seed=8,
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
