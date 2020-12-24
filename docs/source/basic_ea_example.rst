================
Basic EA Example
================

Introduction
============

This tutorial will introduce you to the basic components needed for
building an EA with :mod:`stk`. You can see a basic outline of the
EA in :class:`.EvolutionaryAlgorithm`.

You can get all the code associated with this tutorial by running::

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

An important class in our EA is also the :class:`.MoleculeRecord`.
The :class:`.MoleculeRecord` holds a molecule being considered by the
EA, as well as additional relevant information, such as its fitness
value.

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

See Also
    :ref:`placing-and-retrieving-molecules-from-a-database`



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

    # Go through 50 generations of the EA.
    for generation in ea.get_generations(50):
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
    for generation in ea.get_generations(50):
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
    )

    # Go through 50 generations of the EA.
    generations = []
    for generation in ea.get_generations(50):
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
user must answer:

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
want to minimize the molecular complexity, so that molecules made
by the EA look at least somewhat reasonable. :mod:`rdkit` provides a
function called :func:`BertzCT`, which returns a measure of
molecular complexity. In addition to this, we will also count the
number of rings of size less than 5, as an additional measure of
complexity

.. code-block:: python

    from rdkit.Chem.GraphDescriptors import BertzCT

    def get_complexity(molecule):
        num_bad_rings = sum(
            1 for ring in rdkit.GetSymmSSSR(molecule) if len(ring) < 5
        )
        # Multiply by 10 and raise to the power of 2 to increase the
        # penalty for having many small rings. These numbers were
        # chosen by trial and error, so do don't worry about them
        # too much.
        return BertzCT(molecule) + 10*num_bad_rings**2

Now we can combine the rigidity and complexity into a single fitness
value, that the EA can optimize. There are multiple way to do this, but
an easy thing to do is just divide the rigidity by the complexity

.. code-block:: python

    def get_fitness_value(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        # Multiply by 100 just to scale the values up a bit, which
        # makes for nicer plots later.
        return 100*(
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
molecules. In this, example we will use two files from
https://github.com/lukasturcani/basic_ea, ``bromos.txt`` and
``fluoros.txt``. Each file contains the SMILES strings of buildings
blocks, holding the respective functional groups. The building
blocks in these files are randomly generated molecular graphs.
We can define a function which will load the building blocks from
these files

.. code-block:: python

    def get_building_blocks(path, functional_group_factory):
        with open(path, 'r') as f:
            content = f.readlines()

        for smiles in content:
            molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
            molecule.AddConformer(
                conf=rdkit.Conformer(molecule.GetNumAtoms()),
            )
            rdkit.Kekulize(molecule)
            building_block = stk.BuildingBlock.init_from_rdkit_mol(
                molecule=molecule,
                functional_groups=[functional_group_factory],
            )
            yield building_block.with_position_matrix(
                position_matrix=get_position_matrix(building_block),
            )


    def get_position_matrix(molecule):
        generator = np.random.RandomState(4)
        position_matrix = generator.uniform(
            low=-500,
            high=500,
            size=(molecule.get_num_atoms(), 3),
        )
        molecule = molecule.with_position_matrix(position_matrix)
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit.Compute2DCoords(rdkit_molecule)
        try:
            rdkit.MMFFOptimizeMolecule(rdkit_molecule)
        except Exception:
            pass
        return rdkit_molecule.GetConformer().GetPositions()

Once these functions are defined, we can use :func:`get_building_block`
to generate our building blocks

.. code-block:: python

    import pathlib

    fluoros = tuple(get_building_blocks(
        # Assume that fluoros.txt is in the same folder as this
        # code.
        path=pathlib.Path(__file__).parent / 'fluoros.txt',
        functional_group_factory=stk.FluoroFactory(),
    ))
    bromos = tuple(get_building_blocks(
        # Assume that bromos.txt is in the same folder as this
        # code.
        path=pathlib.Path(__file__).parent / 'bromos.txt',
        functional_group_factory=stk.BromoFactory(),
    ))

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


    initial_population = tuple(get_initial_population(fluoros, bromos))

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
            # Substitutes a building block with a fluoro group with
            # a random building block in fluoros.
            stk.RandomBuildingBlock(
                building_blocks=fluoros,
                is_replaceable=is_fluoro,
                random_seed=generator.randint(0, 1000),
            ),
            # Substitutes a building block with a fluoro group with
            # a similar building block in fluoros.
            stk.SimilarBuildingBlock(
                building_blocks=fluoros,
                is_replaceable=is_fluoro,
                random_seed=generator.randint(0, 1000),
            ),
            # Substitutes a building block with a bromo group with
            # a random building block in bromos.
            stk.RandomBuildingBlock(
                building_blocks=bromos,
                is_replaceable=is_bromo,
                random_seed=generator.randint(0, 1000),
            ),
            # Substitutes a building block with a bromo group with
            # a similar building block in bromos.
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
    )

Defining a Database
-------------------

The last thing we need to do is define the database. The default
database of :mod:`stk` is MongoDB, which can be used with
:class:`.ConstructedMoleculeMongoDb`. Before using this class, you
need to have a MongoDB database, if you want to install one locally
you can see how here__.

__ https://docs.mongodb.com/manual/installation/

Note that this is easy to do, and well worth the minimal effort it
requires to setup. Obviously, if you really don't want to
use the database, you do not have to create it, and you can remove
references to it in your EA loop.

Assuming everything is setup, we can create our database instance

.. code-block:: python

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
    from rdkit import RDLogger
    import pymongo
    import numpy as np
    import itertools as it
    import logging
    import pathlib

    rdkit_logger = RDLogger.logger()
    rdkit_logger.setLevel(RDLogger.CRITICAL)
    logger = logging.getLogger(__name__)


    def get_building_blocks(path, functional_group_factory):
        with open(path, 'r') as f:
            content = f.readlines()

        for smiles in content:
            molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
            molecule.AddConformer(
                conf=rdkit.Conformer(molecule.GetNumAtoms()),
            )
            rdkit.Kekulize(molecule)
            building_block = stk.BuildingBlock.init_from_rdkit_mol(
                molecule=molecule,
                functional_groups=[functional_group_factory],
            )
            yield building_block.with_position_matrix(
                position_matrix=get_position_matrix(building_block),
            )


    def get_position_matrix(molecule):
        generator = np.random.RandomState(4)
        position_matrix = generator.uniform(
            low=-500,
            high=500,
            size=(molecule.get_num_atoms(), 3),
        )
        molecule = molecule.with_position_matrix(position_matrix)
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit.Compute2DCoords(rdkit_molecule)
        try:
            rdkit.MMFFOptimizeMolecule(rdkit_molecule)
        except Exception:
            pass
        return rdkit_molecule.GetConformer().GetPositions()


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
        num_rotatable_bonds = rdkit.CalcNumRotatableBonds(molecule)
        # Add 1 to the denominator to prevent division by 0.
        return 1 / (num_rotatable_bonds + 1)


    def get_complexity(molecule):
        num_bad_rings = sum(
            1 for ring in rdkit.GetSymmSSSR(molecule) if len(ring) < 5
        )
        return BertzCT(molecule) + 10*num_bad_rings**2


    def get_fitness_value(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return 100*(
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


    def write(molecule, path):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        rdkit_molecule = rdkit.RemoveHs(rdkit_molecule)
        building_block = stk.BuildingBlock.init_from_rdkit_mol(
            molecule=rdkit_molecule,
        )
        building_block.with_position_matrix(
            position_matrix=get_position_matrix(building_block),
        ).write(path)


    def main():
        logging.basicConfig(level=logging.INFO)

        # Use a random seed to get reproducible results.
        random_seed = 4
        generator = np.random.RandomState(random_seed)

        logger.info('Making building blocks.')

        # Load the building block databases.
        fluoros = tuple(get_building_blocks(
            path=pathlib.Path(__file__).parent / 'fluoros.txt',
            functional_group_factory=stk.FluoroFactory(),
        ))
        bromos = tuple(get_building_blocks(
            path=pathlib.Path(__file__).parent / 'bromos.txt',
            functional_group_factory=stk.BromoFactory(),
        ))

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
        for generation in ea.get_generations(50):
            for record in generation.get_molecule_records():
                db.put(record.get_molecule())
            generations.append(generation)

        # Write the final population.
        for i, record in enumerate(generation.get_molecule_records()):
            write(record.get_molecule(), f'final_{i}.mol')

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

.. image:: https://i.imgur.com/9Difk6R.png

which shows us that the EA was pretty good at improving the fitness
value. Another thing to look at is the plot for the number of
rotatable bonds

.. image:: https://i.imgur.com/QJKTTEx.png


Clearly, our EA was able to minimize the number of rotatable
bonds to a low value across all members of the population.

We can also compare the molecules in the initial population

.. image:: https://i.imgur.com/C9Gisxf.png

to those in the final population

.. image:: https://i.imgur.com/5lq42FZ.png

where the hydrogen atoms have been left out for clarity. When
considering that these were chosen out of a search space of 1,000,000
randomly constructed molecular graphs, they don't look that bad, though
you will probably want to a better measure of synthetic accessibility
in your own EAs.

Next, you can read the intermediate tutorial, which will show you
some additional customization options for the EA.
