.. stk documentation master file, created by
   sphinx-quickstart on Mon Nov 13 11:15:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ``stk``'s documentation!
===================================

GitHub: https://www.github.com/lukasturcani/stk

Slack: https://t.co/LCPzWhvsVO

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Overview
--------

``stk`` is a Python library for building, manipulating, analyzing and
automatic design of molecules, including a genetic algorithm.

For quick navigation through the modules use :ref:`modindex`.
For a closer look at the genetic algorithm, start with
`Automated Molecular Design with Genetic Algorithms`_.

Basic Examples
--------------

Linear Polymer
..............

The core function of ``stk`` is to assemble molecules. Here is an example
of how a simple, linear polymer can be built. Starting with three monomers,
such as

.. image:: figures/monomers.png

each can be loaded into a :class:`.StructUnit2` object from a molecular
structure file:

.. code-block:: python

    monomer1 = StructUnit2('monomer1.mol', ['bromine'])
    monomer2 = StructUnit2('monomer2.mol', ['bromine'])
    monomer3 = StructUnit2('monomer3.mol', ['bromine'])

The first argument is the path to the structure file and the second
argument specifies the reactive functional groups of the monomer.

To assemble a polymer only a single line of code is required:

.. code-block:: python

    polymer = Polymer([monomer1, monomer2, monomer3], Linear('ABC', [0, 1, 0.5], n=3))

Simply create a :class:`.Polymer` object by giving it a list of
monomers and a topology object, in this
case :class:`.Linear`. The topology object defines the the structure of
the polymer being assembled. The repeating unit is ``'ABC'``, the orientation
of the each monomer along the chain is parallel, anti-parallel and random, respectively, and
the number of repeating units is ``3``.

The assembled polymer can be written to a file

.. code-block:: python

    polymer.write('polymer.mol')

and will look like this.

.. image:: figures/polymer.png

Notice that the functional group has disappeared and been replaced by
new bonds between the monomers. The new bonds seem a little stretched,
so we can optimize the structure using an :class:`.Optimizer`
defined in :mod:`~.calculator.optimization`, in this case
:class:`.MMFF`

.. code-block:: python

    mmff = MMFF()
    mmff.optimize(polymer)

Again, the polymer can be written to a file

.. code-block:: python

    polymer.write('polymer_opt.mol')

and viewed

.. image:: figures/polymer_opt.png

``stk`` also allows you to save the Python objects themselves in the JSON
format.

.. code-block:: python

    polymer.dump('polymer.json')

This allows you to restore a :class:`.Molecule` object from a previous session.

.. code-block:: python

    polymer = Molecule.load('polymer.json')
    polymer # < Polymer NOT Molecule object >

:meth:`.Molecule.load` allows you to load any dumped ``stk`` JSON object
regardless of its class and it will returned an object of the correct
class to you.

Molecular Cages
...............

Molecular cages are relatively
exotic molecules that look like, yes, cages. Here is an example:

.. image:: figures/molecular_cage.png

Despite their apparent complexity, assembling a molecular cage is
extremely straightforward. In fact, it is done in exactly the same way
as a polymer.

First we define the building blocks of the cage:

.. code-block:: python

    bb1 = StructUnit2('bb1.mol', ['amine'])
    bb2 = StructUnit3('bb2.mol', ['aldehyde'])

Here is what they look like:

.. image:: figures/cage_building_blocks.png

Notice a slight difference, while the first building building block still
uses the class :class:`.StructUnit2`, the second uses the class
:class:`.StructUnit3`. The reason is that the first building block has
2 functional groups while the second has 3 functional groups. Each class
defines a slightly different set of operations for manipulating the
positions of the building blocks when assembling the cage. This is
important so that the building blocks are placed exactly how we want them
when constructing a molecule. However, all of this happens behind the scenes
and a cage can be constructed through a simple one-liner:

.. code-block:: python

    cage = Cage([bb1, bb2], FourPlusSix())

Notice that this is exactly the same as the polymer example. To generate
a cage, we simply create a :class:`.Cage` object. We initialize it with
a list of building blocks and provide a topology instance, in this case
:class:`.FourPlusSix`. Unlike the :class:`.Linear` class, the :class:`.FourPlusSix`
does not require any additional arguments. Our assembled cage looks like
this

.. image:: figures/cage.png

If we want to create a cage with a different topology but using the same
building blocks, we provide a different topology instance:

.. code-block:: python

    cage2 = Cage([bb1, bb2], EightPlusTwelve())

The assembled cage looks like a cube:

.. image:: figures/cage2.png

Notice that the building blocks are the same, only the shape has changed.
This is because a different topology instance was provided during initialization.

Here is a third example:

.. code-block:: python

    cage3 = Cage([bb1, bb2], Dodecahedron())

.. image:: figures/cage3.png

While we assembled some cages, the constructed structures are not
particularly realistic. We can optimize the geometry using an
:class:`.Optimizer`:

.. code-block:: python

    macromodel = MacroModelForceField(
        macromodel_path='/opt/schrodinger2017-4'
    )
    macromodel.optimize(cage)
    macromodel.optimize(cage2)
    macromodel.optimize(cage3)

.. image:: figures/cages_opt.png

In this case the optimizer :class:`.MacroModelForceField` was used.
We could have used :class:`.MMFF` again but chances are the structures
would have been optimized quite poorly. The
:class:`.MacroModelForceField` optimizer
requires a valid ``MacroModel`` installation with a working license.
The argument ``'/opt/schrodinger2017-4'`` is the path to the installation.

Note that there are many more cage topologies available (14+), which
can be found by looking in :mod:`stk.molecular.topologies.cage`.

.. image:: figures/cages_two_plus_three.png

.. image:: figures/cages_two_plus_four.png

.. image:: figures/cages_three_plus_four.png

.. image:: figures/cages_three_plus_three.png

The topologies are organized into submodules based on the building blocks
required to build them. For example, all topologies in the
:mod:`stk.molecular.topologies.cage.two_plus_three` submodule are composed of two and three functionalized building
blocks, all cages in the :mod:`stk.molecular.topologies.cage.two_plus_four` submodule are composed of two and
four functionalized building blocks, cages in :mod:`stk.molecular.topologies.cage.three_plus_three` are composed
of three functionalized building blocks and so on.

All cage topologies also support being built from more than two building
blocks, to produce mixed or multi-component cages.

.. image:: figures/advanced_cage_assembly/multi_cage.png

In addition, cage topologies support a number of optional arguments which allow the many
possible structural isomers to be easily assembled. These topics are
discussed in :doc:`advanced_cage_building`.

Covalent Organic Frameworks
...........................

Just like the other molecules, covalent organic frameworks (COFs) are easy to
construct:

.. code-block:: python

    bb1 = StructUnit2('cof_bb1.mol', ['amine'])
    bb2 = StructUnit3('cof_bb2.mol', ['aldehyde']))
    cof = Periodic([bb1, bb2], Honeycomb())

Where the buliding blocks are:

.. image:: figures/cof_bbs.png

The same pattern is used. First building blocks objects are created using
:class:`.StructUnit3` instances. Then a molecule is assembled by creating
an instance of its class, in this case :class:`.Periodic`. The molecule
being assembled is provided with the building blocks and the topology,
in this case :class:`.Honeycomb`.

Because COFs are periodic structures, if we want to make a finite size
molecule we have to create an "island":

.. code-block:: python

    cof.island([5, 5, 1])

Here we created create a 5 by 5 by 1 grid of the periodic unit cell.
This is a cut-out of our structure:

.. image:: figures/honey.png

Other COF topologies are available in :mod:`.topologies.cof`. For example:

:class:`.Hexagonal`

.. image:: figures/hex.png

:class:`.Square`

.. image:: figures/square.png

:class:`.Kagome`

.. image:: figures/kagome.png

``stk`` also gives tools to build a large number of structural isomers
of each COF. This is done analogously to the organic cage case so reading
:doc:`advanced_cage_building` is recommended.

Other Materials
...............

Other materials defined by ``stk`` can be found in
:mod:`.molecular.topologies`. In addition to this, ``stk`` was designed
to be easy extend to other classes of molecules. For a guide on how this
can be done see, :ref:`extending stk`.


Calculators
-----------

``stk`` makes extensive use of calculator objects. Calculators are a
very general pattern for performing operations on molecules. All
calculators are objects with a special method, used to perform the
calculation. You have already been introduced to one type of calculator,
the :class:`.Optimizer`. An optimizer has a special method called
:meth:`~.Optimizer.optimize`, which optimizes a molecule. Here is
another example

.. code-block:: python

    mol = StructUnit(...)
    polymer = Polymer(...)

    optimizer1 = MMFF()

    # Optimize mol with mmff.
    optimizer1.optimize(mol)
    # Optimize polymer conformer 2 with mmff.
    optimizer2.optimize(polymer, conformer=2)

Other types of calculators will have different special methods used to
perform calculations on molecules. For example, an
:class:`.EnergyCalculator` will define an
:meth:`~.EnergyCalculator.energy` method. This method calculates the
energy of the molecule

.. code-block:: python

    mmff_energy = MMFFEnergy()
    uff_energy = UFFEnergy()

    # Calculate the energy of mol conformer 1 with the mmff force field.
    mol_mmff_energy = mmff_energy.energy(mol, conformer=1)

    # Calculate the energy of polymer conformer 10 with the uff force field.
    polymer_uff_energy = uff_energy.energy(polymer, conformer=10)

All genetic algorithm operations are also implemented through
calculators, take for example a :class:`.Mutator` such as
:class:`RandomTopology`

.. code-block:: python

    topologies = [EightPlusTwelve(), FourPlusSix(), Dodecahedron()]
    random_top = RandomTopology(topologies)

    # Assume bb1 and bb2 are StructUnit2 and StructUnit3 objects.
    cage = Cage([bb1, bb2], SixPlusNine())

    # Perform a mutation.
    mutant = random_top.mutate(cage)

A :class:`.Mutator` is a calculator which implements a
:meth:`~.Mutator.mutate` method. This method takes a molecule and
returns a mutant, *i.e.* a modified version, of that molecule.
In this example, the ``random_top`` is a :class:`.Mutator`, which replaces
the topology of a molecule with a new one, in order to generate the mutant.
In this case, the new topologies would be one of :class:`.EightPlusTwelve`,
:class:`.FourPlusSix` or :class:`.Dodecahedron`.

Calculators often support additional options to modify their behaviour.
For example, calculator of type :class:`.Optimizer` or
:class:`.EnergyCalculator` support caching. This means that if the
same molecule and conformer is supplied to the calculator, it will not
run the optimization or energy calculation again, it will return the
previously calculated value

.. code-block:: python

    # Example of caching with optimizers.
    mmff = MMFF()
    caching_mmff = MMFF(use_cache=True)

    mol = StructUnit(...)

    # Performs an optimization on mol.
    mmff.optimize(mol)
    # Performs a second optimization on mol.
    mmff.optimize(mol)

    # Performs an optimization on mol.
    cached_mmff.optimize(mol)
    # Does not peform an optimization on mol.
    cached_mmff.optimize(mol)

    # Make a non-caching and a caching uff energy calculator.
    uff_energy = UFF()
    caching_uff_energy = UFF(use_cache=True)

    # Calculate the energy twice.
    uff_energy.energy(mol)
    uff_energy.energy(mol)

    # Calculate the energy.
    caching_uff_energy.energy(mol)
    # Does not calculate the energy again, returns the previous value.
    caching_uff_energy.energy(mol)

More information on each type of calculator can be seen by visiting the
module dedicated to that calculator, :mod:`.calculators`.

Combining Calculators
.....................

Calculators can be combined to create complex behaviour on the fly.
For example, we may wish to make a 2 step optimization process. First,
we perform an optimization using the MMFF force field and the run
a MacroModel molecular dynamics conformer search with
:class:`.MacroModelMD`. The obvious way to do
this is to run the two optimizers in sequence

.. code-block:: python

    mmff = MMFF()
    macromodel = MacroModelMD(macromodel_path='/opt/schrodinger2017-4')

    mmff.optimize(mol)
    macromodel.optimize(mol)

However, there is a better way! We can use an optimizer called
:class:`.OptimizerSequence`. The :meth:`~OptimizerSequence.optimize`
method of this optimizer calls the :class:`~Optimizer.optimize` methods
of the optimizers it was initialized with

.. code-block:: python

    opt_sequence = OptimizerSequence(mmff, macromodel)
    # Optimize with mmff and then with macromodel.
    opt_sequence.optimize(mol)

This pattern is quite common and powerful. For example, we can take
three different :class:`.Mutator` objects. Each of these defines
a different :meth:`~.Mutator.mutate` method. We want to apply one of
these mutations at random. We can simply use :class:`.RandomMutation`

.. code-block:: python

    random_bb = RandomBuildingBlock(...)
    similar_bb = SimilarBulidingBlock(...)
    random_topology = RandomTopology(...)

    random_mutation = RandomMutation(
        random_bb,
        similar_bb,
        random_topology
    )

    # Use one of the mutate() methods of random_bb, similar_bb and
    # random_topology at random.
    mutant1 = random_mutation.mutate(mol)
    # The next call use a different mutation to the call above.
    mutant2 = random_mutation.mutate(mol)

The :meth:`.RandomMutation.mutate` method randomly selects a
:class:`.Mutator` it was initialized with to carry out the mutation
on its behalf.

Making New Calculators
......................

New calculators can be added very simply and they can be defined in
user code. Specific guidelines are defined in the module
dedicated to each type of calculator in :mod:`.calculators`.

A simple example of adding a new calculator in user code.

.. code-block:: python

    class NewEnergyCalculator(EnergyCalculator):
        def energy(self, mol, conformer=-1):
            return 15

    energy_calc = NewEnergyCalculator()
    # Calculate the energy with the new calculator.
    energy_calc.energy(mol)


Note that you can also modify the behaviour of existing calculators in
your code, some common patterns are described in
:doc:`modifying_calculators`.

Saving Calculator Results
.........................

If you want to save the results which your calculator found you should
set ``use_cache=True`` and then you can simply use :module:`pickle`
to dump and load the calculator object. You can inspect the saved
values by looking at the :attr:`cache` attribute.

Calculating Molecular Properties
................................

``stk`` provides a variety of methods to calculate molecular properties.
What methods can be used depends on what kind of molecule object is created.
Note that all molecules can use methods defined in :class:`.Molecule`.

Here are some examples:

.. code-block:: python

    # Calculate the maximum diameter of of a molecule.
    mol.max_diameter()
    # Calculate the cavity size of a cage molecule.
    cage.cavity_size()
    # Calculate the mean RMSD between the building blocks original
    # structure and their structure inside an assembled macromolecule.
    macro_mol.bb_distortion()
    # Get the center of mass of a molecule.
    mol.center_of_mass()
    # Get the atomic symbol of atom with id of 13.
    mol.atom_symbol(13)

Geometric Manipulations
.......................

In addition to molecular property calculation, ``stk`` provides tools to
rotate and translate molecules. These tools are particularly useful when
defining the assembly process of a new class of molecules.

.. code-block:: python

    # Change the position of a molecule.
    mol.set_position([1, 2, 3])
    # Get a matrix holding the position of every atom in the molecule.
    mol.position_matrix()
    # Use a matrix to set the position of every atom in the molecule.
    mol.set_position_from_matrix(some_matrix)
    # Rotate the molecule along by pi radians about the vector (1, 1, 3)
    mol.rotate(np.pi, [1, 1, 2])

Dealing with Multiple Conformers
--------------------------------

Every :class:`.Molecule`, be it a :class:`.StructUnit` or a :class:`.MacroMolecule`
supports multiple conformers. These are stored in the underlying
:mod:`rdkit` object held in :class:`.Molecule.mol`.

Adding a new conformer to a :class:`.StructUnit` is simple

.. code-block:: python

    bb1 = StructUnit2('bb1_conf1.mol', ['amine']) # Loads molecule into conformer 0.
    bb1.update_from_mol('bb1_conf2.mol', conformer=1) # Load into conformer 1.

Lets take a look

    .. image:: figures/bb1_confs.png

Well, one is clearly better than the other. At least they are easy to
recognize.

Lets try this with a another building block as well

.. code-block:: python

    bb2 = StructUnit3('bb2_conf1.mol', ['aldehyde'])
    bb2 = bb2.update_from_mol('bb2_conf2.mol')

.. image:: figures/bb2_confs.png

Also pretty distinguishable. This is good, it makes it easy to show
that you can pick which conformers to use for assembly with :class:`.MacroMolecule`

.. code-block:: python

    cage = Cage([bb1, bb2], FourPlusSix(), bb_conformers=[0, 0])
    cage.add_conformer([0, 1])
    cage.add_conformer([1, 0])
    cage.add_conformer([1, 1])

These are the conformers we produced!

.. image:: figures/1.png

.. image:: figures/2.png

.. image:: figures/3.png

.. image:: figures/4.png

``stk`` will use whichever :class:`.StructUnit` conformers you want, when
constructing conformers of the :class:`.MacroMolecule` itself.

Most methods and functions in ``stk`` also support conformers.

.. code-block:: python

    mol.position_matrix(conformer=1)
    mol.energy.rdkit('uff', conformer=0)
    mol.cavity_size(conformer=1)

    # And so on ...

Dealing with Multiple Molecules
...............................

When batches of molecules are created, it is often desirable to optimize
them all at once. By placing the molecules in a :class:`.Population`
instance, all molecules can be optimized in parallel.

.. code-block:: python

    mmff = MMFF()
    pop = Population(cage, cage2, cage3)
    pop.optimize(mmff)


In addition to this, the :class:`.Population` class provides some handy
tools to assemble large amounts of molecules at a time. For example if
we want to create every possible cage from a set of building blocks and
topologies:

.. code-block:: python

    bbs1 = [StructUnit2(path, ['aldehyde']) for path in ('1.mol', '2.mol', '3.mol')]
    bbs2 = [StructUnit3(path, ['amine']) for path in ('4.mol', '5.mol', '6.mol')]
    # Create 18 Cage molecules.
    pop2 = Population.init_all(Cage, [bbs1, bbs2], [FourPlusSix(), EightPlusTwelve()])

Or if we want to select building blocks at random and create 5 cages:

.. code-block:: python

    pop3 = Population.init_random(Cage, [bbs1, bbs2], [FourPlusSix(), EightPlusTwelve()], 5)

The population can also be used to calculate statistics across all
molecules. For example, if you want to know the average cavity size of
your cages:

.. code-block:: python

    pop.mean(lambda x: x.cavity_size())

Finally, we can use the population to write structure files of
all the molecules to a folder

.. code-block:: python

    pop.write('dirname')

or we can dump the population to a JSON file

.. code-block:: python

    pop.dump('pop.json')

and recover the ``stk`` objects later

.. code-block:: python

    pop.load('pop.json', Molecule.from_dict)
    pop[0] # <StructUnit2 at 0x01283>
    pop[1] # <Periodic at 0x12498>

Automated Molecular Design with Genetic Algorithms
..................................................

``stk`` includes a genetic algorithm which
can be used to evolve molecules that fulfil user defined design criteria.
The genetic algorithm can be run from the command line using::

    $ python -m stk.ga input_file.py

The input file is a simple python script which defines the
calculators the genetic algorithm should use, as well as some optional
parameters.

The genetic algorithm automatically works with any molecules that ``stk``
uses, both :class:`.StructUnit` and :class:`.MacroMolecule` objects.

Take for example the following input file which runs the GA on polymers
and selects building blocks which have the most atoms.

.. code-block:: python

    # big_monomers.py

    # ################################################################
    # Import stuff.
    # ################################################################

    import stk

    # ################################################################
    # Number of CPU cores to use.
    # ################################################################

    processes = 1

    # ################################################################
    # Disable debug population dumps.
    # ################################################################

    pop_dumps = False

    # ################################################################
    # Population size.
    # ################################################################

    pop_size = 20

    # ################################################################
    # Number of generations to create.
    # ################################################################

    num_generations = 50

    # ################################################################
    # Number of mutation operations to perform each generation.
    # ################################################################

    num_mutations = 5

    # ################################################################
    # Number of crossover operations to perform each generation.
    # ################################################################

    num_crossovers = 3

    # ################################################################
    # Population initialization function.
    # ################################################################

    # First create the StructUnit objects of the monomers. In
    # this simple example, they are simple hydro-carbon chains terminated
    # with either two amine or two aldehyde functional groups.
    amines = []
    aldehydes = []
    for i in range(1, 50):
        # There are 2 different backbones, alkane, alkene.
        backbones = ['C'*i, '='.join('C'*i)]

        for backbone in backbones:
            amine_smiles = 'N' + backbone + 'N'
            amine_monomer = stk.StructUnit2.smiles_init(
                                             smiles=amine_smiles,
                                             functional_groups=['amine'])
            amines.append(amine_monomer)

            aldehyde_smiles = 'O=C' + backbone + 'C=O'
            aldehyde_monomer = stk.StructUnit2.smiles_init(
                                            smiles=aldehyde_smiles,
                                            functional_groups=['aldehyde'])
            aldehydes.append(aldehyde_monomer)

    # List of the possible topologies the polymers can have. In this
    # simple example a single topology will suffice.
    polymer_topologies = [stk.Linear(repeating_unit='AB',
                                     orientation=[0, 0],
                                     n=1)]

    # Finally, define the function which is called to create the
    # initial GA population. The "NAME" key holds the name of an
    # initialization function of the GAPopulation class. The remaining
    # key-value pairs are additional parameters for that function.
    # Default parameters can be specified if desired.
    init_func = {'NAME': 'init_random',
                 'macromol_class': stk.Polymer,
                 'building_blocks': [amines, aldehydes],
                 'topologies': polymer_topologies}

    # ################################################################
    # Selection function for selecting the next generation.
    # ################################################################

    # The "NAME" key holds the name of a method in the Selection class.
    # The remaining key-value pairs are additional parameters for that
    # method. The "population" parameter is not specified and any
    # default parameters can be specified if desired.
    generational_select_func = {'NAME': 'stochastic_sampling',
                                'use_rank': True}

    # ################################################################
    # Selection function for selecting parents.
    # ################################################################

    # The "NAME" key holds the name of a method in the Selection class.
    # The remaining key-value pairs are additional parameters for that
    # method. The "population" parameter is not specified and any
    # default parameters can be specified if desired.
    crossover_select_func = {'NAME': 'crossover_roulette',
                             'truncation': 5}

    # ################################################################
    # Selection function for selecting molecules for mutation.
    # ################################################################

    # The "NAME" key holds the name of a method in the Selection class.
    # The remaining key-value pairs are additional parameters for that
    # method. The "population" parameter is not specified and any
    # default parameters can be specified if desired.
    mutation_select_func = {'NAME': 'stochastic_sampling',
                            'truncation': 3,
                            'duplicates': True}

    # ################################################################
    # Crossover functions.
    # ################################################################

    # In this example there only one crossover function will be used,
    # but can use multiple if we add them into the list.

    # The "NAME" key holds the name of a method in the Crossover class.
    # The "macro_mol" parameters are not specified and any default
    # parameters can be specified optionally. All other parameters
    # must be specified as key-value pairs.
    crossover_funcs = [{'NAME': 'genetic_recombination',
                        'key': lambda x: x.func_group_infos[0].name}]

    # ################################################################
    # Mutation functions.
    # ################################################################

    # The "NAME" key holds the name of a method in the Mutation class.
    # The "macro_mol" parameter is not specified and any default
    # parameters can be specified optionally. All other parameters
    # must be specified as key-value pairs.

    # mutation_func1 substitutes the amine building block in the
    # molecule with one from the list "amines", defined earlier.
    mutation_func1 = {'NAME': 'random_bb',
                      'mols': amines,
                      'key': lambda x: x.func_group_infos[0].name == 'amine'}

    # mutation_func2 substitutes the aldehyde building block in the
    # molecule with one from the list "aldehydes", defined earlier.
    mutation_func2 = {'NAME': 'random_bb',
                      'mols': aldehydes,
                      'key': lambda x: x.func_group_infos[0].name == 'aldehyde'}

    # Here we tell the GA which mutation functions we want to use.
    mutation_funcs = [mutation_func1, mutation_func2]

    # ################################################################
    # Optimization function.
    # ################################################################

    # In this example, for the sake of speed, we will not optimize
    # constructed polymers.

    # The "NAME" key holds the name of a function
    # defined in the optimization.py module. The "mol" parameter is
    # not specified and any default parameters can be specified
    # optionally. All other parameters must be specified as key-value
    # pairs.

    # Note that "do_not_optimize" is the name of a function
    # in this module. It simply returns immediately and as a result does
    # no work.
    opt_func = {'NAME': 'do_not_optimize'}

    # ################################################################
    # Fitness function.
    # ################################################################

    # The "NAME" key holds the name of a function defined in the
    # "fitness.py" module. The "macro_mol" parameter is not
    # specified and any default parameters can be specified if desired.
    # The remaining parameters must be specified as key-value pairs.
    fitness_func = {'NAME': 'building_block_atoms'}


Running the genetic algorithm with this input file::

    $ python -m stk.ga big_monomers.py

will produce the following directory structure::

    |-- stk_ga_runs
    |   |-- 0
    |   |   |-- counters
    |   |   |   |-- gen_1_crossover_counter.png
    |   |   |   |-- gen_1_mutation_counter.png
    |   |   |   |-- gen_1_selection_counter.png
    |   |   |   |-- ...
    |   |   |
    |   |   |-- final_pop
    |   |   |   |-- 150.mol
    |   |   |   |-- 2160.mol
    |   |   |   |-- 9471.mol
    |   |   |   |-- ...
    |   |   |
    |   |   |-- big_monomers.py
    |   |   |-- database.json
    |   |   |-- progress.json
    |   |   |-- errors.log
    |   |   |-- progress.log
    |   |   |-- epp.png
    |   |   |-- epp.csv
    |   |   |-- output.tgz

A glance at the evolutionary progress plot in ``epp.png`` will show us
how well our GA did.

.. image:: figures/epp.png


Running the genetic algorithm again::

    $ python -m stk.ga big_monomers.py

will add a second subfolder with the same structure::

    |-- stk_ga_runs
    |   |-- 0
    |   |   |-- counters
    |   |   |   |-- gen_1_crossover_counter.png
    |   |   |   |-- gen_1_mutation_counter.png
    |   |   |   |-- gen_1_selection_counter.png
    |   |   |   |-- ...
    |   |   |
    |   |   |-- final_pop
    |   |   |   |-- 150.mol
    |   |   |   |-- 2160.mol
    |   |   |   |-- 9471.mol
    |   |   |   |-- ...
    |   |   |
    |   |   |-- big_monomers.py
    |   |   |-- database.json
    |   |   |-- progress.json
    |   |   |-- errors.log
    |   |   |-- progress.log
    |   |   |-- epp.png
    |   |   |-- epp.csv
    |   |   |-- output.tgz
    |
    |   |-- 1
    |   |   |-- counters
    |   |   |   |-- gen_1_crossover_counter.png
    |   |   |   |-- gen_1_mutation_counter.png
    |   |   |   |-- gen_1_selection_counter.png
    |   |   |   |-- ...
    |   |   |
    |   |   |-- final_pop
    |   |   |   |-- 41.mol
    |   |   |   |-- 221.mol
    |   |   |   |-- 1391.mol
    |   |   |   |-- ...
    |   |   |
    |   |   |-- big_monomers.py
    |   |   |-- database.json
    |   |   |-- progress.json
    |   |   |-- errors.log
    |   |   |-- progress.log
    |   |   |-- epp.png
    |   |   |-- epp.csv
    |   |   |-- output.tgz

The genetic algorithm can also be run multiple times in a row::

    $ python -m stk.ga -l 5 big_monomers.py

which will run the GA 5 separate times adding 5 more subfolders to the
directory structure::

    |-- stk_ga_runs
    |   |-- 0
    |   |   |-- ...
    |   |
    |   |-- 1
    |   |   |-- ...
    |   |
    |   |-- 2
    |   |   |-- ...
    |   |
    |   |-- 3
    |   |   |-- ...
    |   |
    |   |-- 4
    |   |   |-- ...
    |   |
    |   |-- 5
    |   |   |-- ...
    |   |
    |   |-- 6
            |-- ...

The benefit of using the ``-l`` option is that the molecular cache is
not reset between each run. This means that a molecule which was constructed,
optimized and had its fitness value calculated in the first run will
not need to be re-constructed, re-optimized or have fitness value
re-calculated in any of the subsequent runs. The cached version
of the molecule will be used.

However, the molecular cache be pre-loaded even when the ``-l`` option is
not used. The input file allows the definition of a variable
:attr:`~.GAInput.databases`, for example

.. code-block:: python

    # some_inputfile.py

    databases = ['path/to/some/population1.json',
                 'path/to/some/other/population2.json']

If this variable is defined in the input file, any molecules present
in the population dump files will be loaded into the cache before the
GA starts. This means that any molecules found in these populations
will not be re-optimized, or have their fitness value re-calculated
if discovered by the GA.

The output of a single GA consists of a number of files and
directories. The ``counter`` directory holds ``.png`` files showing
how frequently a member of the population was selected for mutation,
crossover and generational selection. For example, this is a
crossover counter

.. image:: figures/counter_example.png

It shows that molecule ``168`` was selected twice for crossover, while molecules
``60``, ``170``, ``167``, ``63`` and ``153`` were selected once. The
remaining molecules did not participate in crossover operations
in that generation.

If you do want want the GA to produce these files, simply tell the
GA in the input file

.. code-block:: python

    # some_input_file.py

    counters = False

The ``final_pop`` directory holds the ``.mol`` files holding the
structures of the last generation of molecules.
The file ``big_monomers.py`` is a copy of the input file. The ``database.json``
file is a population dump file which holds every molecule produced by
the GA during the run. This can be turned off with

.. code-block:: python

    # some_input_file.py

    database_dump = False

``progress.json`` is also a population dump file. This population holds
every generation of the GA as a subpopulation. This is quite useful
if you want to analyse the output of the GA generation-wise. It
can be turned off with

.. code-block:: python

    # some_input_file.py

    progress_dump = False

``errors.log`` is a file which contains every exception and its
traceback encountered by the GA during its run.

``progress.log`` is a file which lists which molecules make up each
generation, and their respective fitness values.

``epp.csv`` is a dumped :class:`pandas.DataFrame` which holds the data used
make ``epp.png``.

``output.tgz`` is a tarred and compressed copy of the output folder for
the run.
This means if you want to share you entire run output you can just
share this file. It can be turned off with

.. code-block:: python

    # some_input_file.png

    tar_output = False

Finally, when running the GA the progress will be printed into
stderr. The message should be relatively straightforward. One thing
to note though is that during optimizations and fitness function
calculation, any molecule that has already been optimized and had its
fitness value calculated will be skipped, so that duplicate work is
avoided. Briefly an example message

::

    ======================================================================

    17:42:20 - INFO - stk.ga.mutation - Using random_bb.

    ======================================================================

show the time, the level of the message which can be, in order of
priority DEBUG, INFO, WARNING, ERROR or CRITICAL, the module where
the message originated and finally the message itself.

.. _`extending stk`:

Extending ``stk``
-----------------

Each module of ``stk`` has its own guidelines for adding new functionality.
However, in almost all cases adding new features to ``stk`` only involves
defining a simple function in the appropriate module or a method in the
appropriate class.

    * :ref:`adding macromolecules`
    * :ref:`adding topologies`
    * :ref:`adding functional groups`
    * :ref:`adding optimization functions`
    * :ref:`adding energy functions`
    * :ref:`adding mutation functions`
    * :ref:`adding crossover functions`
    * :ref:`adding fitness functions`
    * :ref:`adding normalization functions`
    * :ref:`adding selection functions`
    * :ref:`adding exit functions`


Further Reading
---------------

    * :ref:`macromolecular assembly`
    * :ref:`cof assembly`
    * :doc:`advanced_cage_building`
    * :doc:`caching`
    * https://chemrxiv.org/articles/STK_A_Python_Toolkit_for_Supramolecular_Assembly/6127826

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
