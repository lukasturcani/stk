==============
Basic Examples
==============


Creating Building Blocks
========================

There are a number of ways to create :class:`.BuildingBlock` molecules,
but the two most common ways are direct creation through SMILES
strings and loading them from molecule structure files.


.. code-block:: python

    import stk

    bb1 = stk.BuildingBlock('NCCCN')
    bb2 = stk.BuildingBlock.init_from_file('path/to/file.mol')


Look at the documentation of :class:`.BuildingBlock` to see other
available initialization methods. They are easy to find, as they all
begin with :mod:`init`.

Specifying Multiple Functional Groups
=====================================

When you create a :class:`.BuildingBlock`, you need also specify
which atoms are modified during construction of a
:class:`.ConstructedMolecule`. This is achieved by providing the
:class:`.BuildingBlock` with
:mod:`~.functional_groups.functional_group` instances. To save you
the pain of creating function groups one by one, you can use a
:mod:`~.functional_group_factory`. If you have a building block
with bromo groups, and you want the bromo groups to be modified
during construction, you would use a :class:`.BromoFactory`

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

The ``bb`` in the example above, would have two :class:`.Bromo`
functional groups. Similarly, if you have a building block with
aldehyde groups


.. code-block:: python

    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactor()])

In this example, ``bb2`` will have two :class:`.Aldehyde` functional
groups. Finally, if you have both aldehdye and bromo groups on a
molecule, and you want both to be modified during construction,
you would use both of the factories

.. code-block:: python

    bb3 = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFacotry(), stk.BromoFactory()],
    )

In the example above, ``b3`` has one :class:`.Bromo` and one
:class:`.Aldehyde` functional group.


Constructing Molecules
======================

To construct molecules, you need to create a new
:class:`.ConstructedMolecule`. The required input consists of
building block molecules and a topology graph. :class:`.BuildingBlock`
molecules which undergo reactions during construction need to have
the functional groups which undergo them specified

.. code-block:: python

    import stk

    # React the amine functional group during construction.
    bb1 = stk.BuildingBlock('NCCN', ['amine'])
    # React the aldehyde functional group during construction.
    bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
    # Build a bunch of cages, each has the same building blocks but
    # different topology graphs.
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2],
        topology_graph=stk.cage.FourPlusSix()
    )
    cage2 = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2],
        topology_graph=stk.cage.EightPlusTwelve()
    )
    cage3 = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2],
        topology_graph=stk.cage.TwentyPlusThirty()
    )

Created :class:`.ConstructedMolecule` instances can be used as building
blocks for others

.. code-block:: python

    guest = stk.BuildingBlock('BrBr')
    host_guest_complex = stk.ConstructedMolecule(
        building_blocks=[cage, guest],
        topology_graph=stk.host_guest.Complex()
    )

Look at the documentation of the topology graph to see any requirements
associated with building :class:`.ConstructedMolecule` instances
with it as well as usage examples.


Using Multiple Building Blocks
==============================


You can have as many different building blocks in any
:class:`.ConstructedMolecule` as there are vertices, however you will
often have to assign them manually. This is because there are often
multiple equally valid ways to place the building blocks on the
vertrices and ``stk`` has no way to figure out which one you want.

.. code-block:: python

    import stk

    bb1 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
    bb2 = stk.BuildingBlock('O=CC(Cl)(C=O)C=O', ['aldehyde'])
    bb3 = stk.BuildingBlock('NCCN', ['amine'])
    bb4 = stk.BuildingBlock('NCC(Cl)N', ['amine'])
    bb5 = stk.BuildingBlock('NCCCCN', ['amine'])

    tetrahedron = stk.cage.FourPlusSix()
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2, bb3, bb4, bb5],
        topology_graph=tetrahedron,
        building_block_vertices={
            bb1: tetrahedron.vertices[:2],
            bb2: tetrahedron.vertices[2:4],
            bb3: tetrahedron.vertices[4:5],
            bb4: tetrahedron.vertices[5:6],
            bb5: tetrahedron.vertices[6:]
        }
    )


Optimizing Molecular Structures
===============================

Molecules used by ``stk`` can be optimized both before and after
construction. Optimization is performed by optimizers_, which implement
the :meth:`~.Optimizer.optimize` method.

.. _optimizers: stk.calculators.optimization.optimizers.html

.. code-block:: python

    import stk

    # Optimize with the MMFF force field.
    mmff = stk.MMFF()
    bb = stk.BuildingBlock('BrCCBr', ['bromine'])
    mmff.optimize(bb)

    polymer = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 15)
    )
    mmff.optimize(polymer)

    # Optimize with UFF.
    uff = stk.UFF()
    uff.optimize(bb)
    uff.optimize(polymer)


All optimizers support the *use_cache* option, which prevents a
second optimization from being run on the same molecule twice

.. code-block:: python

    mmff = stk.MMFF()
    bb = stk.BuildingBlock('CCCC')
    # Runs an optimization.
    mmff.optimize(bb)
    # Run an optimization again.
    mmff.optimize(bb)

    caching_mmff = stk.MMFF(use_cache=True)
    # Runs an optimization.
    caching_mmff.optimize(bb)
    # Does not run an optimization, returns immediately.
    caching_mmff.optimize(bb)

An important optimizer to take note of is the
:class:`.OptimizerSequence`, which allows you to chain multiple
optimizers in one go. For example, you may wish to embed a molecule
with the ETKDG algorithm before running a UFF optimization

.. code-block:: python

    sequence = stk.OptimizerSequence(
        stk.ETKDG(),
        stk.UFF()
    )
    # Embed with ETKDG then do a UFF optimization.
    sequence.optimize(bb)


Calculating Molecular Energy
============================

The energy of molecules can be calculated with `energy calculators`_.

.. _`energy calculators`: stk.calculators.energy.energy_calculators.html

All energy calculators define the :meth:`~.EnergyCalculator.get_energy`
method, which is used to calculate the energy

.. code-block:: python

    import stk

    mmff = stk.MMFFEnergy()
    bb = stk.BuildingBlock('BrCCBr', ['bromine'])
    bb_energy = mmff.get_energy(bb)

    polymer = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 15)
    )
    polymer_energy = mmff.get_energy(polymer)

Much like optimizers, energy calculators support a *use_cache*
option. If this is turned on the energy calculator will not
calculate the energy of a molecule twice. Instead, if the same
molecule is passed a second time, the previous value will be returned
from memory

.. code-block:: python

    mmff = stk.MMFF()
    caching_mmff = stk.MMFFEnergy(use_cache=True)
    bb_energy = caching_mmff.get_energy(bb)
    mmff.optimize(bb)
    # bb_energy2 is equal to bb_energy even though the structure
    # changed, since the old value was returned from memory.
    bb_energy2 = caching_mmff.get_energy(bb)


Using Multiple Functional Groups
================================

Sometimes you may have a :class:`.BuildingBlock` where multiple
functional groups should be used during construction. To allow this
you can pass the names of any number of functional groups to the
:class:`.BuildingBlock`

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('ICCBr', ['bromine', 'iodine'])
    polymer = stk.ConstructedMolecule(
        bulding_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 10)
    )


Removing Extra Functional Groups
================================

Some building blocks may have too many functional groups to be used
with a topology graph. For example a :class:`.Linear` topology graph
requires two functional groups on a :class:`.BuildingBlock`.
If there are more than two it is not clear which of the functional
groups should be used to form bonds during construction, as they are
all equivalent to ``stk``. However, it is possible to resolve this
ambguity by removing the extra functional groups

.. code-block:: python

    import stk

    # This buildingb block has 5 functional groups.
    bb = stk.BuildingBlock('BrCC(Br)CC(Br)CC(Br)CC(Br)CC', ['bromine'])
    # Keep only the first 2 functional groups.
    bb.func_groups = bb.func_groups[:2]
    # bb can now be used to construct linear polymers.
    polymer = stk.ConstructedMolecule(
        bulding_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 10)
    )

Constructing Molecules in Bulk
==============================

If you have a bunch of building blocks and topologies and you want
to create all possible :class:`.ConstructedMolecule` instances, use
:meth:`.Population.init_all`.

If you want to create a random population of a finite size, use
:meth:`.Population.init_random`.


Optimizing Molecules in Bulk
============================

To optimize a bunch of molecules at the same time, in parallel across
your CPU cores, place them into a :class:`.Population`

.. code-block:: python

    import stk

    bb1 = stk.BuildingBlock('NCCCN')
    bb2 = stk.BuildingBlock('BrCCBr', ['bromine'])
    polymer = stk.ConstructedMolecule(
        building_blocks=[bb2],
        topology_graph=stk.polymer.Linear('A', 12)
    )
    pop = stk.Population(bb1, bb2, polymer)

    mmff = stk.MMFF()
    # Optimization done on all molecules in parallel.
    pop.optimize(mmff)


Converting between :mod:`rdkit` Molecules
=========================================

:mod:`rdkit` is a popular and powerful cheminformatics library which
provides many useful tools for molecular analysis. To take advantage
of it ``stk`` can convert its molecules into those of :mod:`rdkit`

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('ICCBr', ['bromine', 'iodine'])
    bb_rdkit = bb.to_rdkit_mol()

    polymer = stk.ConstructedMolecule(
        bulding_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 10)
    )
    polymer_rdkit = polymer.to_rdkit_mol()

``stk`` also allows you to update the structure of molecules from
:mod:`rdkit` molecules

.. code-block:: python

    # Update structure of bb to match structure of bb_rdkit.
    bb.update_from_rdkit_mol(bb_rdkit)
    # Update structure of polymer to match structure of polymer_rdkit.
    polymer.update_from_rdkit_mol(polymer_rdkit)


Finally, ``stk`` allows initialization of new building blocks from
:mod:`rdkit` molecules directly

.. code-block:: python

    bb2 = stk.BuildingBlock.init_from_rdkit_mol(bb_rdkit)


Saving Molecules
================

The simplest way to save molecules is to write them to a file

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('ICCBr', ['bromine', 'iodine'])
    bb.write('bb.mol')

    polymer = stk.ConstructedMolecule(
        bulding_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 10)
    )
    polymer.write('polymer.mol')

However, writing to regular chemical file format is lossy. This is
because :class:`.BuildingBlock` and :class:`.ConstructedMolecule`
contain data that is not held by these files. :class:`.BuildingBlock`
holds information about which functional groups are to be used
during construction and :class:`.ConstructedMolecule` holds information
about which building blocks molecules and topology graph were used to
construct it. While a :class:`.BuildingBlock` can be initialized from
chemical file formats a :class:`.ConstructedMolecule` cannot, as there
is no way to recover the information regarding building blocks and
the topology graph.

If you wish to be able to fully recover :class:`.BuildingBlock` and
:class:`.ConstructedMolecule` instances, you can write them as
JSON files

.. code-block:: python

    bb.dump('bb.dump')
    polymer.dump('polymer.dump')

These can then be loaded in a later session

.. code-block:: python

    recovered_bb = stk.Molecule.load('bb.dump')
    recovered_polymer = stk.Molecule.load('polymer.dump')

Using the :meth:`.Molecule.load` method means you do not have to know
if the molecule in the file was a :class:`.BuildingBlock` or a
:class:`.ConstructedMolecule`, the correct class be determined
for you automatically.

You can write and dump molecules in bulk with a :class:`.Population`

.. code-block:: python

    pop = stk.Population(bb, polymer)
    pop.write('my_pop')
    pop.dump('population.dump')

and load it too

.. code-block:: python

    recovered_pop = stk.Population.load('population.dump')


Caching
=======

A powerful but dangerous feature of ``stk`` use the molecular cache.
Every constructor for a :class:`.BuildingBlock` or a
:class:`.ConstructedMolecule` has a *use_cache* option. If
this is set to ``True`` and an identical molecule has already been
loaded or constructed by ``stk``, then ``stk`` will not create a
new instance but return the existing one from memory. In the case of
:class:`.BuildingBlock` it can mean that the memory footprint is
dramatically reduced, while in the case of
:class:`.ConstructedMolecule` it means that expensive constructions do
not have to be repeated.

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('BrCCBr', ['bromine'])
    # bb and bb2 are separate instances because when bb was created it
    # was not added to the cache.
    bb2 = stk.BuildingBlock('BrCCBr', ['bromine'], use_cache=True)
    # bb2 and bb3 are different names for the same object, as
    # both bb2 and bb3 had use_cache set to True.
    bb3 = stk.BuildingBlock('BrCCBr', ['bromine'], use_cache=True)
    # bb4 is a completely new instance.
    bb4 = stk.BuildingBlock('BrCCBr', ['bromine'])


    polymer = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 5)
    )
    # polymer2 is a separate instance.
    polymer2 = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 5),
        use_cache=True
    )
    # polymer2 and polymer3 are different names for the same object as
    # both polymer2 and polymer3 had use_cache set to True.
    # No construction was carried out when making polymer3, it was
    # returned directly from memory.
    polymer3 = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 5),
        use_cache=True
    )
    # polymer4 is a completely new instance, construction was
    # done from scratch.
    polymer4 = stk.ConstructedMolecule(
        building_blocks=[bb],
        topology_graph=stk.polymer.Linear('A', 5)
    )
