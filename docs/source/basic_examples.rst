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
available initialization methods.


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
construction. Optimization is performed by optimizers, which implement
the :meth:`~.Optimizer.optimize` method. Check here_ for a list of
available optimizers.

.. _here: stk.calculators.optimizers.html

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

    caching_mff = stk.MMFF(use_cache=True)
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



Using Multiple Functional Groups
================================


Removing Extra Functional Groups
================================


Constructing Molecules in Bulk
=============================


Optimizing Molecules in Bulk
============================


Saving Molecules
================


Caching
=======
