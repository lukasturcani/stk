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
the functional groups which undergo the reaction specified.

.. code-block:: python

    import stk

    # React the amine functional group during construction.
    bb1 = stk.BuildingBlock('NCCN', ['amine'])
    # React the aldehyde functional group during construction.
    bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2],
        topology_graph=stk.cage.FourPlusSix()
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



Optimizing Molecules
====================


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
