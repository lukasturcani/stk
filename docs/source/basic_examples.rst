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

The ``bb``, in the example above, would have two :class:`.Bromo`
functional groups. Similarly, if you have a building block with
aldehyde groups

.. code-block:: python

    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])

In this example, ``bb2`` will have two :class:`.Aldehyde` functional
groups. Finally, if you have both aldehyde and bromo groups on a
molecule, and you want both to be modified during construction,
you would use both of the factories

.. code-block:: python

    bb3 = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFactory(), stk.BromoFactory()],
    )

In the example above, ``bb3`` has one :class:`.Bromo` and one
:class:`.Aldehyde` functional group.

Constructing Molecules
======================

To construct molecules, you need to create a new
:class:`.ConstructedMolecule`. The required input consists of
a :class:`.TopologyGraph`, which, in turn,  requires
:class:`.BuildingBlock` instances.

.. code-block:: python

    import stk

    # React the amine functional groups during construction.
    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    # React the aldehyde functional groups during construction.
    bb2 = stk.BuildingBlock('O=CC(C=O)C=O', [stk.AldehydeFactory()])
    # Build a polymer.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=12,
        ),
    )

    # Build a longer polymer.
    longer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=23,
        ),
    )

Each topology graph requires different input parameters,
for example, organic cage topology graphs only require the
:class:`.BuildingBlock` instances.

.. code-block:: python

    cage = stk.ConstructedMolecule(stk.cage.FourPlusSix((bb1, bb2)))


Read the documentation for each kind of :class:`.TopologyGraph`, for
more examples on how to initialize it, and to see what optional
parameters you have available.

Specifying Functional Groups Individually
=========================================

If you want to be more precise about which functional groups get
created, you can provide them directly to the :class:`.BuildingBlock`.
For example, if you have multiple bromo groups on a molecule, but
you only want to use one during construction

.. code-block:: python

    import stk

    bb = stk.BuildingBlock(
        smiles='BrCCCBr',
        functional_groups=[
            stk.Bromo(
                # The number is the atom's id.
                bromine=stk.Br(0),
                atom=stk.C(1),
                # bonders are atoms which have bonds added during
                # construction.
                bonders=(stk.C(1), ),
                # deleters are atoms which are deleted during
                # construction.
                deleters=(stk.Br(0), ),
            ),
        ],
    )

When creating a :class:`.Bromo` functional group, you have to
specify things like which atoms have bonds added during construction,
and which ones are removed during construction. These are specified by
the `bonders` and `deleters` parameters, respectively. You can add
as many functional groups to :class:`.BuildingBlock` as you like
in this way, and you can mix different types of
:mod:`~.functional_groups.functional_group`. You can even mix
a :mod:`~.functional_groups.functional_group` instances with
:mod:`~.functional_group_factory` instances.

Changing Bonder and Deleter Atoms in Functional Group Factories
===============================================================

In the previous example, you saw that during creation of a
:class:`.Bromo` instance, you can specify which atoms have bonds
added during construction, and which atoms are deleted during
construction. You might like to customize this in the functional groups
created by a :mod:`~.functional_group_factory`.

Take, for example, a :class:`.CarboxylicAcid` functional group. There
are two likely ways you would like to modify
this group, ``C(=O)O``, during construction. In the first way, you want
to add a bond to the carbon atom, and delete the ``OH`` group, which is
treated as a leaving group. This is what
:class:`.CarboxylicAcidFactory` will do by default

.. code-block:: python

    import stk

    bb = stk.BuildingBlock(
        smiles='O=C(O)CCC(=O)O',
        functional_groups=[stk.CarboxylicAcidFactory()],
    )

Here, ``bb`` will have two :class:`.CarboxylicAcid` functional groups.
In each, the deleter atoms will be the oxygen and hydrogen atom of
the ``OH`` group, and the bonder atom will be the carbon atom.

Now, the second way you might want to modify a carobxylic acid group,
is to only delete the hydrogen atom of the ``OH`` group during
construction, and add a bond to the oxygen atom of the
``OH`` group. This means the hydrogen atom is the deleter atom and
the oxygen atom is the bonder atom. You can tell the
:class:`.CarboxylicAcidFactory` to create :class:`.CarboxylicAcid`
instances of this kind

.. code-block:: python

    bb2 = stk.BuildingBlock(
        smiles='O=C(O)CCC(=O)O',
        functional_groups=[
            stk.CarboxylicAcidFactory(
                # Atom number 3 corresponds to the OH oxygen atom in a
                # carboxylic acid group. THIS IS NOT THE ATOM'S ID IN
                # THE MOLECULE.
                bonders=(3, ),
                # Atom number 4 corresponds to the hydrogen atom in a
                # carboxylic acid group. THIS IS NOT THE ATOM'S ID IN
                # THE MOLECULE.
                deleters=(4, ),
            ),
        ],
    )

Here, ``bb2`` will also have two :class:`.CarboxylicAcid` functional
groups. In each, the deleter atom will be the hydrogen of the
``OH`` group and the bonder atom will be the oxygen atom of the
``OH`` group.

You might be wondering: "How do I know which number to use for
which atom in the functional group, so that I can specify the correct
atoms to be the bonders or deleters?" The docstring of
:class:`.CarboxylicAcidFactory` will tell you which number corresponds
to which atom in the functional group. The same is true for any
other :mod:`~.functional_group_factory`. Note that the number you
provide to the factory, is not the id of the atom found in the
molecule!!

Using RDKit to Optimize Molecular Structures
============================================

Molecules used by :mod:`stk` can be structure optimized both before and
after construction. One easy way to do is, is with the
:mod:`rdkit` library. You can optimize any :mod:`stk`
:class:`.Molecule`, such as a :class:`.BuildingBlock`

.. code-block:: python

    import stk
    import rdkit.Chem.AllChem as rdkit

    bb = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

    # Optimize with the MMFF force field.

    rdkit_bb = bb.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_bb)
    rdkit.MMFFOptimizeMolecule(rdkit_bb)

    # stk molecules are immutable. with_position_matrix returns a
    # a clone, holding the new position matrix.
    bb = bb.with_position_matrix(
        position_matrix=rdkit_bb.GetConformer().GetPositions(),
    )

or a :class:`.ConstructedMolecule`

.. code-block:: python

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 15),
    )

    # Optimize with the MMFF force field.

    rdkit_polymer = polymer.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_polymer)
    rdkit.MMFFOptimizeMolecule(rdkit_polymer)

    # stk molecules are immutable. with_position_matrix returns a
    # a clone, holding the new position matrix.
    polymer = polymer.with_position_matrix(
        position_matrix=rdkit_polymer.GetConformer().GetPositions(),
    )

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


Writing Molecular Files
=======================

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


Placing and Retrieving Molecules From a Database
================================================

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

Placing and Retrieving Molecular Property Values From a Database
================================================================

Extending stk
=============

There are a lot of way to extend :mod:`stk`, for example by adding
new functional groups, topology graphs, mutation operations and so on.
However, because every part of :mod:`stk` is built around an abstract
base class, all you need to do is find the appropriate abstract base
class, and create a new subclass for it. The abstract base class
will provide documentation and examples on how to create a subclass.
You can easily find the abstract base classes by looking at the
sidebar.
