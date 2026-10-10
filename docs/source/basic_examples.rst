==============
Basic Examples
==============


.. note::

    The molecules displayed are interactive renderings.


Creating Building Blocks
========================

There are a number of ways to create :class:`.BuildingBlock` molecules,
but the two most common ways are direct creation through SMILES
strings and loading them from molecule structure files.


.. testsetup:: creating-building-blocks

    import stk
    import os
    import pathlib

    path = pathlib.Path('path/to/file.mol')
    os.makedirs(path.parent, exist_ok=True)
    stk.MolWriter().write(
        molecule=stk.BuildingBlock('NCCCN'),
        path=str(path),
    )


.. testcode:: creating-building-blocks

    import stk

    bb1 = stk.BuildingBlock('NCCCN')
    bb2 = stk.BuildingBlock.init_from_file('path/to/file.mol')

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCCN')

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb1.get_atoms(),
                bb1.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb1.get_bonds()
        ),
    )

.. testcode:: creating-building-blocks
    :hide:

    assert stk.Smiles().get_key(bb1) == 'NCCCN'
    assert stk.Smiles().get_key(bb2) == 'NCCCN'

.. testcleanup:: creating-building-blocks

    import shutil

    shutil.rmtree('path')

Look at the documentation of :class:`.BuildingBlock` to see other
available initialization methods. They are easy to find, as they all
begin with :mod:`init`.

Specifying Multiple Functional Groups
=====================================

When you create a :class:`.BuildingBlock`, you also need to specify
which atoms are modified during construction of a
:class:`.ConstructedMolecule`. This is achieved by providing the
:class:`.BuildingBlock` with
:mod:`~.functional_groups.functional_group` instances. To save you
the pain of creating function groups one by one, you can use a
:mod:`~.functional_group_factory`. If you have a building block
with bromo groups, and you want the bromo groups to be modified
during construction, you would use a :class:`.BromoFactory`

.. testcode:: specifying-multiple-functional-groups

    import stk

    bb = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

.. testcode:: specifying-multiple-functional-groups
    :hide:

    assert all(
        isinstance(fg, stk.Bromo) for fg in bb.get_functional_groups()
    )
    assert bb.get_num_functional_groups() == 2

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb.get_atoms(),
                bb.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb.get_bonds()
        ),
    )

The ``bb``, in the example above, would have two :class:`.Bromo`
functional groups. Similarly, if you have a building block with
aldehyde groups

.. testcode:: specifying-multiple-functional-groups

    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])

.. testcode:: specifying-multiple-functional-groups
    :hide:

    assert all(
        isinstance(fg, stk.Aldehyde)
        for fg in bb2.get_functional_groups()
    )
    assert bb2.get_num_functional_groups() == 2

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb.get_atoms(),
                bb.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb.get_bonds()
        ),
    )

In this example, ``bb2`` will have two :class:`.Aldehyde` functional
groups. Finally, if you have both aldehyde and bromo groups on a
molecule, and you want both to be modified during construction,
you would use both of the factories

.. testcode:: specifying-multiple-functional-groups

    bb3 = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFactory(), stk.BromoFactory()],
    )

.. testcode:: specifying-multiple-functional-groups
    :hide:

    assert (
        set(map(type, bb3.get_functional_groups()))
        == {stk.Aldehyde, stk.Bromo}
    )

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFactory(), stk.BromoFactory()],
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb.get_atoms(),
                bb.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb.get_bonds()
        ),
    )

In the example above, ``bb3`` has one :class:`.Bromo` and one
:class:`.Aldehyde` functional group.

Constructing Molecules
======================

To construct molecules, you need to create a new
:class:`.ConstructedMolecule`. The required input consists of
a :class:`.TopologyGraph`, which, in turn,  requires
:class:`.BuildingBlock` instances.

.. testcode:: constructing-molecules

    import stk

    # React the amine functional groups during construction.
    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    # React the aldehyde functional groups during construction.
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    # Build a polymer.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=4,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=4,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                polymer.get_atoms(),
                polymer.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in polymer.get_bonds()
        ),
    )

.. testcode:: constructing-molecules

    # Build a longer polymer.
    longer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=8,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=8,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                polymer.get_atoms(),
                polymer.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in polymer.get_bonds()
        ),
    )


.. testcode:: constructing-molecules
    :hide:

    assert polymer.get_num_building_block(bb1) == 4
    assert polymer.get_num_building_block(bb2) == 4
    assert longer.get_num_building_block(bb1) == 8
    assert longer.get_num_building_block(bb2) == 8


Each topology graph requires different input parameters.
For example, organic cage topology graphs only require the
:class:`.BuildingBlock` instances.

.. testcode:: constructing-molecules

    # The cage requires a building block with 3 functional groups.
    cage_bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)CC=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, cage_bb2)),
    )


.. testcode:: constructing-molecules
    :hide:

    assert cage.get_num_building_block(bb1) == 6
    assert cage.get_num_building_block(cage_bb2) == 4

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)CC=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                cage.get_atoms(),
                cage.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in cage.get_bonds()
        ),
    )

Note that this structure is far from ideal, the next example shows
how to make one with shorter bond lengths!


Read the documentation for each kind of :class:`.TopologyGraph`, for
more examples on how to initialize it, and to see what optional
parameters you have available.

Using Built-in Optimizers During Construction
=============================================

All :class:`.TopologyGraph` instances take an `optimizer` argument,
which provides efficient optimization of :mod:`stk` structures from
their `expanded` form. No optimization will be performed with the
:class:`.NullOptimizer`.

:class:`.Collapser` performs rigid translations of the building blocks
toward the centroid of the :class:`.ConstructedMolecule` until steric
clashes occur.

.. testcode:: using-built-in-optimizers-during-construction

    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=3,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

.. testcode:: using-built-in-optimizers-during-construction
    :hide:

    assert polymer.get_num_building_block(bb1) == 3
    assert polymer.get_num_building_block(bb2) == 3

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=3,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                polymer.get_atoms(),
                polymer.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in polymer.get_bonds()
        ),
    )


Similarly, :class:`.MCHammer` performs rigid translations of the
building blocks either toward the centroid of the
:class:`.ConstructedMolecule` or along the bonds formed during
construction following a Metropolis Monte Carlo algorithm with
simplified potential energy terms for the long bonds and nonbonded
interactions.

.. testcode:: using-built-in-optimizers-during-construction

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)CC=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(bb1, bb2),
            optimizer=stk.MCHammer(),
        ),
    )

.. testcode:: using-built-in-optimizers-during-construction
    :hide:

    assert cage.get_num_building_block(bb1) == 6
    assert cage.get_num_building_block(bb2) == 4

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)CC=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(bb1, bb2),
            optimizer=stk.MCHammer(),
        ),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                cage.get_atoms(),
                cage.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in cage.get_bonds()
        ),
    )

.. seealso::

    The :class:`.Collapser` and :class:`.MCHammer` optimizers use the
    algorithms from https://github.com/andrewtarzia/MCHammer.
    :mod:`stk` returns the final molecule only but further
    visualisation of the full trajectory and properties can be
    performed using the :mod:`MCHammer` code explicitly. This is useful
    for determining optimal optimization parameters, for which safe
    options are provided by default in :mod:`stk`.

Using RDKit to Optimize Molecular Structures
============================================

Molecules used by :mod:`stk` can be structure optimized both before and
after construction. One easy way to do is, is with the
:mod:`rdkit` library. You can optimize any :mod:`stk`
:class:`.Molecule`, such as a :class:`.BuildingBlock`

.. testcode:: using-rdkit-to-optimize-molecular-structures

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

.. testcode:: using-rdkit-to-optimize-molecular-structures
   :hide:

   import numpy as np

   assert np.all(np.equal(
        bb.get_position_matrix(),
        rdkit_bb.GetConformer().GetPositions(),
    ))


.. moldoc::

    import moldoc.molecule as molecule
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

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                bb.get_atoms(),
                bb.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in bb.get_bonds()
        ),
    )

or a :class:`.ConstructedMolecule`

.. testcode:: using-rdkit-to-optimize-molecular-structures

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 8),
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

.. testcode:: using-rdkit-to-optimize-molecular-structures
    :hide:

    assert np.all(np.equal(
        polymer.get_position_matrix(),
        rdkit_polymer.GetConformer().GetPositions(),
    ))


.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import rdkit.Chem.AllChem as rdkit

    bb = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 8),
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

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                polymer.get_atoms(),
                polymer.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in polymer.get_bonds()
        ),
    )

Writing Molecular Files
=======================

The simplest way to save molecules is to write them to a file.
This works with any :class:`.Molecule`, including both the
:class:`.BuildingBlock`

.. testcode:: writing-molecular-files

    import stk

    bb = stk.BuildingBlock(
        smiles='ICCBr',
        functional_groups=[stk.BromoFactory(), stk.IodoFactory()],
    )
    stk.MolWriter().write(bb, 'bb.mol')

.. testcode:: writing-molecular-files
    :hide:

    _loaded_bb = stk.BuildingBlock.init_from_file('bb.mol')
    assert stk.Smiles().get_key(_loaded_bb) == 'BrCCI'


and the :class:`.ConstructedMolecule`

.. testcode:: writing-molecular-files

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 10),
    )
    stk.MolWriter().write(polymer, 'polymer.mol')


.. testcode:: writing-molecular-files
    :hide:

    _loaded_polymer = stk.BuildingBlock.init_from_file('polymer.mol')
    assert (
        stk.Smiles().get_key(_loaded_polymer) == 'Br' + 'CC' * 10 + 'I'
    )

.. testcleanup:: writing-molecular-files

    import os

    os.remove('bb.mol')
    os.remove('polymer.mol')

You can see what file formats are supported by reading the
documentation for :meth:`~.Molecule.write`.


Visualising Molecular Datasets using chemiscope
===============================================


Requirements
------------

:mod:`chemiscope` makes it easy for you to write a ``.json`` or ``.json.gz``
containing :mod:`stk` molecules and their properties; see an example__.

__ https://chemiscope.org/docs/examples/9-showing_custom_bonds.html

To get :mod:`.chemiscope`, you can install it with pip::

  $ pip install chemiscope

Writing to file
---------------

:mod:`chemiscope` can be used through the `web-interface`__ or can be embedded
into read-the-docs pages through a sphinx package, like for `toy models`__.

__ https://chemiscope.org/

__ https://cgmodels.readthedocs.io/en/latest/cg_model_jul2023_2p3_ton.html

Either way, you need to write a ``.json`` or ``.json.gz`` file

.. testcode:: saving-to-chemiscope

    import stk
    import chemiscope

    structures = [
        stk.BuildingBlock(smiles="NCCN"),
        stk.BuildingBlock(smiles="NCCCN"),
        # A mostly optimised cage molecule.
        stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles="NCCN",
                        functional_groups=[stk.PrimaryAminoFactory()],
                    ),
                    stk.BuildingBlock(
                        smiles="O=CC(C=O)C=O",
                        functional_groups=[stk.AldehydeFactory()],
                    ),
                ),
                optimizer=stk.MCHammer(),
            ),
        )
    ]

    # Write their properties to a dictionary.
    properties = {
        "num_atoms": [molecule.get_num_atoms() for molecule in structures],
        "num_bonds": [
            len(list(molecule.get_bonds())) for molecule in structures
        ],
    }

    # Define stk bonding.
    shape_dict = chemiscope.convert_stk_bonds_as_shapes(
        frames=structures,
        bond_color="#fc5500",
        bond_radius=0.12,
    )

    # Write the shape string for settings to turn them on automatically.
    shape_string = ",".join(shape_dict.keys())

    # Write to file.
    chemiscope.write_input(
        path="stk_example.json.gz",
        frames=structures,
        properties=properties,
        meta=dict(name="A name."),
        settings=chemiscope.quick_settings(
            x="num_atoms",
            y="num_bonds",
            color="",
            structure_settings={
                "shape": shape_string,
                "atoms": True,
                "bonds": False,
                "spaceFilling": False,
            },
        ),
        shapes=shape_dict,
    )

.. testcleanup:: saving-to-chemiscope

    import pathlib

    pathlib.Path("stk_example.json.gz").unlink()

.. seealso::

    chemiscope: https://chemiscope.org/


Specifying Functional Groups Individually
=========================================

If you want to be more precise about which functional groups get
created, you can provide them directly to the :class:`.BuildingBlock`.
For example, if you have multiple bromo groups on a molecule, but
you only want to use one during construction

.. testcode:: specifying-functional-groups-individually

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

.. testcode:: specifying-functional-groups-individually
    :hide:

    _fg, = bb.get_functional_groups()
    assert type(_fg) is stk.Bromo
    assert set(_fg.get_atom_ids()) == {0, 1}

When creating a :class:`.Bromo` functional group, you have to
specify things like which atoms have bonds added during construction,
and which ones are removed during construction. These are specified by
the `bonders` and `deleters` parameters, respectively. You can add
as many functional groups to :class:`.BuildingBlock` as you like
in this way, and you can mix different types of
:mod:`~.functional_groups.functional_group`. You can even mix
a :mod:`~.functional_groups.functional_group` instances with
:mod:`~.functional_group_factory` instances.

Changing Bonder, Placer and Deleter Atoms in Functional Group Factories
=======================================================================

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

.. testcode:: changing-bonder-and-deleter-atoms-in-fg-factories

    import stk

    bb = stk.BuildingBlock(
        smiles='O=C(O)CCC(=O)O',
        functional_groups=[stk.CarboxylicAcidFactory()],
    )

.. testcode:: changing-bonder-and-deleter-atoms-in-fg-factories
    :hide:

    for _fg in bb.get_functional_groups():
        assert isinstance(_fg, stk.CarboxylicAcid)

        _bonder, = _fg.get_bonders()
        assert isinstance(_bonder, stk.C)

        _deleter1, _deleter2 = _fg.get_deleters()
        assert type(_deleter1) is not type(_deleter2)
        assert {type(_deleter1), type(_deleter2)} == {stk.O, stk.H}

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

.. testcode:: changing-bonder-and-deleter-atoms-in-fg-factories

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

.. testcode:: changing-bonder-and-deleter-atoms-in-fg-factories
    :hide:

    for _fg in bb2.get_functional_groups():
        assert isinstance(_fg, stk.CarboxylicAcid)

        _bonder, = _fg.get_bonders()
        assert isinstance(_bonder, stk.O)

        _deleter, = _fg.get_deleters()
        assert isinstance(_deleter, stk.H)

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

Finally, each :mod:`~.functional_group_factory` allows for the user to
provide the placer ids, which tell :mod:`stk` which atoms in a functional
group should be used for aligning and placing the molecule during
construction. In most cases, this is not an important parameter and the
default behaviour (where placers and bonders are equivalent) is
appropropriate. However, in cases where your building block is small,
you may get `NAN` or `division by zero` errors due to the vectors in a
building block overlapping. If this occurs, we recommend using bonder
and deleter atoms as placers like so:

.. testcode:: changing-bonder-and-deleter-atoms-in-fg-factories

    bb2 = stk.BuildingBlock(
        smiles='C(Br)(Br)Br',
        functional_groups=stk.BromoFactory(
            # Defaults.
            bonders=(0,),
            deleters=(1,),
            # Include both atoms.
            placers=(0, 1)
        ),
    )

Handling Molecules with Metal Atoms and Dative Bonds
====================================================

All :mod:`stk` :class:`.Molecule` instances (such as
:class:`.BuildingBlock` and :class:`.ConstructedMolecule`) can contain
metal atoms and handle various coordination reactions.
In order to represent dative bonds in these systems, a bond order of
9 is used.

Furthermore, when working with metal-containing systems, any
:class:`.BuildingBlock` initialization functions that require ETKDG
may fail, because the ETKDG algorithm is liable to fail in these cases.
In cases like this, you probably want to set the position matrix
explicitly, which will mean that ETKDG will not be used.

.. testcode:: handling-molecules-with-metal-atoms-and-dative-bonds

    import stk

    bb = stk.BuildingBlock('[Fe+2]', position_matrix=[[0., 0., 0.]])

.. testcode:: handling-molecules-with-metal-atoms-and-dative-bonds
    :hide:

    import numpy as np

    assert stk.Smiles().get_key(bb) == '[Fe+2]'
    assert np.all(np.equal(
        bb.get_position_matrix(),
        [[0., 0., 0.]]
    ))


If you want to get a more complex position matrix, defining a
function may be a good idea

.. testcode:: handling-molecules-with-metal-atoms-and-dative-bonds

    import rdkit.Chem.AllChem as rdkit


    def get_position_matrix(smiles):
        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        rdkit.EmbedMolecule(molecule, randomSeed=12)
        rdkit.UFFOptimizeMolecule(molecule)
        return molecule.GetConformer().GetPositions()


    smiles = 'CCC[OH]->[Fe+2]'
    bb2 = stk.BuildingBlock(
        smiles=smiles,
        position_matrix=get_position_matrix(smiles),
    )

.. testcode:: handling-molecules-with-metal-atoms-and-dative-bonds
    :hide:

    assert stk.Smiles().get_key(bb2) == smiles
    assert np.all(np.equal(
        bb2.get_position_matrix(),
        get_position_matrix(smiles),
    ))


Finally, :mod:`stk` will also read bonds from ``.mol`` files,
which have a bond order of 9, as dative.

Making Keys for Molecules with Dative Bonds
===========================================

Dative bonds are not defined in an InChI or InChiIKey.
Therefore, when storing metal-containing molecules in a
database, a different key is required. Because dative bonds are
implemented in SMILES, the SMILES string makes a
useful key for metal-containing molecules. You can use the
:class:`.Smiles` key maker for this purpose

.. testsetup:: making-keys-for-molecules-with-dative-bonds

    import stk

    # Change the default database used, so that when a developer runs
    # the doctests locally, their "stk" database is not contaminated.
    _test_database = '_stk_doctest_database'
    _old_molecule_init = stk.MoleculeMongoDb
    stk.MoleculeMongoDb = lambda mongo_client, jsonizer: (
        _old_molecule_init(
            mongo_client=mongo_client,
            jsonizer=jsonizer,
            database=_test_database,
        )
    )

    # Change the database MongoClient will connect to.

    import os
    import pymongo

    _mongo_client = pymongo.MongoClient
    _mongodb_uri = os.environ.get(
        'MONGODB_URI',
        'mongodb://localhost:27017/'
    )
    pymongo.MongoClient = lambda: _mongo_client(_mongodb_uri)

.. testcode:: making-keys-for-molecules-with-dative-bonds

    import stk
    import pymongo

    db = stk.MoleculeMongoDb(
        mongo_client=pymongo.MongoClient(),
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(stk.Smiles(), ),
        ),
    )
    bb = stk.BuildingBlock('BrO->[Fe+2]')
    db.put(bb)
    # Use the Smiles() key maker to get the retrieval SMILES,
    # to make sure it has canonical atom ordering.
    canonical_smiles = stk.Smiles().get_key(bb)
    retrieved_bb = db.get({'SMILES': canonical_smiles})

.. testcode:: making-keys-for-molecules-with-dative-bonds
    :hide:

    _smiles = stk.Smiles()
    assert _smiles.get_key(bb) == _smiles.get_key(retrieved_bb)

.. testcleanup:: making-keys-for-molecules-with-dative-bonds

    pymongo.MongoClient().drop_database(_test_database)
    stk.MoleculeMongoDb = _old_molecule_init
    pymongo.MongoClient = _mongo_client

Creating New Topology Graphs with Existing Vertices
===================================================

The vertex classes that make up topology graphs in :mod:`.stk` can be
accessed to speed up the implemention of new and arbitrary topology
graphs (as shown below). The exact details of how vertices can be used
to implement new topology graphs depends on the topology graph, so read
that documentation for further examples. For metal complexes, you would
read the documentation of :class:`.MetalComplex`.

.. testcode:: creating-new-topology-graphs-with-existing-vertices

    import stk
    import numpy as np

    class NewMetalComplexTopology(stk.metal_complex.MetalComplex):

        _metal_vertex_prototypes = (
            stk.metal_complex.MetalVertex(0, (1., 0., 0.)),
        )
        _ligand_vertex_prototypes = (
            stk.metal_complex.MonoDentateLigandVertex(1, (2., 0., 0.)),
            stk.metal_complex.MonoDentateLigandVertex(2, (0., 0., 0.)),
        )

        # Define Edges below.
        _edge_prototypes = (
            stk.Edge(
                id=0,
                vertex1=_metal_vertex_prototypes[0],
                vertex2=_ligand_vertex_prototypes[0],
            ),
            stk.Edge(
                id=1,
                vertex1=_metal_vertex_prototypes[0],
                vertex2=_ligand_vertex_prototypes[1],
            ),
        )

    # Build new metal complex.
    ligand = stk.BuildingBlock('NCC', [stk.PrimaryAminoFactory()])
    metal = stk.BuildingBlock(
        smiles='[Fe+2]',
        functional_groups=(
            stk.SingleAtom(stk.Fe(0, charge=2))
            for i in range(2)
        ),
        position_matrix=np.array([[0., 0., 0.]]),
    )
    complex = stk.ConstructedMolecule(
        topology_graph=NewMetalComplexTopology(
            metals=metal,
            ligands=ligand,
        )
    )

.. testcode:: creating-new-topology-graphs-with-existing-vertices
    :hide:

    assert complex.get_num_atoms() == (
        (ligand.get_num_atoms()-2)*2 + 1
    )


Extending stk
=============

There are a lot of ways to extend :mod:`stk`, for example by adding
new functional groups, topology graphs, mutation operations and so on.
However, because every part of :mod:`stk` is built around an abstract
base class, all you need to do is find the appropriate abstract base
class, and create a new subclass for it. The abstract base class
will provide documentation and examples on how to create a subclass.
You can easily find the abstract base classes by looking at the
sidebar.
