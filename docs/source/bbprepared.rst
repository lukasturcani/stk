==========================================
Preparing building blocks for construction
==========================================

.. note::

    This requires the installation of ``bbprepared`` as ``pip install bbprepared``.

:mod:`stk` uses the ETKDG algorithm in :mod:`rdkit` to embed 3D conformers to a
SMILES string. However, it is often necessary to have specific building block
conformations to ease larger molecule construction. For this, Andrew developed
:mod:`bbprepared`.

Aligning Binding Sites
======================

For example, a ditopic building block with rotatable bonds is better suited for
construction when both binding sites point in the same direction. This can be
achieved with a few lines of code, where the alignment applies to your
prescribed functional groups.

.. code:: python

    import stk
    import bbprepared

    # stk BuildingBlock.
    bb = stk.BuildingBlock(
        smiles="c1cncc(-c2ccc(-c3cccnc3)cc2)c1",
        functional_groups=stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    )

    # Define a generator function, using rdkit's ETKDG.
    generator = bbprepared.generators.ETKDG(num_confs=30)
    # Get a bbprepared Ensemble containing the 30 conformers.
    ensemble = generator.generate_conformers(bb)

    # Find the conformer with the lowest score, where the score here is the
    # alignment of the binding vecgtors.
    process = bbprepared.DitopicFitter(ensemble=ensemble)
    min_molecule = process.get_minimum()

    # This is now an stk BuildingBlock.
    bb = min_molecule.molecule


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        bb = stk.BuildingBlock.init_from_file("source/_static/bbprep1_in.mol")
    except OSError:
        bb = stk.BuildingBlock.init_from_file("_static/bbprep1_in.mol")

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

Becomes:

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        min_molecule = stk.BuildingBlock.init_from_file(
            "source/_static/bbprep1_out.mol"
        )
    except OSError:
        min_molecule = stk.BuildingBlock.init_from_file(
            "_static/bbprep1_out.mol"
        )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                min_molecule.get_atoms(),
                min_molecule.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in min_molecule.get_bonds()
        ),
    )



Filtering Functional Groups
===========================

Functional group factories can often find undesired binding sites, so bbprepared
makes it simple to iterate through options, leading to two different cages:

.. code:: python

    import stk
    import bbprepared

    bb = stk.BuildingBlock(
        smiles="Brc1ccc(cc1)c2ccc(Br)cc2Br",
        functional_groups=stk.BromoFactory(),
    )

    # Get the two closest.
    bb1 = bbprepared.ClosestFGs().modify(
        building_block=bb,
        desired_functional_groups=2,
    )

    # And the two furthest.
    bb2 = bbprepared.FurthestFGs().modify(
        building_block=bb,
        desired_functional_groups=2,
    )

    # And use them both in construction.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=[bb1, bb2],
            repeating_unit='AB',
            num_repeating_units=1,
        ),
    )


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        bb = stk.BuildingBlock.init_from_file(
            "source/_static/bbprep2_in.mol"
        )
    except OSError:
        bb = stk.BuildingBlock.init_from_file(
            "_static/bbprep2_in.mol"
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


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        polymer = stk.BuildingBlock.init_from_file(
            "source/_static/bbprep2_out.mol"
        )
    except OSError:
        polymer = stk.BuildingBlock.init_from_file(
            "_static/bbprep2_out.mol"
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


Other Capabilities
==================

:mod:`bbprepared` also allows you to ``Planarfy`` building blocks, find the
lowest energy conformer with various methods and scan bonds, angles and torsions.

Please see the
`bbprepared documentation <https://bbprepared.readthedocs.io/en/latest/index.html>`_
and `recipes <https://bbprepared.readthedocs.io/en/latest/recipes.html>`_.
