import stk


def func_groups(building_blocks):
    fgs = (fg for bb in building_blocks for fg in bb.func_groups)
    yield from stk.dedupe(fgs, key=lambda fg: fg.fg_type.name)


def _test_reaction(
    reactor,
    atom_change_per_reaction,
    bond_change_per_reaction,
    costruction_bonds_per_reaction,
    expected_construction_bond_order
):
    mol = reactor._mol
    num_start_atoms = len(mol.atoms)
    num_start_bonds = len(mol.bonds)
    edge_clones = reactor._mol._edge_clones
    original_fgs = tuple(mol.func_groups)
    assert len(original_fgs) == 12
    reacted_fgs = []
    deleters = set()

    assert len(mol.construction_bonds) == 0

    periodicities = [
        (0, 0, 0),
        (1, 0, -1)
    ]

    num_expected_periodic_bonds = 0
    degrees = {}
    for i, edge in enumerate(edge_clones):
        for fg in edge.get_func_groups():
            for atom in fg.deleters:
                degree = 0
                for bond in mol.bonds:
                    if bond.atom1 is atom or bond.atom2 is atom:
                        degree += 1
                degrees[atom] = degree

            deleters.update(fg.deleters)

        start_bonds = len(mol.bonds)
        reactor.add_reaction(
            func_groups=edge.get_func_groups(),
            periodicity=periodicities[i % 2]
        )
        if i % 2 == 1:
            num_expected_periodic_bonds += start_bonds - len(mol.bonds)

        reacted_fgs.extend(edge.get_func_groups())
    reactor.finalize()
    assert len(reacted_fgs) == 10

    # Make sure that unreacted functional groups are unchanged.
    num_original_atoms = {
        fg.fg_type.name: len(fg.atoms)
        for fg in func_groups(mol.building_block_vertices.keys())
    }
    num_original_bonders = {
        fg.fg_type.name: len(fg.bonders)
        for fg in func_groups(mol.building_block_vertices.keys())
    }
    num_original_deleters = {
        fg.fg_type.name: len(fg.deleters)
        for fg in func_groups(mol.building_block_vertices.keys())
    }

    num_unreacted_fgs = 0
    for fg in original_fgs:
        if fg not in reacted_fgs:
            fg_name = fg.fg_type.name
            assert len(fg.atoms) == num_original_atoms[fg_name]
            assert len(fg.bonders) == num_original_bonders[fg_name]
            assert len(fg.deleters) == num_original_deleters[fg_name]
            num_unreacted_fgs += 1
    assert num_unreacted_fgs == 2

    # Make sure the deleter atoms got purged from the functional
    # groups.
    for fg in reacted_fgs:
        assert not fg.deleters
        fg_name = fg.fg_type.name
        num_expected_atoms = (
            num_original_atoms[fg_name] -
            num_original_deleters[fg_name]
        )
        assert len(fg.atoms) == num_expected_atoms
        assert all(atom not in deleters for atom in fg.atoms)

        assert len(fg.bonders) == num_original_bonders[fg_name]
        assert all(atom not in deleters for atom in fg.bonders)

    # Make sure the deleters atoms are not present in the molecule.
    assert all(atom not in deleters for atom in mol.atoms)
    assert (
        len(mol.atoms) == num_start_atoms + atom_change_per_reaction*5
    )

    # Make sure the correct number of construction bonds was made.
    assert (
        len(mol.construction_bonds) == costruction_bonds_per_reaction*5
    )

    # Make sure all constructed bonds have the correct bond oder.
    for bond in mol.construction_bonds:
        assert bond.order == expected_construction_bond_order(bond)

    # Make sure the correct number of bonds is left.
    assert (
        len(mol.bonds) == num_start_bonds - bond_change_per_reaction*5
    )

    # Make sure construction bonds are shared with bonds.
    bonds = set(mol.bonds)
    assert all(bond in bonds for bond in mol.construction_bonds)

    # Make sure the correct amount of bonds is periodic
    num_periodic_bonds = sum(
        1 for bond in mol.bonds if bond.is_periodic()
    )
    for bond in mol.bonds:
        if bond.is_periodic():
            assert bond.periodicity == (1, 0, -1)
        else:
            assert bond.periodicity == (0, 0, 0)
    assert num_expected_periodic_bonds == num_periodic_bonds


def test_react_any_single(make_reactor, amine2):
    reactor = make_reactor(
        building_blocks=[amine2, amine2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )
    _test_reaction(
        reactor=reactor,
        atom_change_per_reaction=-4,
        bond_change_per_reaction=-3,
        costruction_bonds_per_reaction=1,
        expected_construction_bond_order=lambda bond: 1
    )


def test_react_any_double(make_reactor, amine2, aldehyde2):
    reactor = make_reactor(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )
    _test_reaction(
        reactor=reactor,
        atom_change_per_reaction=-3,
        bond_change_per_reaction=-2,
        costruction_bonds_per_reaction=1,
        expected_construction_bond_order=lambda bond: 2
    )


def test_react_diol_with_dihalogen(
    make_reactor,
    diol2,
    difluorene_dibromine
):
    reactor = make_reactor(
        building_blocks=[diol2, difluorene_dibromine],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )
    _test_reaction(
        reactor=reactor,
        atom_change_per_reaction=-4,
        bond_change_per_reaction=-2,
        costruction_bonds_per_reaction=2,
        expected_construction_bond_order=lambda bond: 1
    )


def test_react_boronic_acid_with_diol(
    make_reactor,
    boronic_acid2,
    diol2
):
    reactor = make_reactor(
        building_blocks=[boronic_acid2, diol2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )
    _test_reaction(
        reactor=reactor,
        atom_change_per_reaction=-6,
        bond_change_per_reaction=-4,
        costruction_bonds_per_reaction=2,
        expected_construction_bond_order=lambda bond: 1
    )


def test_react_ring_amine_with_ring_amine(make_reactor, ring_amine):
    reactor = make_reactor(
        building_blocks=[ring_amine, ring_amine],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )
    _test_reaction(
        reactor=reactor,
        atom_change_per_reaction=3,
        bond_change_per_reaction=6,
        costruction_bonds_per_reaction=12,
        expected_construction_bond_order=lambda bond: 1
    )
