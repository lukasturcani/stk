
def is_equivalent_atom(atom1, atom2):
    return (
        atom1.id == atom2.id
        and atom1.charge == atom2.charge
        and atom1.__class__ is atom2.__class__
    )


def is_equivalent_bond(bond1, bond2):
    return (
        bond1.__class__ is bond2.__class__
        and bond1.order == bond2.order
        and is_equivalent_atom(bond1.atom1, bond2.atom1)
        and is_equivalent_atom(bond1.atom2, bond2.atom2)
        and bond1.periodicity == bond2.periodicity
    )


def is_equivalent_fg(fg1, fg2):
    equivalent_atoms = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_atom_ids(), fg2.get_atom_ids())
    )
    equivalent_bonders = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_bonder_ids(), fg2.get_bonder_ids())
    )
    equivalent_deleters = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_deleter_ids(), fg2.get_deleter_ids())
    )
    return (
        fg1.__class__ is fg2.__class__
        and equivalent_atoms
        and equivalent_bonders
        and equivalent_deleters
    )


def is_equivalent_building_block(building_block1, building_block2):
    atoms = it.zip_longest(
        building_block1.atoms,
        building_block2.atoms,
    )
    equivalent_atoms = all(
        is_equivalent_atom(a1, a2) for a1, a2 in atoms
    )

    bonds = it.zip_longest(
        building_block1.bonds,
        building_block2.bonds,
    )
    equivalent_bonds = all(
        is_equivalent_bond(b1, b2) for b1, b2 in bonds
    )

    fgs = it.zip_longest(
        building_block1.func_groups,
        building_block2.func_groups,
    )
    equivalent_fgs = all(
        is_equivalent_fg(fg1, fg2) for fg1, fg2 in fgs
    )
    return equivalent_atoms and equivalent_bonds and equivalent_fgs


class TestInitFromFile:
    @pytest.fixture(
        params=[
            'building_block.mol',
            'building_block.pdb',
        ],
    )
    def filename(self, request):
        return request.param

    def test(self, tmpdir, filename, building_block):
        path = str(tmpdir / filename)
        building_block.write(path)

        loaded = stk.BuildingBlock.init_from_file(
            path=path,
            functional_groups={
                fg.fg_type.name for fg in building_block.func_groups
            },
        )

        atoms = it.zip_longest(building_block.atoms, loaded.atoms)
        for a1, a2 in atoms:
            assert is_equivalent_atom(a1, a2)

        bonds = it.zip_longest(building_block.bonds, loaded.bonds)
        for b1, b2 in bonds:
            assert is_equivalent_bond(b1, b2)

        fgs = it.zip_longest(
            building_block.func_groups,
            loaded.func_groups
        )
        for fg1, fg2 in fgs:
            assert is_equivalent_fg(fg1, fg2)
