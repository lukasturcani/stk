import stk
import itertools as it


def is_fg_match(fg1, fg2):
    return (
        fg1.atoms == fg2.atoms
        and fg1.bonders == fg2.bonders
        and fg1.deleters == fg2.deleters
    )


def test_get_functional_groups(amine2, aldehyde3):
    amine = stk.fg_types['amine']
    aldehyde = stk.fg_types['aldehyde']

    for i, fg in enumerate(amine.get_functional_groups(amine2)):
        assert any(
            is_fg_match(fg, other) for other in amine2.func_groups
        )
        assert len(fg.atoms) == 3
        assert len(fg.bonders) == 1
        assert fg.bonders[0].__class__ is stk.N
        assert len(fg.deleters) == 2
        assert all(a.__class__ is stk.H for a in fg.deleters)
    assert i == 1
    assert not list(aldehyde.get_functional_groups(amine2))

    for i, fg in enumerate(aldehyde.get_functional_groups(aldehyde3)):
        assert any(
            is_fg_match(fg, other) for other in aldehyde3.func_groups
        )
        assert len(fg.atoms) == 3
        assert len(fg.bonders) == 1
        assert fg.bonders[0].__class__ is stk.C
        assert len(fg.deleters) == 1
        assert fg.deleters[0].__class__ is stk.O
    assert i == 2
    assert not list(amine.get_functional_groups(aldehyde3))


def test_clone(aldehyde3, hydrogen, carbon):
    fg = aldehyde3.func_groups[0]
    fg_a0, fg_a1 = fg.atoms[:2]
    clone = fg.clone({
        fg_a0: hydrogen,
        fg_a1: carbon
    })

    for a0, a1 in it.zip_longest(fg.atoms, clone.atoms):
        if a0 is fg_a0:
            assert a1 is hydrogen
        elif a0 is fg_a1:
            assert a1 is carbon
        else:
            assert a0 is a1

    for a0, a1 in it.zip_longest(fg.bonders, clone.bonders):
        if a0 is fg_a0:
            assert a1 is hydrogen
        elif a0 is fg_a1:
            assert a1 is carbon
        else:
            assert a0 is a1

    for a0, a1 in it.zip_longest(fg.deleters, clone.deleters):
        if a0 is fg_a0:
            assert a1 is hydrogen
        elif a0 is fg_a1:
            assert a1 is carbon
        else:
            assert a0 is a1


def test_get_atom_ids(aldehyde3):
    for fg in aldehyde3.func_groups:
        atoms = it.zip_longest(fg.get_atom_ids(), fg.atoms)
        for atom, atom_id in atoms:
            assert atom.id == atom_id


def test_get_bonder_ids(aldehyde3):
    for fg in aldehyde3.func_groups:
        bonders = it.zip_longest(fg.get_bonder_ids(), fg.bonders)
        for bonder, bonder_id in bonders:
            assert bonder.id == bonder_id


def test_get_deleter_ids(aldehyde3):
    for fg in aldehyde3.func_groups:
        deleters = it.zip_longest(fg.get_deleter_ids(), fg.deleters)
        for deleter, deleter_id in deleters:
            assert deleter.id == deleter_id
