import itertools as it


def test_clone(functional_group):
    clone = functional_group.clone()
    atoms = it.zip_longest(
        functional_group.atoms,
        clone.atoms,
    )
    for a1, a2 in atoms:
        assert a1 is not a2
        assert a1.id == a2.id

    bonders = it.zip_longest(
        functional_group.bonders,
        clone.bonders,
    )
    for b1, b2 in bonders:
        assert b1 is not b2
        assert b1.id == b2.id

    deleters = it.zip_longest(
        functional_group.deleters,
        clone.deleters,
    )
    for d1, d2 in deleters:
        assert d1 is not d2
        assert d1.id == d2.id
