

def test_initialization(
    bond,
    periodic_bond,
    hydrogen,
    carbon,
    lithium,
    chlorine
):

    assert bond.atom1 is hydrogen
    assert bond.atom2 is carbon
    assert bond.order == 2
    assert bond.attr1 == 1
    assert bond.attr2 == '2'
    assert bond._attr3 == 12.2
    assert not hasattr(bond, 'direction')
    assert not hasattr(bond, 'attr10')
    assert not hasattr(bond, 'attr20')
    assert not hasattr(bond, '_attr30')

    assert periodic_bond.atom1 is lithium
    assert periodic_bond.atom2 is chlorine
    assert periodic_bond.order == 21
    assert periodic_bond.direction == [1, 0, -1]
    assert periodic_bond.attr10 == 16
    assert periodic_bond.attr20 == '26'
    assert periodic_bond._attr30 == 126.2
    assert not hasattr(periodic_bond, 'attr1')
    assert not hasattr(periodic_bond, 'attr2')
    assert not hasattr(periodic_bond, '_attr3')


def test_clone(bond, periodic_bond):
    bond_clone = bond.clone()

    assert bond_clone is not bond
    assert bond_clone.__class__ is bond.__class__
    assert bond_clone.atom1 is bond.atom1
    assert bond_clone.atom2 is bond.atom2
    assert bond_clone.order == bond.order
    assert bond_clone.attr1 == bond.attr1
    assert bond_clone.attr2 == bond.attr2
    assert not hasattr(bond_clone, '_attr3')

    mapped_clone = bond.clone({
        bond.atom1: bond.atom2,
        bond.atom2: bond.atom1
    })

    assert mapped_clone is not bond
    assert mapped_clone is not bond_clone
    assert mapped_clone.__class__ is bond.__class__
    assert mapped_clone.atom1 is bond.atom2
    assert mapped_clone.atom2 is bond.atom1
    assert mapped_clone.order == bond.order
    assert mapped_clone.attr1 == bond.attr1
    assert mapped_clone.attr2 == bond.attr2
    assert not hasattr(mapped_clone, '_attr3')

    periodic_clone = periodic_bond.clone()
    assert periodic_clone is not periodic_bond
    assert periodic_clone.__class__ is periodic_bond.__class__
    assert periodic_clone.atom1 is periodic_bond.atom1
    assert periodic_clone.atom2 is periodic_bond.atom2
    assert periodic_clone.direction == periodic_bond.direction
    assert periodic_clone.order == periodic_bond.order
    assert periodic_clone.attr10 == periodic_bond.attr10
    assert periodic_clone.attr20 == periodic_bond.attr20
    assert not hasattr(periodic_clone, '_attr30')

    mapped_periodic_clone = periodic_bond.clone({
        periodic_bond.atom1: periodic_bond.atom2
    })
    assert mapped_periodic_clone is not periodic_bond
    assert mapped_periodic_clone is not periodic_clone
    assert mapped_periodic_clone.__class__ is periodic_bond.__class__
    assert mapped_periodic_clone.atom1 is periodic_bond.atom2
    assert mapped_periodic_clone.atom2 is periodic_bond.atom2
    assert mapped_periodic_clone.direction == periodic_bond.direction
    assert mapped_periodic_clone.order == periodic_bond.order
    assert mapped_periodic_clone.attr10 == periodic_bond.attr10
    assert mapped_periodic_clone.attr20 == periodic_bond.attr20
    assert not hasattr(mapped_periodic_clone, '_attr30')


def test_is_periodic(bond, periodic_bond):
    assert not bond.is_periodic()
    assert periodic_bond.is_periodic()


def test_str(bond, periodic_bond):
    assert f'{bond}' == 'Bond(H(1), C(333), 2)'
    assert f'{periodic_bond}' == 'PeriodicBond(Li(121), Cl(786), 21)'


def test_repr(bond, periodic_bond):
    bond_repr = (
        "Bond(H(1, charge=-3, attr1='12', attr2=32), "
        "C(333, attr2=None, attr12='abc'), "
        "2, "
        "attr1=1, attr2='2')"
    )
    assert f'{bond!r}' == bond_repr

    periodic_repr = (
        "Bond(Li(121, charge=2, alpha=232, beta='bvc'), "
        "Cl(786, a=9, b='bbb'), "
        "21, "
        "periodicity=(1, 0, -1), "
        "attr10=16, attr20='26')"
    )
    assert f'{periodic_bond!r}' == periodic_repr
