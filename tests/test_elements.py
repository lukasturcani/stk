import stk


def test_initialization():
    for i in range(1, 119):
        a1 = stk.Atom(
            id=1,
            atomic_number=i,
            attr1=12,
            attr2='hi',
            _private_attr=123
        )
        a2 = stk.Atom._elements[i](2, -3, attr3='bye', attr4=321.12)

        assert a1.__class__ is a2.__class__

        assert a1.id == 1
        assert a1.charge == 0
        assert a1.attr1 == 12
        assert a1.attr2 == 'hi'
        assert a1._private_attr == 123
        assert not hasattr(a1, 'attr3')
        assert not hasattr(a1, 'attr4')

        assert a2.id == 2
        assert a2.charge == -3
        assert a2.attr3 == 'bye'
        assert a2.attr4 == 321.12
        assert not hasattr(a2, 'attr1')
        assert not hasattr(a2, 'attr2')


def test_clone(carbon):

    carbon_clone = carbon.clone()

    assert carbon_clone is not carbon
    assert carbon_clone.__class__ is carbon.__class__
    assert carbon_clone.id == carbon.id
    assert carbon_clone.charge == carbon.charge
    assert carbon_clone.attr2 == carbon.attr2
    assert carbon_clone.attr12 == carbon.attr12
    assert not hasattr(carbon_clone, '_pattr')


def test_str():
    for i in range(1, 119):
        atom = stk.Atom(
            id=1,
            atomic_number=i,
            charge=32,
            attr1=12,
            attr2='hi'
        )
        assert f'{atom}' == f'{atom.__class__.__name__}(1)'


def test_repr():
    for i in range(1, 119):
        atom = stk.Atom(
            id=1,
            atomic_number=i,
            attr1=12,
            attr2='hi',
            _private_attr=222
        )
        atom_repr = (
            f'{atom.__class__.__name__}('
            f"{atom.id}, attr1=12, attr2='hi'"
            ')'
        )
        assert f'{atom!r}' == atom_repr

        charged = stk.Atom(
            id=1,
            atomic_number=i,
            charge=-12,
            attr1=12,
            attr2='hi',
            _private_attr=222
        )
        charged_repr = (
            f'{atom.__class__.__name__}('
            f"{atom.id}, charge=-12, attr1=12, attr2='hi'"
            ')'
        )
        f'{atom.__class__.__name__}({atom.id})'
        assert f'{charged!r}' == charged_repr
