import stk


def test_initialization():
    for i in range(1, 119):
        a1 = stk.Atom(id=1, atomic_number=i, attr1=12, attr2='hi')
        a2 = stk.Atom._elements[i](id=2, attr3='bye', attr4=321.12)

        assert a1.__class__ is a2.__class__

        assert a1.attr1 == 12
        assert a1.attr2 == 'hi'
        assert not hasattr(a1, 'attr3')
        assert not hasattr(a1, 'attr4')

        assert a2.attr3 == 'bye'
        assert a2.attr4 == 321.12
        assert not hasattr(a2, 'attr1')
        assert not hasattr(a2, 'attr2')


def test_clone():
    carbon = stk.C(
        id=109,
        attr1='something',
        attr2=123,
        attr3=None,
        attr4=13.322,
        _attr5='private'
    )

    carbon_clone = carbon.clone()

    assert carbon_clone is not carbon
    assert carbon_clone.__class__ == carbon.__class__
    assert carbon_clone.attr1 == carbon.attr1
    assert carbon_clone.attr2 == carbon.attr2
    assert carbon_clone.attr3 == carbon.attr3
    assert carbon_clone.attr4 == carbon.attr4
    assert carbon._attr5 == 'private'
    assert not hasattr(carbon_clone, '_attr5')
