def test_init(get_atom, id, charge):
    atom = get_atom(id, charge)
    assert atom.get_id() == id
    assert atom.get_charge() == charge
