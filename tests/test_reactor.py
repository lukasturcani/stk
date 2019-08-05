import stk


def test_default_reaction(make_reactor, amine2, aldehyde2):
    reactor = make_reactor(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )

    # Carry out the reactions.
    mol = reactor._mol
    num_start_atoms = len(mol.atoms)
    edge_clones = reactor._mol._edge_clones
    topology_graph = reactor._mol.topology_graph
    reacted_fgs = []
    deleters = set()
    for fgs in topology_graph._get_bonded_fgs(mol, edge_clones):
        for fg in fgs:
            deleters.update(fg.deleters)

        reactor.add_reaction(*fgs)
        reacted_fgs.extend(fgs)
    reactor.finalize()

    # Make sure the deleter atoms got purged from the functional
    # groups.
    for fg in reacted_fgs:
        assert not fg.deleters

        for atom in fg.atoms:
            assert atom not in deleters

        for atom in fg.bonders:
            assert atom not in deleters

    # Make sure the deleters atoms are not present in the molecule.
    for atom in mol.atoms:
        assert atom not in deleters
    assert len(mol.atoms) == num_start_atoms - len(deleters)


def test_diol_with_dihalogen(reactor):
    ...
