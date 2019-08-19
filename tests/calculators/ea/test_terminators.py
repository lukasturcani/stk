import stk


def test_any_terminator():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(20)))

    terminator = stk.AnyTerminator(
        stk.NumGenerations(7),
        stk.NumGenerations(15)
    )
    assert not terminator.terminate(pop1)
    assert terminator.terminate(pop2)
    assert terminator.terminate(pop3)


def test_all_terminators():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(20)))

    terminator = stk.AllTerminators(
        stk.NumGenerations(7),
        stk.NumGenerations(15)
    )
    assert not terminator.terminate(pop1)
    assert not terminator.terminate(pop2)
    assert terminator.terminate(pop3)


def test_num_generations():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(15)))

    terminator = stk.NumGenerations(10)

    assert not terminator.terminate(pop1)
    assert terminator.terminate(pop2)
    assert terminator.terminate(pop3)


def test_molecule_present(amine2):
    pop1 = stk.Population(stk.Population())
    pop2 = stk.Population(stk.Population(amine2))

    terminator = stk.MoleculePresent(amine2)

    assert not terminator.terminate(pop1)
    assert terminator.terminate(pop2)


def test_fitness_plateau():
    bbs = [
        stk.BuildingBlock.__new__(stk.BuildingBlock)
        for i in range(5*10)
    ]
    for i, bb in enumerate(bbs):
        bb._key = i

    pop = stk.Population(
        *(stk.Population(*bbs[i:i+5]) for i in range(0, len(bbs), 5))
    )
    for i, mol in enumerate(pop):
        mol.fitness = i

    terminator = stk.FitnessPlateau(2)

    assert not terminator.terminate(pop)

    for mol in pop:
        mol.fitness = 2
    assert not terminator.terminate(pop)

    pop2 = stk.Population(*pop)
    pop3 = stk.Population(pop2, pop2)
    assert terminator.terminate(pop3)
