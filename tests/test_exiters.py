import stk


def test_any_exiter():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(20)))

    exiter = stk.AnyExiter(
        stk.NumGenerations(7),
        stk.NumGenerations(15)
    )
    assert not exiter.exit(pop1)
    assert exiter.exit(pop2)
    assert exiter.exit(pop3)


def test_all_exiters():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(20)))

    exiter = stk.AllExiters(
        stk.NumGenerations(7),
        stk.NumGenerations(15)
    )
    assert not exiter.exit(pop1)
    assert not exiter.exit(pop2)
    assert exiter.exit(pop3)


def test_num_generations():
    pop1 = stk.Population(*(stk.Population() for i in range(5)))
    pop2 = stk.Population(*(stk.Population() for i in range(10)))
    pop3 = stk.Population(*(stk.Population() for i in range(15)))

    exiter = stk.NumGenerations(10)

    assert not exiter.exit(pop1)
    assert exiter.exit(pop2)
    assert exiter.exit(pop3)


def test_molecule_present(amine2):
    pop1 = stk.Population(stk.Population())
    pop2 = stk.Population(stk.Population(amine2))

    exiter = stk.MoleculePresent(amine2)

    assert not exiter.exit(pop1)
    assert exiter.exit(pop2)


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

    exiter = stk.FitnessPlateau(2)

    assert not exiter.exit(pop)

    for mol in pop:
        mol.fitness = 2
    assert not exiter.exit(pop)

    pop2 = stk.Population(*pop)
    pop3 = stk.Population(pop2, pop2)
    assert exiter.exit(pop3)
