import stk
import pytest


def _population_members():
    bb1 = stk.BuildingBlock('NC(CCO)CN', ['amine'])
    bb2 = stk.BuildingBlock('[Br]CCCC[Br]', ['bromine'])
    bb3 = stk.BuildingBlock('[I]COCC[I]', ['iodine'])
    bb4 = stk.BuildingBlock('O=CC(C=O)CC=O', ['aldehyde'])

    c1 = stk.ConstructedMolecule(
        building_blocks=[bb2],
        topology_graph=stk.polymer.Linear('A', 3)
    )
    c2 = stk.ConstructedMolecule(
        building_blocks=[bb1, bb4],
        topology_graph=stk.cage.FourPlusSix()
    )
    c3 = stk.ConstructedMolecule(
        building_blocks=[bb1, bb4],
        topology_graph=stk.cage.EightPlusTwelve()
    )
    c4 = stk.ConstructedMolecule(
        building_blocks=[bb2, bb3],
        topology_graph=stk.polymer.Linear('AB', 3)
    )

    return bb1, bb2, bb3, bb4, c1, c2, c3, c4


@pytest.fixture('session')
def population():
    bb1, bb2, bb3, bb4, c1, c2, c3, c4 = _population_members()
    return stk.Population(
        c3,
        stk.BuildingBlock('NCCCN'),
        bb1,
        stk.BuildingBlock('NCCCN', ['amine']),
        c1,
        stk.BuildingBlock('O=CCC=O'),
        c1,
        stk.BuildingBlock('O=CCC=O', ['aldehyde']),
        stk.Population(
            c1,
            stk.BuildingBlock('[Br]CC[Br]'),
            bb1,
            bb2,
            bb1
        ),
        c2,
        stk.Population(
            bb1,
            stk.BuildingBlock('CCCC'),
            stk.Population(
                bb3,
                stk.BuildingBlock('NNNN'),
                c4
            )
        )
    )


@pytest.fixture(scope='function')
def tmp_population():
    bb1, bb2, bb3, bb4, c1, c2, c3, c4 = _population_members()
    return stk.Population(
        c3,
        stk.BuildingBlock('NCCCN'),
        bb1,
        stk.BuildingBlock('NCCCN', ['amine']),
        c1,
        stk.BuildingBlock('O=CCC=O'),
        c1,
        stk.BuildingBlock('O=CCC=O', ['aldehyde']),
        stk.Population(
            c1,
            stk.BuildingBlock('[Br]CC[Br]'),
            bb1,
            bb2,
            bb1
        ),
        c2,
        stk.Population(
            bb1,
            stk.BuildingBlock('CCCC'),
            stk.Population(
                bb3,
                stk.BuildingBlock('NNNN'),
                c4
            )
        )
    )


@pytest.fixture('session')
def generation():
    pop = stk.Population(
        *(stk.BuildingBlock('C') for i in range(10))
    )
    for i, mol in enumerate(pop):
        mol.fitness = i
    return pop


@pytest.fixture('session')
def progress():
    pop = stk.Population()
    for i in range(15):
        subpop = stk.Population(
            *(stk.BuildingBlock('C') for j in range(5))
        )
        pop.subpopulations.append(subpop)

    for i, mol in enumerate(pop):
        mol.fitness = i

    return pop
