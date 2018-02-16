from ..ga import Crossover, Population
from ..molecular import Molecule
from os.path import join

pop = Population.load(join('data', 'crossover', 'molecules.json'),
                      Molecule.from_dict)
p1 = pop[0]
p2 = pop[1]


def test_bb_lk_exchange():
    offspring = Crossover.bb_lk_exchange(None, p1, p2)
    for p in (p1, p2):
        assert p not in offspring
