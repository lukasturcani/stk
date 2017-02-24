from ..ga import Crossover, Population
from os.path import join

Population.load(join('data', 'crossover', 'molecules.json'))

def test_bb_lk_exchange():
    offspring = Crossover.bb_lk_exchange(None, p1, p2)
    for p in (p1, p2):
        assert p not in offspring
