from ..ga import cage, cage_target, cage_c60
from ..molecular import Molecule
from os.path import join

cagemol = Molecule.load(join('data', 'fitness', 'cage.json'))

def test_cage():
    cage(cagemol, 8)

    assert cagemol.unscaled_fitness == [1,1,1,1]

def test_cage_target():
    ...

def test_cage_c60():
    ...
