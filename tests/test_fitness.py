from ..ga import cage, cage_target, cage_c60
from ..molecular import Molecule
from os.path import join
import numpy as np

cagemol = Molecule.load(join('data', 'fitness', 'cage.json'))

def test_cage():
    assert np.allclose(cage(cagemol),
        [17.32137896, 14.87621275, 0.71956529, -133.76115919], atol=1)
