from ..molecular import Cage, Molecule
from os.path import join

def test_window_difference():
    path = join('data', 'cage', 'cage.json')
    with open(path, 'r') as f:
        mol = Molecule.load(f.read())

    print(mol.window_difference())
    mol.write('/home/lukas/hi.mol')
    assert False

def test_windows():
    ...
