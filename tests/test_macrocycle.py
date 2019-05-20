import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk


def test_cycle_atoms():
    cycle = rdkit.MolFromSmiles('CCCC1CCCCCCCCC1')
    cycle = rdkit.AddHs(cycle)
    rdkit.EmbedMolecule(cycle)

    cycle_su = stk.MacrocycleStructUnit(cycle, [])
    assert set(cycle_su.cycle_atoms()) == set([3, 4, 5, 6, 7, 8,
                                               9, 10, 11, 12])


def test_cycle_coords(amine2, aldehyde2):
    cycle = stk.Macrocycle([amine2, aldehyde2],
                           stk.Cyclic('AB', [0, 0], 3))

    catoms = cycle.cycle_atoms()
    for i, coords in enumerate(cycle.cycle_coords(), 1):
        assert type(coords[0]) == int
        assert type(coords[1]) == float
        assert type(coords[2]) == float
        assert type(coords[3]) == float
        assert len(coords) == 4
    assert len(catoms) == i