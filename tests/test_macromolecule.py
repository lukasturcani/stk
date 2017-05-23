from types import SimpleNamespace
from os.path import join
from ..molecular import (MacroMolecule, Molecule, FourPlusSix,
                         StructUnit, StructUnit2, StructUnit3)
from ..ga import Population

pop = Population.load(join('data', 'macromolecule', 'mm.json'))


def test_building_block_cores():
    # Check that the yielded rdkit molecules match the cores of the
    # building block molecules.
    macromol = pop[0]
    for i in range(len(macromol.building_blocks)):
        for frag in macromol.building_block_cores(i):
            bb1match = len(frag.GetSubstructMatch(
                           macromol.building_blocks[0]))
            bb2match = len(frag.GetSubstructMatch(
                           macromol.building_blocks[1]))
            nfrag_atoms = frag.GetNumAtoms()
            assert bb1match == nfrag_atoms or bb2match == nfrag_atoms


def test_bb_distortion():
    assert isinstance(pop[0].bb_distortion(), float)


def test_comparison():
    """
    Checks ``==``, ``>``, ``>=``, etc. operators.

    """

    # Generate cages with various fitnesses.
    a = MacroMolecule.testing_init('a', 'a', SimpleNamespace(a=1))
    a.fitness = 1

    b = MacroMolecule.testing_init('b', 'b', SimpleNamespace(b=1))
    b.fitness = 1

    c = MacroMolecule.testing_init('c', 'c', SimpleNamespace(c=1))
    c.fitness = 2

    # Comparison operators should compare their fitness.
    assert not a < b
    assert a <= b
    assert a == b
    assert c > b
    assert c >= a


def test_caching():

    # Make a MacroMolecule the regular way.
    bb1 = StructUnit2(join('data', 'struct_unit2', 'amine.mol2'))
    bb2 = StructUnit3(join('data', 'struct_unit3', 'amine.mol2'))
    mol1 = MacroMolecule([bb1, bb2], FourPlusSix())

    # Make a MacroMolecule using JSON.
    mol2 = pop[0]

    assert mol1 is not mol2

    # Remake the MacroMolecules.
    mol3 = Molecule.fromdict(mol1.json())
    mol4 = MacroMolecule(mol2.building_blocks,
                         mol2.topology.__class__())

    # Confirm they are cached.
    assert mol1 is mol3
    assert mol1 is not mol4
    assert mol2 is mol4
    assert mol2 is not mol3


def test_json_init():
    og_c = dict(MacroMolecule.cache)
    try:
        # Make a MacroMolecule using JSON.
        MacroMolecule.cache = {}
        mol = pop[0]

        assert mol.fitness is None
        assert all(isinstance(x, StructUnit) for x in
                   mol.building_blocks)
        assert isinstance(mol.key, tuple)
        assert mol.optimized is True
        assert mol.unscaled_fitness == {}
        assert mol.bonder_ids == [
                              6, 15, 24, 59, 68, 77, 112, 121, 130,
                              165, 174, 183, 219, 222, 252, 255, 285,
                              288, 318, 321, 351, 354, 384, 387]
        assert mol.energy.__class__.__name__ == 'Energy'
        assert mol.topology.__class__.__name__ == 'FourPlusSix'
        assert len(mol.mol.GetAtoms()) == 410
        assert mol.bonds_made == 12
        assert set(mol.bb_counter.values()) == {4, 6}
        assert mol.progress_params == {}
    finally:
        MacroMolecule.cache = og_c
