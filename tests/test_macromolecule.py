from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ..classes import MacroMolecule, FourPlusSix, BuildingBlock, Linker, Cage, EightPlusTwelve


bb_file = next(x for x in get_mol_file() 
                                    if 'amine3f_14.mol' in x)
lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

bb = BuildingBlock(bb_file)
lk = Linker(lk_file)    
building_blocks = (bb, lk)
mol = MacroMolecule(building_blocks, FourPlusSix, 'you_can_delete_this.mol')

def test_caching():
    """
    Cages created with same arguments should return the same instance.

    """
    lk2 = Linker(lk_file)
    bb2 = BuildingBlock(bb_file)
    building_blocks2 = (lk2, bb2)
    mol2 = MacroMolecule(building_blocks2, FourPlusSix, 'you_can_delete_this2.mol')
    
    assert mol is mol2

    lk2 = Linker(lk_file)
    bb2 = BuildingBlock(bb_file)
    building_blocks2 = (lk2, bb2)
    mol2 = MacroMolecule(building_blocks2, EightPlusTwelve, 'you_can_delete_this2.mol')    
    assert mol is not mol2
    
    
    bb3_file = next(x for x in get_mol_file() 
                                    if 'amine3f_5.mol' in x)
    lk3_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_74.mol' in x)
                                        
    bb3 = BuildingBlock(bb3_file)
    lk3 = Linker(lk3_file)
    building_blocks3 = (bb3,lk3)
    mol3 = MacroMolecule(building_blocks3, FourPlusSix, 'you_can_delete_this3.mol')

    assert mol is not mol3
    assert mol2 is not mol3

def test_comparison():
    """
    Checks ``==``, ``>``, ``>=``, etc. operators.    
    
    """
    
    # Generate cages with various fitnesses.
    a = MacroMolecule.testing_init('a','a','a')
    a.fitness = 1
    
    b = MacroMolecule.testing_init('b', 'b', 'b')
    b.fitness = 1
    
    c = MacroMolecule.testing_init('c', 'c', 'c')
    c.fitness = 2
    
    # Comparison operators should compare their fitness.
    assert not a < b
    assert a <= b
    assert a == b
    assert c > b
    assert c >= a

def test_same():
    """
    Tests the `same_cage` method.    
    
    Cages initialized from the same arguments should return ``True`` 
    through this method, even if the ``Cage`` class stops being cached.    
    
    """
    
    a = MacroMolecule.testing_init('a', 'b', 'c')
    b = MacroMolecule.testing_init('a', 'a', 'b')
    c = MacroMolecule.testing_init('a', 'a', 'b')
    d = MacroMolecule.testing_init('a', 'b', 'a')

    assert not a.same(b)
    assert b.same(c)    
    assert c.same(b)
    assert not d.same(c)

def test_get_heavy_as_graph():
    graph = mol.get_heavy_as_graph()
    expected_nodes = 224
    expected_edges = 244

    assert len(graph.nodes()) == expected_nodes
    assert len(graph.edges()) == expected_edges

        
def test_heavy_get_atom_coords():
    
    expected_output = (1, 52, 1)
    output = tuple(int(x) for x in mol.heavy_get_atom_coords(21))
    assert expected_output == output

def test_heavy_distance():

    assert int(mol.heavy_distance(22, 51)) == 62
    
def test_get_atom_distances():
    assert len(list(mol.get_heavy_atom_distances())) == 276
    for x,y,z in mol.get_heavy_atom_distances():
        assert isinstance(x, float)
        assert x >= 0
        assert isinstance(y, int)
        assert isinstance(z, int)

