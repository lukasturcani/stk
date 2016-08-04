from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ..classes import MacroMolecule, FourPlusSix, BuildingBlock, Linker, Cage


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
    
    pop1 = generate_population()
    pop2 = generate_population()
    pop3 = generate_population(offset=True)
    
    for cage1, cage2, cage3 in zip(pop1, pop2, pop3):
        assert cage1 is cage2
        assert cage1 is not cage3
        assert cage2 is not cage3

def test_comparison():
    """
    Checks ``==``, ``>``, ``>=``, etc. operators.    
    
    """
    
    # Generate cages with various fitnesses.
    a = MacroMolecule('a','a','a', 1)
    a.fitness = 1
    
    b = MacroMolecule('b', 'b', 'b', 1)
    b.fitness = 1
    
    c = MacroMolecule('c', 'c', 'c', 1)
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
    
    a = MacroMolecule('a', 'b', 'c', 'd')
    b = MacroMolecule('a', 'a', 'b', 'd')
    c = MacroMolecule('a', 'a', 'b', 'd')
    d = MacroMolecule('a', 'b', 'a', 'd')

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
    
    expected_output = (5, 50, -1)
    output = tuple(int(x) for x in mol.heavy_get_atom_coords(21))
    assert expected_output == output

def test_heavy_distance():

    assert int(mol.heavy_distance(22, 51)) == 54
    
def test_get_atom_distances():
    assert len(list(mol.get_heavy_atom_distances())) == 276
    for x,y,z in mol.get_heavy_atom_distances():
        assert isinstance(x, float)
        assert x >= 0
        assert isinstance(y, int)
        assert isinstance(z, int)

