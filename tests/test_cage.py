from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ..classes import Cage, FourPlusSix, BuildingBlock, Linker

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
    a = Cage('a','a','a', 1)
    a.fitness = 1
    
    b = Cage('b', 'b', 'b', 1)
    b.fitness = 1
    
    c = Cage('c', 'c', 'c', 1)
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
    
    a = Cage('a', 'b', 'c', 'd')
    b = Cage('a', 'a', 'b', 'd')
    c = Cage('a', 'a', 'b', 'd')
    d = Cage('a', 'b', 'a', 'd')

    assert not a.same(b)
    assert b.same(c)    
    assert c.same(b)
    assert not d.same(c)
    
def test_init():
    """
    Ensure that cages initialize attributes successfully.
    
    This function uses the ``aldehyde2f_3.mol`` and
    ``amine3f_14.mol`` files.
    
    Only presence of the attributes is tested here. Testing whether the
    correct cage was built should be done in the topology testing
    module.

    """ 
    
    bb_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)
    lk_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)    
    
    bb = BuildingBlock(bb_file)
    lk = Linker(lk_file)    
    building_blocks = (bb, lk)
    cage = Cage(building_blocks, FourPlusSix, 'you_can_delete_this.mol')
    
    assert hasattr(cage, 'prist_mol_file')
    assert hasattr(cage, 'heavy_mol_file')
    assert hasattr(cage, 'prist_mol')
    assert hasattr(cage, 'heavy_mol')
    assert hasattr(cage, 'prist_smiles')
    assert hasattr(cage, 'heavy_smiles')
    assert hasattr(cage, 'topology')
    
    
    
    
    
    