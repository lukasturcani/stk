from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ..classes import Cage, FourPlusSix

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
    a = Cage('a','a','a')
    a.fitness = 1
    
    b = Cage('b', 'b', 'b')
    b.fitness = 1
    
    c = Cage('c', 'c', 'c')
    c.fitness = 2
    
    # Comparison operators should compare their fitness.
    assert not a < b
    assert a <= b
    assert a == b
    assert c > b
    assert c >= a

def test_same_cage():
    """
    Tests the `same_cage` method.    
    
    Cages initialized from the same arguments should return ``True`` 
    through this method, even if the ``Cage`` class stops being cached.    
    
    """
    
    a = Cage('a', 'b', 'c')
    b = Cage('a', 'a', 'b')
    c = Cage('a', 'a', 'b')
    d = Cage('a', 'b', 'a')

    assert not a.same_cage(b)
    assert b.same_cage(c)    
    assert c.same_cage(b)
    assert not d.same_cage(c)
    
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
    
    cage = Cage(bb_file, lk_file, FourPlusSix, 
                    'you_can_delete_this.mol')
    
    assert hasattr(cage, 'prist_mol_file')
    assert hasattr(cage, 'heavy_mol_file')
    assert hasattr(cage, 'prist_mol')
    assert hasattr(cage, 'heavy_mol')
    assert hasattr(cage, 'prist_smiles')
    assert hasattr(cage, 'heavy_smiles')
    assert hasattr(cage, 'topology')
    
    
    
    
    
    