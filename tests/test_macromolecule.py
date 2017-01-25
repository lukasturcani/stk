from ..classes import MacroMolecule

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
