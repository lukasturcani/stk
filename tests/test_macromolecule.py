import os

from ..classes import MacroMolecule, StructUnit3, StructUnit2
from ..classes.topology import FourPlusSix, EightPlusTwelve

bb_file = os.path.join('data', 'macromolecule', 
                       'macromolecule_bb_amine.mol')
lk_file = os.path.join('data', 'macromolecule', 
                       'macromolecule_lk_aldehyde_1.mol')

def test_caching():
    """
    Cages created with same arguments should return the same instance.

    """
    try:
        bb = StructUnit3(bb_file)
        lk = StructUnit2(lk_file)    
        building_blocks = (bb, lk)
    
        mol = MacroMolecule(building_blocks, FourPlusSix, 
                                     'delete_me0.mol')     
        
        # First init a new MacroMolecule using the same values as prviously.
        # This should yield the same MacroMolecule.
        mol2 = MacroMolecule(building_blocks, FourPlusSix, 
                                     'delete_me1.mol')    
        assert mol is mol2
    
        # Change the order of the building blocks in the tuple, this
        # should still be the same MacroMolecule.
        building_blocks2 = (lk, bb)
        mol2 = MacroMolecule(building_blocks2, FourPlusSix, 
                             'delete_me2.mol')    
        assert mol is  mol2
        
        # Change the topology, this should make a new MacroMolecule.     
        mol2 = MacroMolecule(building_blocks, EightPlusTwelve,
                             'delete_me3.mol')    
        assert mol is not mol2
        
        
        # Change one of the building blocks, this should make a new 
        # MacroMolecule
        lk2_file = os.path.join('data', 'macromolecule',
                                'macromolecule_lk_aldehyde_2.mol')
                                            
        lk2 = StructUnit2(lk2_file)
        building_blocks2 = (bb, lk2)
        mol2 = MacroMolecule(building_blocks2, FourPlusSix, 
                              'delete_me4.mol')
    
        assert mol is not mol2
        
    finally:
        for file in os.listdir():
            if 'delete_me' in file:
                os.remove(file)


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
