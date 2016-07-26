from ..classes import Population, Cage, GATools
import pytest
import itertools as it
from collections import Counter

def generate_population():
    """
    Returns a population of subpopulations and direct members.
   
    """
    #Generate a bunch of cages.
    cages = [Cage.__new__(Cage) for x in range(0,22)]    
    
    # Generate a couple of populations to be used as subpopulations.
    sub1 = Population(*cages[0:4])
    sub2 = Population(*cages[4:9])
    sub3 = Population(*cages[9:14])
    sub4 = Population(*cages[14:19])                      
      
    # Place subpopulations into one another.
    sub1.populations.append(sub3)
    sub2.populations.append(sub4)
    
    # Initialize final population of subpopulations and cages.
    return Population(sub1, sub2, *cages[-3:])
    

def test_init():
    """
    Tests the __init__ method of the Population class. 
    
    """
        
    # A test where only ``Cage`` and ``Population`` instances are used
    # for initialization as well as a single ``GATools`` instnace, as 
    # intended. Should pass.       
    Population(Cage.__new__(Cage), Cage.__new__(Cage),
               Population.__new__(Population), Cage.__new__(Cage), 
               Population.__new__(Population), GATools.__new__(GATools))

    # Only ``Cage`` and ``Population`` instances are used for 
    # initialization, no ``GATools`` instance. Valid, should pass.
    Population(Cage.__new__(Cage), Cage.__new__(Cage),
               Population.__new__(Population), Cage.__new__(Cage), 
               Population.__new__(Population))
    
    # Multiple ``GATools`` instances. Invalid, should raise TypeError.           
    with pytest.raises(TypeError):
       Population(Cage.__new__(Cage), Cage.__new__(Cage), 
                  Population.__new__(Population), Cage.__new__(Cage),
                  Population.__new__(Population), 
                  GATools.__new__(GATools), GATools.__new__(GATools))
                  
    # Non ``GATools``, ``Cage`` or ``Population`` instance used for
    # initialization (``int``). Invalid, should raise TypeError.
    with pytest.raises(TypeError):        
        Population(Cage.__new__(Cage), Cage.__new__(Cage), 
                   Population.__new__(Population), Cage.__new__(Cage),
                   Population.__new__(Population), 
                   GATools.__new__(GATools), 12)

    # Due to a technicality (see __init__ source code) any number of 
    # ``None`` type arguments can be supplied. Valid, should pass.
    Population(*[None for x in range(0,10)])
    
def test_all_members():
    """
    Check that all members, direct and in subpopulations, are returned.

    """    
    
    #Generate a bunch of cages.
    cages = [Cage.__new__(Cage) for x in range(0,22)]    
    
    # Generate a couple of ``Populations`` to be used as subpopulations.
    sub1 = Population(*cages[0:4])
    sub2 = Population(*cages[4:9])
    sub3 = Population(*cages[9:14])
    sub4 = Population(*cages[14:19])                      
      
    # Place subpopulations in one another.
    sub1.populations.append(sub3)
    sub2.populations.append(sub4)
    
    # Initialize main population from subpopulations and cages.
    main = Population(sub1, sub2, *cages[-3:])

    # Place the results of ``all_members`` into a list.
    all_members = [cage for cage in main.all_members()]

    # Check that each generated cage is in `all_members`. Should pass.
    assert all(cage in all_members for cage in cages)  
    
    # Add a cage to `cages`. Now there should be a cage in `cages`, not 
    # present in main. Should fail.
    cages.append(Cage.__new__(Cage))      
    with pytest.raises(AssertionError):
        assert all(cage in all_members for cage in cages)
        
def test_add_members_duplicates():
    """
    Members in population added to `members` of the other.
    
    Duplicate additions should be allowed.    
    
    Because uninitialized cages are used and ``Cages`` compare with the 
    `is` operator on the basis of their linkers and building blocks 
    rather than on the their `__id__` the `is in` operator is not used. 
    This is because the cages have no `bb` or `lk` attributes as they 
    are not initialized.
    
    Instead, lists are consisting of the id values of the cages are 
    generated and those are compared with the `is in` operator.
    
    """
    
    # Create a population to be added and one to be added to.
    receiver = generate_population()
    receiver_ids = [id(cage) for cage in receiver]
    
    supplier = generate_population()
    supplier_ids = [id(cage) for cage in supplier]

    # Cages in `supplier` should not be in `reciever` yet.
    assert all(id_val not in receiver_ids for id_val in supplier_ids)
    # Add all cages in `supplier` to `receiver`.
    receiver.add_members(supplier, duplicates=True)
    # Recalcualte id list from `members` attribute as this is the 
    # intended destination.
    receiver_ids = [id(cage) for cage in receiver.members]
    # Cages in `supplier` should now be held in `members` attribute of 
    # `receiver`.
    assert all(id_val in receiver_ids for id_val in supplier_ids)
                      
    # Add supplier cages a second time, demonstarting that duplicates
    # are allowed.
    receiver.add_members(supplier, duplicates=True)
    
    # Count the frequency of each `__id__` in `receiver.members`.
    count = Counter(receiver.members)
    # Ensure that each of suppliers ids is present twice.
    assert all(count[cage] == 2 for cage in supplier)
    
def test_add_members_no_duplicates():
    """
    Members in population added to `members` of another.
    
    Duplicate additions not allowed.

    Duplicates are found via the ``==`` operator on the `bb` and `lk`
    attributes of a cage. The attributes here are initialized to 
    simple strings - not indicative of the what the attributes will be
    set to in a realistic setting. Likely a `BuildingBlock` class etc.    

    """

    values = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'    
    
    # Generate an initial populaiton, initialize the cage's `bb` and 
    # `lk` attributes.
    receiver = generate_population()
    for index, cage in enumerate(receiver):
        cage.bb = values[index]
        cage.lk = values[index]
        
    # Note size of receiver population.
    receiver_size = len(receiver)        
        
    # Same as above but with another population. Note the objects are 
    # different instances, but their `lk` and `bb` attributes are the 
    # same. As a result they should compare equal to those in 
    # `receiver`.
    supplier_same = generate_population()
    for index, cage in enumerate(supplier_same):
        cage.bb = values[index]
        cage.lk = values[index]
        
    # Add supplier to the receiver. None of the suppliers cages should
    # be added and therefore the len of supplier shoudl be the same as 
    # at the start.
    receiver.add_members(supplier_same)

    assert receiver_size == len(receiver)
    
    # Generate another population. This time the 'bb' and 'lk' of the 
    # cages will have different combinations to the receiver population.
    supplier_different = generate_population()
    for index, cage in enumerate(supplier_different):
        cage.bb = values[index]
        cage.lk = values[index+1]    

    # Add `supplier_different` to `receiver`. All of the cages should be
    # addable as none should be duplicates. The size of the `receiver`
    # population should increase by the size of the `supplier_different`    
    # pop.
    receiver.add_members(supplier_different)
    assert receiver_size + len(supplier_different) == len(receiver)
    
def test_getitem():
    """
    Test that the '[]' operator is working.
        
    """
    
    # [:] should flatten the population, all members from subpopulations
    # should be transferred to the members attribute.    
    pop = generate_population()
    flat_pop = pop[:]
    assert all(cage in flat_pop.members for cage in pop)
    # Verify lack of subpopulations.
    assert not flat_pop.populations
    # An integer index should return a ``Cage`` instance.
    assert type(pop[5]) is Cage
    # Non integer/slice indices are not supported
    with pytest.raises(TypeError):
        pop[5.5]

def test_
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
    
    
    