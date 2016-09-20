from collections import Counter
import os

from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ...classes import (GATools, Selection, 
                        FunctionData, BuildingBlock, Mutation)
from ...classes import  (Linker, Cage, FourPlusSix, 
                         EightPlusTwelve, SixPlusNine, Population)
from ...convenience_functions import flatten

bb_file = next(x for x in get_mol_file() 
                                    if 'amine3f_14.mol' in x)
lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

bb = BuildingBlock(bb_file)
lk = Linker(lk_file)    
building_blocks = (bb, lk)
mol = Cage(building_blocks, FourPlusSix, 
           'you_can_delete_this3.mol')

bb2_file = next(x for x in get_mol_file() 
                                    if 'amine3f_5.mol' in x)
lk2_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_28.mol' in x) 

bb2 = BuildingBlock(bb2_file)
lk2 = Linker(lk2_file)    
building_blocks2 = (bb2, lk2)
mol2 = Cage(building_blocks2, EightPlusTwelve, 
             'you_can_delete_this4.mol')


def test_random_mutation_function_selection():
    # This tests if the selection of a mutation function by the Mutation
    # object is working properly.

    # The original population, the one being mutated, will have only
    # a single topology and bb. This means when checking the mutants,
    # about half should have a new topology and about half should have 
    # a new bb. When no weights are used. With weights the ratio should
    # be different


    # Make a ``GATools`` attribute and give it to the population.
    database = os.path.join(os.getcwd(), 
                            'Database_prec', 'amines3f')
    rand_bb = FunctionData('random_bb', database=database)
            
    rand_topology = FunctionData('random_cage_topology', 
            topologies=[FourPlusSix, EightPlusTwelve, SixPlusNine])
    
    sel = FunctionData('roulette', duplicates=True)
    
    selector = Selection('a', 'a', sel)
    
    # No weights first.    
    
    mutator = Mutation([rand_bb, rand_topology], 20)
    ga_tools = GATools(selector, 'a', mutator, 
                       'rdkit_optimization', 'cage')

    # Generate the original population
    lk_db = os.path.join(os.getcwd(), 'Database_prec', 'aldehydes2f')

    pop1 = Population.init_fixed_bb_cages(bb_file, lk_db, [FourPlusSix],
                                          20, ga_tools)    

    pop2 = pop1.gen_mutants()
    new_top = 0
    new_bb = 0
    for mutant in pop2:
        if type(mutant.topology) != FourPlusSix:
            new_top += 1
        bb = next(x.prist_mol_file for x in mutant.building_blocks if 
                    isinstance(x, BuildingBlock))
        if bb != bb_file:
            new_bb +=1
            
    assert new_top + new_bb == 20
    

    

def test_random_bb():
    
    mol.fitness = 1
    mol2.fitness = 2
    # Make a ``GATools`` attribute and give it to the population.
    database = os.path.join(os.getcwd(), 
                            'Database_prec', 'amines3f')
    rand_bb = FunctionData('random_bb', database=database)
    sel = FunctionData('fittest')
    
    selector = Selection('a', 'a', sel)
    mutator = Mutation([rand_bb], 1)
    ga_tools = GATools(selector, 'a', mutator, 
                       'rdkit_optimization', 'cage')


    pop1 = Population(ga_tools, mol, mol2)
    pop2 = pop1.gen_mutants()

    for mol3 in pop2:
        assert mol3 not in pop1
    
    assert len(pop2) == 1
    bb = next(x for x in pop2[0].building_blocks if isinstance(x, BuildingBlock))
    lk = next(x for x in pop2[0].building_blocks if isinstance(x, Linker))
    assert bb not in mol2.building_blocks
    assert lk in mol2.building_blocks

def test_random_cage_topology():   
    mol.fitness = 1
    mol2.fitness = 1
    # Make a ``GATools`` attribute and give it to the population.
    rand_topology = FunctionData('random_cage_topology', 
            topologies=[FourPlusSix, EightPlusTwelve, SixPlusNine])
    
    sel = FunctionData('roulette', duplicates=True)
    
    selector = Selection('a', 'a', sel)
    mutator = Mutation([rand_topology], 20)
    ga_tools = GATools(selector, 'a', mutator, 
                       'rdkit_optimization', 'cage')


    pop1 = Population(ga_tools, mol, mol2)
    pop2 = pop1.gen_mutants()
    
    # Only new cages should be present in the generated population.
    for mol3 in pop2:
        assert mol3 not in pop1
    
    # I have two cages, each cage can make two mutants. This is
    # because there are 2 topologies for each cage available to make
    # a new mutant.        
    
    assert len(pop2) == 4
    
    # Check that the correct topologies were substituted.
    for mutant in pop2:
        parent = next(x for x in pop1 if 
                    x.building_blocks == mutant.building_blocks)
        assert type(parent.topology) != type(mutant.topology)
