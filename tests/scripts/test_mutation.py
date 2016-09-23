import os
import pickle

from ...classes import (GATools, Selection, 
                        FunctionData, BuildingBlock, Mutation)
from ...classes import  (Linker, FourPlusSix, EightPlusTwelve, 
                         SixPlusNine, Population)

mol_file = os.path.join('data', 'mutation', 'cage1')
mol2_file = os.path.join('data', 'mutation', 'cage2')

with open(mol_file, 'rb') as dump_file:
    mol = pickle.load(dump_file)
    
with open(mol2_file, 'rb') as dump_file:
    mol2 = pickle.load(dump_file)

def test_random_mutation_function_selection():
    # This tests if the selection of a mutation function by the Mutation
    # object is working properly.

    # The original population, the one being mutated, will have only
    # a single topology and bb. This means when checking the mutants,
    # about half should have a new topology and about half should have 
    # a new bb. When no weights are used. With weights the ratio should
    # be different

    # Load the original population
    pop_file = os.path.join('data', 'mutation', 'single_bb_pop')
    pop1 = Population.load(pop_file)
    bb_file = next(x.prist_mol_file for x in pop1[0].building_blocks 
                                        if isinstance(x, BuildingBlock))

    # Create the mutants.
    pop2 = pop1.gen_mutants()
    
    # The original weights dicatate the only the mutation of bbs will
    # take place.
    new_top = 0
    new_bb = 0
    for mutant in pop2:
        if type(mutant.topology) != FourPlusSix:
            new_top += 1
        bb = next(x.prist_mol_file for x in mutant.building_blocks if 
                    isinstance(x, BuildingBlock))
        if bb != bb_file:
            new_bb +=1
            
    assert new_top == 0
    assert new_bb != 0  
    
    # Reverse the weights.
    pop1.ga_tools.mutation.weights = [0,1]
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
            
    assert new_top != 0
    assert new_bb == 0
    
    # Split the weights.
    pop1.ga_tools.mutation.weights = [0.5,0.5]
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
            
    assert new_top != 0
    assert new_bb != 0    
   
def test_random_bb():
    
    mol.fitness = 1
    mol2.fitness = 2
    # Make a ``GATools`` attribute and give it to the population.
    database = os.path.join('data', 'mutation', 'bb_db')
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
