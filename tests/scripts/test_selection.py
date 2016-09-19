from operator import attrgetter
import numpy as np
from collections import Counter

from .test_population import generate_population
from ...classes import GATools, Selection, FunctionData

labels = "abcdefghijklmnopqrstuvwxyz"
# Make a population.
pop1 = generate_population()


        
def test_selection_population_integration():
    """
    Tests to see if selection classes is integrated correctly with the
    population class.
    
    """

    #Set fitnesses and labels

    for index, mol in enumerate(pop1):
        mol.fitness = index
        mol.label = labels[index]

    # Make a ``GATools`` attribute and give it to the population.
    fittest = FunctionData('fittest')
    fittest2 = FunctionData('fittest')
    fittest3 = FunctionData('fittest')
    
    selector = Selection(fittest, fittest2, fittest3)
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = pop1.select('generational')
    pop3 = pop1.select('mating')
    pop4 = pop1.select('mutation')

    pop1_members = sorted(pop1, reverse=True)

    assert list(pop2) == pop1_members
    assert list(pop3) ==  pop1_members
    assert list(pop4) == pop1_members
    
    
def test_fittest():
    fittest = FunctionData('fittest')
    selector = Selection(fittest, fittest, fittest)
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools

    # Set fitnesses. 
    for mol in pop1:
        mol.fitness = np.random.randint(0, 100)

    pop2 = list(pop1.select())
    assert len(pop2) == 22
    
    for i, mol in enumerate(pop2):
        if i == 0:
            continue
        
        mol2 = pop2[i-1]
        assert mol2 >= mol
        
def test_all_combinations():
    all_combs = FunctionData('all_combinations')
    selector = Selection('a', all_combs, 'b',)
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    #Set fitnesses and labels
    for index, mol in enumerate(pop1):
        mol.fitness = 1/(index+1)
        mol.label = labels[index]
    
    pop2 = list(pop1.select('mating'))
    assert len(pop2) == 231
    for parents in pop2:
        assert len(parents) == 2
    

def test_all_combinations_n_fittest():
    all_combs = FunctionData('all_combinations_n_fittest', n=5)
    fittest = FunctionData('fittest')
    
    selector = Selection('a', all_combs, fittest)
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    #Set fitnesses and labels
    for index, mol in enumerate(pop1):
        mol.fitness = 1/(index+1)
        mol.label = labels[index]
    
    pop2 = list(pop1.select('mating'))
    n_fittest = list(pop1.select('mutation'))[:5]
    
    assert len(pop2) == 10
    for parents in pop2:
        assert len(parents) == 2
        assert parents[0] in n_fittest
        assert parents[1] in n_fittest
        
def test_roulette():
    
    # Set fitnesses. 
    for mol in pop1:
        mol.fitness = np.random.randint(1, 100)

    # Calculate probabilities of roulette selection
    total_fitness = sum(x.fitness for x in pop1)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness

    # Conditions 1: elitism = False, truncation = False, 
    #               duplicates = False
    roulette = FunctionData('roulette')
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = list(pop1.select())    
    pop2_dedupe = set(pop2)    
    
    assert len(pop2) == len(pop1)
    assert len(pop2) == len(pop2_dedupe)
    
    # Conditions 2: elitism = False, truncation = False, 
    #               duplicates = True
    roulette = FunctionData('roulette', duplicates=True)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i > 500:
            break
        pop2.append(x)
            
    for ind in pop1:
        assert ind in pop2

    # Conditions 3: elitism = False, truncation = True,
    #               duplicates = False
    roulette = FunctionData('roulette', truncation=5)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = list(pop1.select())
    pop2_dedupe = set(pop2)
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:5]
    assert len(pop2) == 5
    assert len(pop2) == len(pop2_dedupe)
    for ind in pop2:
        assert ind in fittest
        
    # Conditions 4: elitism = False, truncation = True, 
    #               duplicates = True
    
    # re-Calculate probabilities of roulette selection
    total_fitness = sum(x.fitness for x in pop1 if x in fittest)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness

    roulette = FunctionData('roulette', duplicates=True, truncation=5)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i > 500:
            break
        pop2.append(x)
    
    for ind in fittest:
        assert ind in pop2

    # Conditions 5: elitism = True, truncation = False, 
    #               duplicates = False
    roulette = FunctionData('roulette', elitism=3)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = list(pop1.select())
    pop2_dedupe = set(pop2)
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    assert len(pop2) == 22
    assert len(pop2) == len(pop2_dedupe)
    for ind in pop2[:3]:
        assert ind in fittest
        
    # Conditions 6: elitism = True, truncation = False, 
    #               duplicates = True
    roulette = FunctionData('roulette', elitism=3, duplicates=True)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools

    # re-Calculate probabilities of roulette selection
    total_fitness = sum(x.fitness for x in pop1)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness
 
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i > 500:
            break
        if i == 4:
            for ind in fittest:
                assert ind in pop2
        pop2.append(x)
            
    for ind in pop1:
        assert ind in pop2
    
    # Conditions  7: elitism = True, truncation = True,
    #               duplicates = False
    roulette = FunctionData('roulette', elitism=3, truncation=18)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = list(pop1.select())
    pop2_dedupe = set(pop2)
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    assert len(pop2) == 18
    assert len(pop2) == len(pop2_dedupe)
    for ind in pop2[:3]:
        assert ind in fittest        
    
    # Conditions 8: elitism = True, truncation = True,
    #               duplicates = True

    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:10]
    total_fitness = sum(x.fitness for x in pop1 if x in fittest)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness

    roulette = FunctionData('roulette', duplicates=True, truncation=10,
                            elitism=3)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i > 500:
            break
        if i == 3:
            for ind in pop2:
                assert ind in fittest[:3]
        pop2.append(x)
           
    for ind in fittest:
        assert ind in pop2

    
def test_mating_roulette():
    # Set fitnesses. 
    for mol in pop1:
        mol.fitness = np.random.randint(1, 100)

    # No truncation.
    # Calculate probabilities of roulette selection
    total_fitness = sum(x.fitness for x in pop1)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness    

    roulette = FunctionData('mating_roulette')
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, parents in enumerate(pop1.select()):
        if i == 500:
            assert len(parents) == 2
            break
        pop2.extend(parents)
        
    for ind in pop1:
        assert ind in pop2

    

    # With truncation.
    # Calculate probabilities of roulette selection
    fittest = sorted(pop1, key=attrgetter('fitness'), reverse=True)[:10]
    total_fitness = sum(x.fitness for x in pop1 if x in fittest)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness    

    roulette = FunctionData('mating_roulette', truncation=10)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, parents in enumerate(pop1.select()):
        if i == 500:
            assert len(parents) == 2
            break
        pop2.extend(parents)
        
    for ind in fittest:
        assert ind in pop2
        
        
        
        