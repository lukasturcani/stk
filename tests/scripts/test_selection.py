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
        mol.fitness = np.random.randint(1, 25)

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
    
    # Select members of the population
    pop2 = list(pop1.select())    
    pop2_dedupe = set(pop2)    
    
    # No duplicates are allowed so lengths of each population should be
    # the same.
    assert len(pop2) == len(pop1)
    assert len(pop2) == len(pop2_dedupe)
    
    # Uncomment the following lines to check distribution.
#    count1 = Counter(pop1)
#    count2 = Counter(pop2)
#    print(count1, len(count1))
#    print(count2, len(count2))
#    print(pop2)
#    assert False
    
    
    # Conditions 2: elitism = False, truncation = False, 
    #               duplicates = True
    roulette = FunctionData('roulette', duplicates=True)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i > 2000:
            break
        pop2.append(x)

    count = Counter(pop2)
    # The biggest count should be more than 1 as duplicates are allowed.
    assert max(count.values()) > 1
    # No truncation means that every population members should be
    # present at least once.
    pop2_dedupe = set(pop2)
    assert len(pop2_dedupe) == len(pop1)

    # Uncomment the following lines to check distribution.
#    count1 = Counter(pop1)
#    count2 = Counter(pop2)
#    print(count1, len(count1))
#    print(count2, len(count2))
#    assert False

    # Conditions 3: elitism = False, truncation = True,
    #               duplicates = False
    roulette = FunctionData('roulette', truncation=17)
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
    # Uncomment the following lines to check distribution.
#    count1 = Counter(pop1)
#    count2 = Counter(pop2)
#    print(count1, len(count1))
#    print(count2, len(count2))
#    print(pop2)
#    assert False        
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
        if i > 1000:
            break
        pop2.append(x)

    count = Counter(pop2)
    assert list(count.values())[0] > 1
    pop2_dedupe = set(pop2)
    assert len(pop2_dedupe) == 17
    # Uncomment the following lines to check distribution.
#    count1 = Counter(pop1)
#    count2 = Counter(pop2)
#    print(count1, len(count1))
#    print(count2, len(count2))
#    assert False

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
        if i == 4:
            for ind in fittest:
                assert ind in pop2
            break
        pop2.append(x)

    
    # Conditions  7: elitism = True, truncation = True,
    #               duplicates = False
    roulette = FunctionData('roulette', elitism=3, truncation=17)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = list(pop1.select())
    pop2_dedupe = set(pop2)
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:5]
    assert len(pop2) == 5
    assert len(pop2) == len(pop2_dedupe)
    for ind in pop2[:3]:
        assert ind in fittest        
    
    # Conditions 8: elitism = True, truncation = True,
    #               duplicates = True

    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:10]
    total_fitness = sum(x.fitness for x in pop1 if x in fittest)
    for mol in pop1:
        mol.p = mol.fitness / total_fitness

    roulette = FunctionData('roulette', duplicates=True, truncation=12,
                            elitism=3)
    selector = Selection(roulette, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools
    
    pop2 = []
    for i, x in enumerate(pop1.select()):
        if i == 2:
            for ind in pop2:
                assert ind in fittest[:3]
            break
        pop2.append(x)
    
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
        if i == 50:
            assert len(parents) == 2
            break
        pop2.extend(parents)
        
    count = Counter(pop2)
    assert max(count.values()) > 1

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
        if i == 50:
            assert len(parents) == 2
            break
        pop2.extend(parents)
        
    count = Counter(pop2)
    assert max(count.values()) > 1
        
def test_deterministic_sampling():
    # Set fitnesses. 
    for mol in pop1:
        mol.fitness = np.random.randint(1, 25)

    # Conditions 1: elitism = False, truncation = False, 
    #               duplicates = False
    deterministic = FunctionData('deterministic_sampling')
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    
    pop2 = list(pop1.select())
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
            
    for ind in pop1:
        assert ind in pop2
        
    # Conditions 2: elitism = False, truncation = False, 
    #               duplicates = True
    deterministic = FunctionData('deterministic_sampling', 
                                 duplicates=True)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average to verify that
    # multiple yields are taking place.
        
    pop1[0].fitness = 100
    pop2 = list(pop1.select())
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) > 1
    # Check for no truncation.   
    for ind in pop1:
        assert ind in pop2    
   
   # Conditions 3: elitism = False, truncation = True, 
    #               duplicates = False
    deterministic = FunctionData('deterministic_sampling', 
                                 truncation=3)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = list(pop1.select())
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) == 1
    # Check for truncation.   
    assert len(pop2) == len(pop1) - 3   

   # Conditions 4: elitism = False, truncation = True, 
    #               duplicates = True
    deterministic = FunctionData('deterministic_sampling', 
                                 truncation=3, duplicates=True)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = list(pop1.select())
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) > 1
    # Check for truncation.   
    assert min(x.fitness for x in pop1) < min(x.fitness for x in pop2)    
    
   # Conditions 5: elitism = True, truncation = False, 
    #               duplicates = False
    deterministic = FunctionData('deterministic_sampling', 
                                 elitism=3)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = []
    gen = pop1.select()
    for x in range(4):
        pop2.append(next(gen))
    
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) == 1
    
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    for x in fittest:
        assert x in pop2

   # Conditions 6: elitism = True, truncation = False, 
    #               duplicates = True
    deterministic = FunctionData('deterministic_sampling', 
                                 elitism=3,duplicates=True)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = []
    gen = pop1.select()
    for x in range(4):
        pop2.append(next(gen))
    
    for i, ind in enumerate(pop2):
        if i > 0 and i < 3:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) == 2
    
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    for x in fittest:
        assert x in pop2

   # Conditions 7: elitism = True, truncation = True, 
    #               duplicates = False
    deterministic = FunctionData('deterministic_sampling', 
                                 truncation=3, elitism=3, duplicates=False)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = list(pop1.select())
    
    for i, ind in enumerate(pop2):
        if i > 0:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) == 1
    
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    for x in fittest:
        assert x in pop2
    assert len(pop2) == len(pop1) -3

   # Conditions 8: elitism = True, truncation = True, 
    #               duplicates = True
    deterministic = FunctionData('deterministic_sampling', 
                                 elitism=3, duplicates=True, 
                                 truncation=3)
    selector = Selection(deterministic, 'a', 'b')
    ga_tools = GATools(selector, 'a', 'b', 'do_not_optimize', 'cage')
    pop1.ga_tools = ga_tools   
        
    # Make one of the fitness higher than average         
    pop1[0].fitness = 100
    pop2 = list(pop1.select())
    for i, ind in enumerate(pop2):
        if i > 0 and i != 3:
            assert ind.fitness <= pop2[i-1].fitness
    count = Counter(pop2)
    assert max(count.values()) > 1
    
    fittest = sorted(pop1, reverse=True, key=attrgetter('fitness'))[:3]
    for x in fittest:
        assert x in pop2












