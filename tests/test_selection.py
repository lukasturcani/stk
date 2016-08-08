from operator import attrgetter

from .test_population import generate_population
from ..classes import GATools, Selection, FunctionData

labels = "abcdefghijklmnopqrstuvwxyz"
# Make a population.
pop1 = generate_population()


        
def test_selection_population_integration():

    #Set fitnesses and labels

    for index, mol in enumerate(pop1):
        mol.fitness = index
        mol.label = labels[index]

    # Make a ``GATools`` attribute and give it to the population.
    fittest = FunctionData('fittest')
    fittest2 = FunctionData('fittest')
    fittest3 = FunctionData('fittest')
    
    selector = Selection(fittest, fittest2, fittest3)
    ga_tools = GATools(selector, 'a', 'b')
    pop1.ga_tools = ga_tools
    
    pop2 = pop1.select('generational')
    pop3 = pop1.select('mating')
    pop4 = pop1.select('mutation')

    pop1_members = list(pop1)
    pop1_members.sort(key=attrgetter('fitness'), reverse=True)
    assert list(pop2) == pop1_members
    assert list(pop3) ==  pop1_members
    assert list(pop4) == pop1_members
    
    
def test_fittest():
    fittest = FunctionData('fittest')
    selector = Selection(fittest, fittest, fittest)
    ga_tools = GATools(selector, 'a', 'b')
    pop1.ga_tools = ga_tools

    #Set fitnesses and labels
    for index, mol in enumerate(pop1):
        mol.fitness = 1/(index+1)
        mol.label = labels[index]
    
    pop2 = pop1.select()

    
    assert len(list(pop2)) == 22
    for mol in pop2:
        assert mol.label in "abc"
        
def test_all_combinations():
    all_combs = FunctionData('all_combinations')
    selector = Selection('a', all_combs, 'b')
    ga_tools = GATools(selector, 'a', 'b')
    pop1.ga_tools = ga_tools
    
    #Set fitnesses and labels
    for index, mol in enumerate(pop1):
        mol.fitness = 1/(index+1)
        mol.label = labels[index]
    
    pop2 = pop1.select('mating')
    assert list(pop2)
    for parents in pop2:
        assert len(parents) == 2
    




    