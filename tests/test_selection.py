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
    fittest = FunctionData('fittest', size=3)
    fittest2 = FunctionData('fittest', size=4)
    fittest3 = FunctionData('fittest', size=5)
    
    selector = Selection(fittest, fittest2, fittest3)
    ga_tools = GATools(selector, 'a', 'b')
    pop1.ga_tools = ga_tools
    
    pop2 = pop1.select('generational')
    pop3 = pop1.select('mating')
    pop4 = pop1.select('mutation')

    assert len(pop2) == 3
    assert len(pop3) == 4
    assert len(pop4) == 5
    
    assert all(mol in pop3 and mol in pop4 for mol in pop2)
    assert all(mol in pop4 for mol in pop3)
    
def test_fittest():
    fittest = FunctionData('fittest', size=3)
    selector = Selection(fittest, fittest, fittest)
    ga_tools = GATools(selector, 'a', 'b')
    pop1.ga_tools = ga_tools

    #Set fitnesses and labels
    for index, mol in enumerate(pop1):
        mol.fitness = 1/(index+1)
        mol.label = labels[index]
    
    pop2 = pop1.select()

    
    assert len(pop2) == 3
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
    assert not pop2.members
    for parents in pop2.populations:
        assert len(parents.members) == 2
        assert not parents.populations




    