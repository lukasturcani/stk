from ..ga import GAInput
from ..convenience_tools import FunctionData
from ..molecular import FourPlusSix
from os.path import join

def test_init():
    i = GAInput(join('data', 'gainput', 'test.in'))
    assert i.init_func == FunctionData('init_random_cages',
                                        bb_db='path1', lk_db='path2')
    assert i.generational_select_func == FunctionaData(
    'stochastic_sampling', use_rank=True)
    assert i.parent_select_func == FunctionData('crossover_roulette',
    duplicates=True)
    assert i.mutant_select_func == FunctionData('stochastic_sampling',
    duplicates=True)
    assert i.crossover_func = FunctionData('bb_lk_exchange')
    assert i.mutation_func == [FunctionData('cage_random_bb',
    database='path1'),
    FunctionData('cage_random_lk', database='path2')]
    assert i.mutation_weights == [1/4,1/4,1/4,1/4]
    assert i.opt_func == FunctionData('do_not_optimize')
    assert i.fitness_func == FunctionData('random_fitness_vector')
    assert i.normalization_func == [FunctionData('shift_elements',
    indices=[-1]), FunctionData('magnitudes')]
    assert i.num_generations == 5
    assert i.num_mutations == 2
    assert i.num_crossovers == 2
    assert i.pop_size == 5
    assert i.comparison_pops == ['pop1', 'pop2']
