from ..ga import GAInput
from ..convenience_tools import FunctionData
from ..molecular import FourPlusSix
from os.path import join

i = GAInput(join('data', 'gainput', 'test.py'))


def test_init():

    assert i.processes == 1

    assert i.init_func == {'NAME': 'init_random_cages',
                           'bb_db': 'path1',
                           'lk_db': 'path2',
                           'topologies': [FourPlusSix()]}

    assert i.generational_select_func == {
                'NAME': 'stochastic_sampling',
                'use_rank': True}

    assert i.crossover_select_func == {'NAME': 'crossover_roulette'}

    assert i.mutation_select_func == {'NAME': 'stochastic_sampling',
                                      'duplicates': True}

    assert i.crossover_funcs == [{'NAME': 'bb_lk_exchange'}]

    assert i.mutation_funcs == [
        {'NAME': 'cage_random_bb',
         'database': 'path1'},
        {'NAME': 'cage_random_lk',
         'database': 'path2'}]

    assert i.mutation_weights is None
    assert i.opt_func == {'NAME': 'do_not_optimize'}
    assert i.fitness_func == {'NAME': 'random_fitness_vector'}

    assert i.normalization_funcs == [
        {'NAME': 'shift_elements',
         'indices': [-1]},
        {'NAME': 'magnitudes'}]

    assert i.num_generations == 5
    assert i.num_mutations == 2
    assert i.num_crossovers == 2
    assert i.pop_size == 5
    assert i.comparison_pops == ['pop1', 'pop2']
    assert i.exit_func == {'NAME': 'no_exit'}
    assert i.databases == ['1', '2']


def test_crosser():
    crosser = i.crosser()
    assert crosser.funcs == [FunctionData('bb_lk_exchange')]
    assert crosser.num_crossovers == i.num_crossovers
    assert crosser.weights is None


def test_exiter():
    exiter = i.exiter()
    assert exiter.func_data == FunctionData('no_exit')


def test_fitnessor():
    fitnessor = i.fitnessor()
    assert fitnessor.name == 'random_fitness_vector'


def test_selector():
    sel = i.selector()
    assert sel.generational == FunctionData('stochastic_sampling',
                                            use_rank=True)
    assert sel.crossover == FunctionData('crossover_roulette')
    assert sel.mutation == FunctionData('stochastic_sampling',
                                        duplicates=True)


def test_mutator():
    mut = i.mutator()
    assert mut.funcs == [FunctionData('cage_random_bb',
                                      database='path1'),
                         FunctionData('cage_random_lk',
                                      database='path2')]
    assert mut.weights is None
    assert mut.num_mutations == 2


def test_normalizer():
    norm = i.normalizer()
    assert norm.funcs == [FunctionData('shift_elements',
                                       indices=[-1]),
                          FunctionData('magnitudes')]


def test_opter():
    opter = i.opter()
    assert opter.name == 'do_not_optimize'


def test_ga_tools():
    gatools = i.ga_tools()
    assert hasattr(gatools, 'selection')
    assert hasattr(gatools, 'crossover')
    assert hasattr(gatools, 'mutation')
    assert hasattr(gatools, 'normalization')
    assert hasattr(gatools, 'optimization')
    assert hasattr(gatools, 'fitness')
    assert hasattr(gatools, 'input')
    assert hasattr(gatools, 'exit')
