import stk


def test_init(ga_input):

    assert ga_input.processes == 1

    assert ga_input.init_func == {
                           'NAME': 'init_random_cages',
                           'bb_db': 'path1',
                           'lk_db': 'path2',
                           'topologies': [stk.FourPlusSix()]}

    assert ga_input.generational_select_func == {
                'NAME': 'stochastic_sampling',
                'use_rank': True}

    assert ga_input.crossover_select_func == {
                                        'NAME': 'crossover_roulette'}

    assert ga_input.mutation_select_func == {
                                      'NAME': 'stochastic_sampling',
                                      'duplicates': True}

    assert ga_input.crossover_funcs == [{'NAME': 'bb_lk_exchange'}]

    assert ga_input.mutation_funcs == [
        {'NAME': 'cage_random_bb',
         'database': 'path1'},
        {'NAME': 'cage_random_lk',
         'database': 'path2'}]

    assert ga_input.mutation_weights is None
    assert ga_input.opt_func == {'NAME': 'do_not_optimize'}
    assert ga_input.fitness_func == {'NAME': 'random_fitness_vector'}

    assert ga_input.normalization_funcs == [
        {'NAME': 'shift_elements',
         'indices': [-1]},
        {'NAME': 'magnitudes'}]

    assert ga_input.num_generations == 5
    assert ga_input.num_mutations == 2
    assert ga_input.num_crossovers == 2
    assert ga_input.pop_size == 5
    assert ga_input.comparison_pops == ['pop1', 'pop2']
    assert ga_input.exit_func == {'NAME': 'no_exit'}
    assert ga_input.databases == ['1', '2']


def test_crosser(ga_input):
    crosser = ga_input.crosser()
    assert crosser.funcs == [stk.FunctionData('bb_lk_exchange')]
    assert crosser.num_crossovers == ga_input.num_crossovers
    assert crosser.weights is None


def test_exiter(ga_input):
    exiter = ga_input.exiter()
    assert exiter.func_data == stk.FunctionData('no_exit')


def test_fitnessor(ga_input):
    fitnessor = ga_input.fitnessor()
    assert fitnessor.name == 'random_fitness_vector'


def test_selector(ga_input):
    sel = ga_input.selector()
    assert sel.generational == stk.FunctionData(
                                            'stochastic_sampling',
                                            use_rank=True)
    assert sel.crossover == stk.FunctionData('crossover_roulette')
    assert sel.mutation == stk.FunctionData(
                                        'stochastic_sampling',
                                        duplicates=True)


def test_mutator(ga_input):
    mut = ga_input.mutator()
    assert mut.funcs == [stk.FunctionData(
                                      'cage_random_bb',
                                      database='path1'),
                         stk.FunctionData(
                                      'cage_random_lk',
                                      database='path2')]
    assert mut.weights is None
    assert mut.num_mutations == 2


def test_normalizer(ga_input):
    norm = ga_input.normalizer()
    assert norm.funcs == [stk.FunctionData(
                                       'shift_elements',
                                       indices=[-1]),
                          stk.FunctionData('magnitudes')]


def test_opter(ga_input):
    opter = ga_input.opter()
    assert opter.name == 'do_not_optimize'


def test_ga_tools(ga_input):
    gatools = ga_input.ga_tools()
    assert hasattr(gatools, 'selection')
    assert hasattr(gatools, 'crossover')
    assert hasattr(gatools, 'mutation')
    assert hasattr(gatools, 'normalization')
    assert hasattr(gatools, 'fitness')
    assert hasattr(gatools, 'input')
    assert hasattr(gatools, 'exit')
