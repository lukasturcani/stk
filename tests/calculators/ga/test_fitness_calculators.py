import stk
import pytest


class _Thing:
    ...


def test_property_vector(tmp_amine2):
    fitness_calculator = stk.PropertyVector(
        lambda mol: 1,
        lambda mol: 2,
        lambda mol: 3,
        lambda mol: 4
    )
    assert fitness_calculator.get_fitness(tmp_amine2) == [1, 2, 3, 4]


def test_is_caching():
    assert stk.PropertyVector(use_cache=True).is_caching()
    assert not stk.PropertyVector(use_cache=False).is_caching()


def test_add_to_cache(tmp_amine2):
    calculator = stk.PropertyVector(use_cache=True)
    fitness = object()
    calculator.add_to_cache(tmp_amine2, fitness)
    assert calculator.get_fitness(tmp_amine2) is fitness


def test_cache_use(tmp_amine2):
    calls = 0
    things = [_Thing(), _Thing(), _Thing(), _Thing()]

    def prop(mol):
        nonlocal calls
        nonlocal things
        calls += 1
        r = things[0]
        things.remove(r)
        return r

    calc = stk.PropertyVector(
        prop,
        use_cache=False
    )
    fst = calc.get_fitness(tmp_amine2)
    assert calls == 1
    assert calc.get_fitness(tmp_amine2) is not fst
    assert calls == 2

    # Test that the cache is being filled when use_cache is True.
    calc = stk.PropertyVector(
        prop,
        use_cache=True
    )
    snd = calc.get_fitness(tmp_amine2)
    assert calls == 3
    assert calc.get_fitness(tmp_amine2) is snd
    assert calls == 3


def test_attribute_creation(tmp_amine2):
    calc = stk.PropertyVector(lambda mol: 1)
    assert not hasattr(tmp_amine2, 'fitness')
    calc.get_fitness(tmp_amine2)
    assert tmp_amine2.fitness == [1]


def test_raising_fitness_calculator(tmp_amine2):
    fitness_calculator = stk.PropertyVector(lambda m: 1)
    never_raiser = stk.RaisingFitnessCalculator(fitness_calculator, 0)
    never_raiser.get_fitness(tmp_amine2)

    always_raiser = stk.RaisingFitnessCalculator(fitness_calculator, 1)
    with pytest.raises(stk.RaisingFitnessCalculatorError):
        always_raiser.get_fitness(tmp_amine2)
