"""
Tests for normalization functions.

"""

import stk


def test_power(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Power(2)

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: 1,
        tmp_aldehyde2: 2,
        tmp_bromine2: 3,
    })
    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == 1
    assert normalized[tmp_aldehyde2] == 4
    assert normalized[tmp_bromine2] == 9

    pop.set_fitness_values_from_dict({
        tmp_amine2: [1, 2, 3],
        tmp_aldehyde2: [4, 5, 6],
        tmp_bromine2: [7, 8, 9],
    })

    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2].tolist() == [1, 4, 9]
    assert normalized[tmp_aldehyde2].tolist() == [16, 25, 36]
    assert normalized[tmp_bromine2].tolist() == [49, 64, 81]


def test_multiply(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Multiply(2)

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: 1,
        tmp_aldehyde2: 2,
        tmp_bromine2: 3,
    })
    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == 2
    assert normalized[tmp_aldehyde2] == 4
    assert normalized[tmp_bromine2] == 6

    pop.set_fitness_values_from_dict({
        tmp_amine2: [1, 2, 3],
        tmp_aldehyde2: [4, 5, 6],
        tmp_bromine2: [7, 8, 9],
    })

    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == [2, 4, 6]
    assert normalized[tmp_aldehyde2] == [8, 10, 12]
    assert normalized[tmp_bromine2] == [14, 16, 18]


def test_sum(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Sum()

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: [1, 2, 3],
        tmp_aldehyde2: [4, 5, 6],
        tmp_bromine2: [7, 8, 9],
    })
    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == 6
    assert normalized[tmp_aldehyde2] == 15
    assert normalized[tmp_bromine2] == 24


def test_divide_by_mean(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.DivideByMean()

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: 1,
        tmp_aldehyde2: 2,
        tmp_bromine2: 3,
    })
    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == 0.5
    assert normalized[tmp_aldehyde2] == 1
    assert normalized[tmp_bromine2] == 1.5

    stk.EAPopulation({
        tmp_amine2: [1, 10, 100],
        tmp_aldehyde2: [2, 20, 200],
        tmp_bromine2: [3, 30, 300],
    })

    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2].tolist() == [0.5, 0.5, 0.5]
    assert normalized[tmp_aldehyde2].tolist() == [1, 1, 1]
    assert normalized[tmp_bromine2].tolist() == [1.5, 1.5, 1.5]


def test_shift_up(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.ShiftUp()

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: 1,
        tmp_aldehyde2: -2,
        tmp_bromine2: 3,
    })
    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2] == 4
    assert normalized[tmp_aldehyde2] == 1
    assert normalized[tmp_bromine2] == 6

    pop.set_fitness_values_from_dict({
        tmp_amine2: [1, -5, 5],
        tmp_aldehyde2: [3, -10, 2],
        tmp_bromine2: [2, 20, 1],
    })

    normalized = normalizer.normalize(pop)

    assert normalized[tmp_amine2].tolist() == [1, 6, 5]
    assert normalized[tmp_aldehyde2].tolist() == [3, 1, 2]
    assert normalized[tmp_bromine2].tolist() == [2, 31, 1]


def test_normalizer_sequence(
    tmp_amine2,
    tmp_aldehyde2,
    tmp_bromine2,
):
    sequence = stk.Sequence(
        stk.Power(2),
        stk.Sum()
    )

    pop = stk.EAPopulation(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    pop.set_fitness_values_from_dict({
        tmp_amine2: [1, 2, 3],
        tmp_aldehyde2: [4, 5, 6],
        tmp_bromine2: [7, 8, 9],
    })
    normalized = sequence.normalize(pop)

    assert normalized[tmp_amine2] == 14
    assert normalized[tmp_aldehyde2] == 77
    assert normalized[tmp_bromine2] == 194
