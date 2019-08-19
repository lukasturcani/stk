"""
Tests for normalization functions.

"""

import stk


def test_power(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Power(2)
    normalizer._handle_failed = False

    tmp_amine2.fitness = 1
    tmp_aldehyde2.fitness = 2
    tmp_bromine2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 1
    assert tmp_aldehyde2.fitness == 4
    assert tmp_bromine2.fitness == 9

    tmp_amine2.fitness = [1, 2, 3]
    tmp_aldehyde2.fitness = [4, 5, 6]
    tmp_bromine2.fitness = [7, 8, 9]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [1, 4, 9]
    assert tmp_aldehyde2.fitness.tolist() == [16, 25, 36]
    assert tmp_bromine2.fitness.tolist() == [49, 64, 81]


def test_multiply(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Multiply(2)
    normalizer._handle_failed = False

    tmp_amine2.fitness = 1
    tmp_aldehyde2.fitness = 2
    tmp_bromine2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 2
    assert tmp_aldehyde2.fitness == 4
    assert tmp_bromine2.fitness == 6

    tmp_amine2.fitness = [1, 2, 3]
    tmp_aldehyde2.fitness = [4, 5, 6]
    tmp_bromine2.fitness = [7, 8, 9]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [2, 4, 6]
    assert tmp_aldehyde2.fitness.tolist() == [8, 10, 12]
    assert tmp_bromine2.fitness.tolist() == [14, 16, 18]


def test_sum(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.Sum()
    normalizer._handle_failed = False

    tmp_amine2.fitness = [1, 2, 3]
    tmp_aldehyde2.fitness = [4, 5, 6]
    tmp_bromine2.fitness = [7, 8, 9]

    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 6
    assert tmp_aldehyde2.fitness == 15
    assert tmp_bromine2.fitness == 24


def test_scale_by_mean(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.ScaleByMean()
    normalizer._handle_failed = False

    tmp_amine2.fitness = 1
    tmp_aldehyde2.fitness = 2
    tmp_bromine2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 0.5
    assert tmp_aldehyde2.fitness == 1
    assert tmp_bromine2.fitness == 1.5

    tmp_amine2.fitness = [1, 10, 100]
    tmp_aldehyde2.fitness = [2, 20, 200]
    tmp_bromine2.fitness = [3, 30, 300]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [0.5, 0.5, 0.5]
    assert tmp_aldehyde2.fitness.tolist() == [1, 1, 1]
    assert tmp_bromine2.fitness.tolist() == [1.5, 1.5, 1.5]


def test_shift_up(tmp_amine2, tmp_aldehyde2, tmp_bromine2):
    normalizer = stk.ShiftUp()
    normalizer._handle_failed = False

    tmp_amine2.fitness = 1
    tmp_aldehyde2.fitness = -2
    tmp_bromine2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 4
    assert tmp_aldehyde2.fitness == 1
    assert tmp_bromine2.fitness == 6

    tmp_amine2.fitness = [1, -5, 5]
    tmp_aldehyde2.fitness = [3, -10, 2]
    tmp_bromine2.fitness = [2, 20, 1]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [1, 6, 5]
    assert tmp_aldehyde2.fitness.tolist() == [3, 1, 2]
    assert tmp_bromine2.fitness.tolist() == [2, 31, 1]


def test_normalizer_sequence(
    tmp_amine2,
    tmp_aldehyde2,
    tmp_bromine2
):
    sequence = stk.NormalizerSequence(
        stk.Power(2),
        stk.Sum()
    )
    for norm in sequence._normalizers:
        norm._handle_failed = False

    tmp_amine2.fitness = [1, 2, 3]
    tmp_aldehyde2.fitness = [4, 5, 6]
    tmp_bromine2.fitness = [7, 8, 9]
    pop = stk.Population(tmp_amine2, tmp_aldehyde2, tmp_bromine2)

    sequence.normalize(pop)

    assert tmp_amine2.fitness == 14
    assert tmp_aldehyde2.fitness == 77
    assert tmp_bromine2.fitness == 194
