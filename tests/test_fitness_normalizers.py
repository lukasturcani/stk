"""
Tests for normalization functions.

"""

import stk


def test_power(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2):
    normalizer = stk.Power(2)

    tmp_amine2.fitness = 1
    tmp_amine2_alt1.fitness = 2
    tmp_amine2_alt2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 1
    assert tmp_amine2_alt1.fitness == 4
    assert tmp_amine2_alt2.fitness == 9

    tmp_amine2.fitness = [1, 2, 3]
    tmp_amine2_alt1.fitness = [4, 5, 6]
    tmp_amine2_alt2.fitness = [7, 8, 9]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [2, 4, 9]
    assert tmp_amine2_alt1.fitness.tolist() == [16, 25, 36]
    assert tmp_amine2_alt2.fitness.tolist() == [49, 64, 81]


def test_multiply(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2):
    normalizer = stk.Multiply(2)

    tmp_amine2.fitness = 1
    tmp_amine2_alt1.fitness = 2
    tmp_amine2_alt2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 2
    assert tmp_amine2_alt1.fitness == 4
    assert tmp_amine2_alt2.fitness == 6

    tmp_amine2.fitness = [1, 2, 3]
    tmp_amine2_alt1.fitness = [4, 5, 6]
    tmp_amine2_alt2.fitness = [7, 8, 9]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [2, 4, 6]
    assert tmp_amine2_alt1.fitness.tolist() == [8, 10, 12]
    assert tmp_amine2_alt2.fitness.tolist() == [14, 16, 18]


def test_sum(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2):
    normalizer = stk.Sum()

    tmp_amine2.fitness = [1, 2, 3]
    tmp_amine2_alt1.fitness = [4, 5, 6]
    tmp_amine2_alt2.fitness = [7, 8, 9]

    pop = stk.Population(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 6
    assert tmp_amine2_alt1.fitness == 15
    assert tmp_amine2_alt2.fitness == 24


def test_scale_by_mean(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2):
    normalizer = stk.ScaleByMean()

    tmp_amine2.fitness = 1
    tmp_amine2_alt1.fitness = 2
    tmp_amine2_alt2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 0.5
    assert tmp_amine2_alt1.fitness == 1
    assert tmp_amine2_alt2.fitness == 1.5

    tmp_amine2.fitness = [1, 2, 3]
    tmp_amine2_alt1.fitness = [10, 20, 30]
    tmp_amine2_alt2.fitness = [100, 200, 300]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [0.5, 0.5, 0.5]
    assert tmp_amine2_alt1.fitness.tolist() == [1, 1, 1]
    assert tmp_amine2_alt2.fitness.tolist() == [1.5, 1.5, 1.5]


def test_shift_up(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2):
    normalizer = stk.ShiftUp()

    tmp_amine2.fitness = 1
    tmp_amine2_alt1.fitness = -2
    tmp_amine2_alt2.fitness = 3

    pop = stk.Population(tmp_amine2, tmp_amine2_alt1, tmp_amine2_alt2)
    normalizer.normalize(pop)

    assert tmp_amine2.fitness == 4
    assert tmp_amine2_alt1.fitness == 1
    assert tmp_amine2_alt2.fitness == 6

    tmp_amine2.fitness = [1, -5, 5]
    tmp_amine2_alt1.fitness = [3, -10, 2]
    tmp_amine2_alt2.fitness = [2, 20, 1]

    normalizer.normalize(pop)

    assert tmp_amine2.fitness.tolist() == [1, 6, 5]
    assert tmp_amine2_alt1.fitness.tolist() == [3, 1, 2]
    assert tmp_amine2_alt2.fitness.tolist() == [2, 31, 1]
