import stk


def test_bb_lk_exchange(pop):
    m1, m2 = pop[:2]
    parent_tops = {m1.topology.__class__, m2.topology.__class__}
    offspring_pop = stk.Crossover.bb_lk_exchange(None, m1, m2)
    for offspring in offspring_pop:

        assert offspring.topology.__class__ in parent_tops
        bb1, bb2 = offspring.building_blocks

        if bb1 in m1.building_blocks:
            assert bb2 not in m1.building_blocks
            assert bb2 in m2.building_blocks
        else:
            assert bb2 in m1.building_blocks
            assert bb2 not in m2.building_blocks
