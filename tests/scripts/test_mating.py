from collections import Counter
import os

from ...classes import FunctionData, Mating
from ...classes import  FourPlusSix, EightPlusTwelve, Population
from ...convenience_functions import flatten

pop_file = os.path.join('data', 'mating', 'test_pop')
pop1 = Population.load(pop_file)

def test_bb_lk_exchange():
        
        bb_lk_exchange = FunctionData('bb_lk_exchange')
        mator = Mating(bb_lk_exchange, 1)
        pop1.ga_tools.mating = mator

        
        pop2 = pop1.gen_offspring()

        for mol3 in pop2:
            assert mol3 not in pop1

        all_topologies = {type(x.topology) for x in pop2}
        assert EightPlusTwelve in all_topologies and FourPlusSix in all_topologies
        
        all_building_blocks = list(flatten([x.building_blocks for x in pop2]))
        all_bb_count = Counter(all_building_blocks)
        assert all(val == 2 for val in all_bb_count.values())
    