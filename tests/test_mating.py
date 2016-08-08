from collections import Counter


from .test_population import generate_population
from .test_struct_unit import get_mol_file
from ..classes import GATools, Selection, FunctionData, BuildingBlock, Mating
from ..classes import  Linker, Cage, FourPlusSix, EightPlusTwelve, Population
from ..convenience_functions import flatten

bb_file = next(x for x in get_mol_file() 
                                    if 'amine3f_14.mol' in x)
lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

bb = BuildingBlock(bb_file)
lk = Linker(lk_file)    
building_blocks = (bb, lk)
mol = Cage(building_blocks, FourPlusSix, 'you_can_delete_this3.mol')

bb2_file = next(x for x in get_mol_file() 
                                    if 'amine3f_5.mol' in x)
lk2_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_28.mol' in x) 

bb2 = BuildingBlock(bb2_file)
lk2 = Linker(lk2_file)    
building_blocks2 = (bb2, lk2)
mol2 = Cage(building_blocks2, EightPlusTwelve, 'you_can_delete_this4.mol')

def test_bb_lk_exchange():
    
    
        # Make a ``GATools`` attribute and give it to the population.
        all_combs = FunctionData('all_combinations')
        bb_lk_exchange = FunctionData('bb_lk_exchange')

        
        selector = Selection('a', all_combs, 'a')
        mator = Mating(bb_lk_exchange, 1)
        ga_tools = GATools(selector, mator, 'b')


        pop1 = Population(ga_tools, mol, mol2)
        pop2 = pop1.gen_offspring()

        for mol3 in pop2:
            assert mol3 not in pop1
    
        all_topologies = {type(x.topology) for x in pop2}
        assert EightPlusTwelve in all_topologies and FourPlusSix in all_topologies
        
        all_building_blocks = list(flatten([x.building_blocks for x in pop2]))
        all_bb_count = Counter(all_building_blocks)
        assert all(val == 2 for val in all_bb_count.values())
    