from ..classes import Polymer, Linker, BlockCopolymer
from .test_struct_unit import get_mol_file

def test_init():
    lk1_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    lk2_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_12.mol' in x)
    lk3_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_74.mol' in x)
    lk4_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_96.mol' in x)
                                        
    lk1 = Linker(lk1_file)
    lk2 = Linker(lk2_file)
    lk3 = Linker(lk3_file)
    lk4 = Linker(lk4_file)
    
    building_blocks = (lk1, lk2, lk3, lk4)    
    
    polymer = Polymer(building_blocks, BlockCopolymer, 'you_can_delete_this2.mol', 
                      topology_args= ["AAABBACDABCDD"])
    