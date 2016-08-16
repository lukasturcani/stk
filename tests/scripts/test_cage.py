from .test_struct_unit import get_mol_file
from ...classes import Cage, FourPlusSix, BuildingBlock, Linker



def test_init():
    """
    Ensure that cages initialize attributes successfully.
    
    This function uses the ``aldehyde2f_3.mol`` and
    ``amine3f_14.mol`` files.
    
    Only presence of the attributes is tested here. Testing whether the
    correct cage was built should be done in the topology testing
    module.

    """ 
    
    bb_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)
    lk_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)    
    
    bb = BuildingBlock(bb_file)
    lk = Linker(lk_file)    
    building_blocks = (bb, lk)
    cage = Cage(building_blocks, FourPlusSix, 
                                            'you_can_delete_this.mol')
    
    assert hasattr(cage, 'prist_mol_file')
    assert hasattr(cage, 'heavy_mol_file')
    assert hasattr(cage, 'prist_mol')
    assert hasattr(cage, 'heavy_mol')
    assert hasattr(cage, 'prist_smiles')
    assert hasattr(cage, 'heavy_smiles')
    assert hasattr(cage, 'topology')
    assert hasattr(cage, 'optimized')