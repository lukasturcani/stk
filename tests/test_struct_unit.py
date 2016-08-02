import os
import rdkit

from ..classes import StructUnit, FGInfo


def get_mol_file():
    # The following lines first create a directory tree starting from
    # the current working directory. A generator expression then
    # iterates through the directory, where the ``if`` condition ensures
    # that the desired ``.mol`` file is found. It is the one in the 
    # ``data`` directory. Finally the full path of the ``.mol`` file is 
    # generated using the ``os.path.join`` function. This approach means
    # that the test should work on any machine as it does not depend on
    # absolute paths to find the ``.mol`` file. It also means that the 
    # test does not need to be run from a specific directory.    
    mol_file = os.walk(os.getcwd())    
    for x in mol_file:
        if 'data' in x[0]:
            for y in x[2]:
                if '.mol' in y and 'HEAVY' not in y:
                    yield os.path.join(x[0], y)



def test_init():
    """
    Ensures that StructUnit instances are initiated correctly.
    
    This function only checks that the attributes exist and that
    their values of of the correct type. Checking that correct values
    are initialized is done in other tests.

    """     
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
     
    # Check that heavy attributes were created by the initializer.
    assert hasattr(struct_unit, 'heavy_mol')
    assert hasattr(struct_unit, 'heavy_mol_file')
    assert hasattr(struct_unit, 'heavy_smiles')
    
    assert isinstance(struct_unit.func_grp, FGInfo)
    assert isinstance(struct_unit.prist_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.heavy_mol, rdkit.Chem.rdchem.Mol)
    assert isinstance(struct_unit.prist_smiles, str)
    assert isinstance(struct_unit.heavy_smiles, str)
    assert isinstance(struct_unit.prist_mol_file, str)
    assert isinstance(struct_unit.heavy_mol_file, str)
    
def test_find_functional_group_atoms():
    """
    Make sure correct atoms are found in the functional groups.
    
    This function uses the ``aldehyde2f_3.mol`` file.

    """    
    # These are the expected atom ids for the test molecule.
    expected = ((1, 0, 12), (10, 11, 19))    
    
    # Initializing the test molecule.    
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
        
    func_grp_atoms = struct_unit.find_functional_group_atoms()

    assert func_grp_atoms == expected
    
def test_shift_heavy_mol():
    """
    Ensure that shifting the position of a ``StructUnit`` works.    
    
    This function uses the ``aldehyde2f_3.mol`` file.    
    
    """
    # Initializing the test molecule.    
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)

    # Shifting the same molecule twice should return two 
    # ``rdkit.Chem.rdchem.Mol`` instances with conformers describing the
    # same atomic positions. Furthermore, the original conformer should
    # in the ``StructUnit`` instance should be unchanged.
    og_conformer = struct_unit.heavy_mol.GetConformer()

    shift_size = 10    
    a = struct_unit.shift_heavy_mol(shift_size, shift_size, shift_size)
    a_conformer = a.GetConformer()
    
    b = struct_unit.shift_heavy_mol(shift_size, shift_size, shift_size)
    b_conformer = b.GetConformer()
    
    # Check that the same atom coords are present in `a` and `b`. Also
    # Check that these are different to the coords in the original.
    for atom in a.GetAtoms():
        atom_id = atom.GetIdx()
      
        og_coords = og_conformer.GetAtomPosition(atom_id)
        atom1_coords = a_conformer.GetAtomPosition(atom_id)
        atom2_coords = b_conformer.GetAtomPosition(atom_id)
        
        assert atom1_coords.x == atom2_coords.x
        assert atom1_coords.y == atom2_coords.y
        assert atom1_coords.z == atom2_coords.z
        
        # Checking that the coords are differnt to original. By the 
        # correct amount.        
        assert atom1_coords.x == og_coords.x + shift_size
        assert atom1_coords.y == og_coords.y + shift_size
        assert atom1_coords.z == og_coords.z + shift_size
        
def test_get_heavy_coords():
    """
    Make sure the correct output is provided.

    This function uses the ``aldehyde2f_3.mol`` file.
    
    """
    
    exp_output = [(-1.4115803404500527, 1.8779299561981242, 
                  1.1334757839928156), (-0.8352914448038408, 
                  0.921636471037717, 0.4640498148275688), 
                  (-0.9934336108721664, -0.46310619741690934, 
                  0.9297809098773657), (0.009983650030028485, 
                  -1.374989887966107, 0.5029433827643824), 
                  (0.5857664717779867, -1.1253267727860947, 
                  -0.7814930422117815), (1.0940416957920942, 
                  0.27300219443714546, -0.8475140299159067), 
                  (0.1636399430227658, 1.3016922000662716, 
                  -0.4879216899401005), (2.183353421827768, 
                  0.6152269876890724, -1.4230929234365053), 
                  (-0.9792559513888446, -0.43625964749933865, 
                  2.045856432494676), (-1.9980521717049462, 
                  -0.8115881557059094, 0.5844968489801275), 
                  (0.7184667366156154, -1.6378536276273608, 
                  1.2121664696477477), (-0.0922488825818627,
                  -1.345866937866273, -1.623993788008611), 
                  (1.4793159165999186, -1.7897896888108131, 
                  -0.8712425966393919), (0.06887317227154656, 
                  2.1882763097889426, -1.0468793931272447), 
                  (-3.3836844645263615, 1.0794453705615819, 
                  -0.30623877663805055), (3.3901058583903567, 
                  0.7275714258999577, 0.5156065973329162)]

    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
    output = list(struct_unit.get_heavy_coords())
    assert exp_output == output
    
def test_amine_substitution():
    """
    Ensure that the amine functional group is correctly replaced.    
    
    """
    exp_smiles = ("[H][C]([H])([H])[C]1=[C]([Rh])[C](=[O])[C]2=[C]"
                  "([C]1=[O])[N]1[C](=[C]2[C]([H])([H])[O][C](=[O])"
                  "[Rh])[C]([H])([H])[C]([H])([Rh])[C]1([H])[H]")
    mol_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)   
    mol = StructUnit(mol_file)
    assert mol.heavy_smiles == exp_smiles

def test_aldehyde_substitution():
    """
    Ensure that the aldehyde functional group is correctly replaced.    
    
    """
    exp_smiles = ("[H][N]1[C](=[N][Y])[C]([H])([H])[N]([H])[C]([H])"
                  "([H])[C]1=[N][Y]")
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)   
    mol = StructUnit(mol_file)
    assert mol.heavy_smiles == exp_smiles
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    