from ..classes import StructUnit
import os

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
    
    This function uses the ``aldehyde2f_3.mol`` file.

    """     
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
     
    # Check that heavy attributes were created by the initializer.
    assert hasattr(struct_unit, 'heavy_mol')
    assert hasattr(struct_unit, 'heavy_mol_file')
    assert hasattr(struct_unit, 'heavy_smiles')
    
    # Check that the contents of the `heavy_mol_file` and 
    # `heavy_mol_smiles` are correct.
    
    # First create the path where the expected content of 
    # `heavy_mol_file` is located.
    expected_mol_file_dir = os.path.dirname(mol_file)
    expected_mol_file = os.path.join(expected_mol_file_dir, 
                                     'expected_mol_file')
    
    # Load the data and compare.
    with open(expected_mol_file, 'r') as exp_file:
        exp_content = exp_file.read()
    with open(struct_unit.heavy_mol_file, 'r') as actual_file:
        actual_content = actual_file.read()
    
    assert exp_content == actual_content

    # Check that  generated and expected SMILES are the same.
    expected_smiles_file = os.path.join(expected_mol_file_dir, 
                                   'expected_smiles')
    
    with open(expected_smiles_file, 'r') as exp_smiles_file:
        exp_smiles = exp_smiles_file.read()
        
    assert exp_smiles == struct_unit.heavy_smiles
    
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
    
    exp_output = [(-1.8887, 1.7568, -2.7688), 
                  (-2.3678, 1.2888, -1.7419), (-1.765, 0.3766, -0.9175), 
                  (-0.4846, 0.4733, -0.6949), (0.4424, 1.5847, -1.2336), 
                  (1.6024, 1.7609, -0.3356), (2.3737, 0.5097, -0.1851), 
                  (1.4665, -0.6127, 0.3644), (0.1328, -0.4971, 0.0561), 
                  (1.8655, -1.6442, 1.0535), (3.2048, -1.8132, 1.2835), 
                  (4.1063, -1.5803, 0.486), (-3.3891, 1.5689, -1.4308), 
                  (0.7669, 1.3216, -2.246), (-0.0945, 2.5366, -1.2628), 
                  (2.2151, 2.4716, -0.7355), (2.7962, 0.192, -1.1443), 
                  (3.1861, 0.7105, 0.5183), (-0.4645, -1.2363, 0.3984), 
                  (3.4403, -2.2326, 2.277)]

    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)
    struct_unit = StructUnit(mol_file)
    output = list(struct_unit.get_heavy_coords())
    
    assert exp_output == output
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    