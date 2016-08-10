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

    """
    mol_file = next(x for x in get_mol_file() 
                                        if 'amine3f_14.mol' in x)   
    mol = StructUnit(mol_file)
    for x,y,z in mol.get_heavy_coords():
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert isinstance(z, float)
    
    assert len(list(mol.get_heavy_coords())) == 32
    
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
    exp_smiles = ("[H][N]1/[C](=[N]/[Y])[C]([H])([H])[N]([H])[C]([H])"
                    "([H])/[C]1=[N]\[Y]")
    mol_file = next(x for x in get_mol_file() 
                                        if 'aldehyde2f_3.mol' in x)   
    mol = StructUnit(mol_file)
    print(mol.heavy_smiles)
    assert mol.heavy_smiles == exp_smiles
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    