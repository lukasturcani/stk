from os.path import join
import numpy as np
import itertools as it
from scipy.spatial.distance import euclidean
import rdkit.Chem as chem

from ..molecular import StructUnit
from ..convenience_tools import normalize_vector

data_dir = join('data', 'struct_unit', 'amine.mol')
mol = StructUnit(data_dir)
conf = mol.mol.GetConformer()

def test_init():
    f = join('data', 'struct_unit', 'amine2.mol2')
    mol2 = StructUnit(f)
    mol3 = StructUnit(f, 'aldehyde')
    print(mol2.func_grp.name)
    assert len(mol2.functional_group_atoms()) == 3
    assert len(mol3.functional_group_atoms()) == 2

def test_all_bonder_distances():
    conf_ds = []
    for a1, a2 in it.combinations(mol.bonder_ids, 2):
        a1_coord = conf.GetAtomPosition(a1)
        a2_coord = conf.GetAtomPosition(a2)
        conf_ds.append(euclidean(a1_coord, a2_coord))

    assert len(conf_ds) == sum(1 for _ in mol.all_bonder_distances())
    assert np.allclose(list(x for *_, x in mol.all_bonder_distances()),
                       conf_ds, atol=1e-8)

def test_bonder_centroid():
    position = np.array([0.,0.,0.])
    for id_ in mol.bonder_ids:
        position += np.array([*conf.GetAtomPosition(id_)])

    assert np.allclose((position/len(mol.bonder_ids)),
                       mol.bonder_centroid(), atol=1e-8)

def test_bonder_direction_vectors():
    vs = []
    for atom1_id, atom2_id in it.combinations(mol.bonder_ids, 2):
        p1 = np.array([*conf.GetAtomPosition(atom1_id)])
        p2 = np.array([*conf.GetAtomPosition(atom2_id)])

        vs.append(normalize_vector(p1-p2))

    assert len(vs) == sum(1 for _ in mol.bonder_direction_vectors())
    assert np.allclose(
        list(x for *_, x in mol.bonder_direction_vectors()), vs,
            atol=1e-8)

def test_bonder_position_matrix():
        pos_array = np.array([])

        for atom_id in mol.bonder_ids:
            pos_vect = np.array([*conf.GetAtomPosition(atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        m = np.matrix(pos_array.reshape(-1,3).T)

        assert np.allclose(m, mol.bonder_position_matrix(), atol=1e-8)

def test_centroid_centroid_dir_vector():
    c1 = mol.bonder_centroid()
    c2 = mol.centroid()
    assert np.allclose(normalize_vector(c2-c1),
                       mol.centroid_centroid_dir_vector(),
                       atol=1e-8)

def test_functional_group_atoms():
        func_grp_mol = chem.MolFromSmarts(mol.func_grp.fg_smarts)
        assert (mol.mol.GetSubstructMatches(func_grp_mol)  ==
                mol.functional_group_atoms())

def test_json_init():
    bb1 = Molecule.load(('{\n    "mol_block": "\\n     RDKit      '
    '    3D\\n\\n  0  0  0  0  0  0  0  0  0  0999 V3000\\nM  V30'
    ' BEGIN CTAB\\nM  V30 COUNTS 37 39 0 0 0\\nM  V30 BEGIN ATOM\\nM'
    '  V30 1 C 3.1747 2.1211 -2.1979 0\\nM  V30 2 C 2.9418 2.6800 -0'
    '.7878 0\\nM  V30 3 O 2.0164 1.8702 -0.0715 0\\nM  V30 4 C 0.6792'
    ' 2.0573 -0.1762 0\\nM  V30 5 O 0.1723 2.9421 -0.8721 0\\nM  V30'
    ' 6 C -0.1130 1.1117 0.6436 0\\nM  V30 7 C -1.3766 1.4228 1.0308'
    ' 0\\nM  V30 8 N -2.1101 2.5298 0.7135 0\\nM  V30 9 N -2.1025 0.'
    '6126 1.8592 0\\nM  V30 10 C -1.5958 -0.5329 2.4033 0\\nM  V30 1'
    '1 N -2.5051 -1.1161 3.2347 0\\nM  V30 12 C -0.3664 -0.9908 2.08'
    '12 0\\nM  V30 13 C 0.5796 -0.1946 1.1427 0\\nM  V30 14 C 1.1651'
    ' -1.2009 0.1328 0\\nM  V30 15 C 1.4006 -0.9712 -1.2326 0\\nM  V'
    '30 16 C 1.9778 -1.9505 -2.0426 0\\nM  V30 17 C 2.3235 -3.1888 -'
    '1.5054 0\\nM  V30 18 C 2.0736 -3.4524 -0.1607 0\\nM  V30 19 C 1'
    '.4891 -2.4717 0.6397 0\\nM  V30 20 O 1.2100 -2.8331 1.9436 0\\n'
    'M  V30 21 C 0.1697 -2.2572 2.6262 0\\nM  V30 22 O -0.2146 -2.86'
    '02 3.6283 0\\nM  V30 23 H 3.5339 1.0926 -2.1558 0\\nM  V30 24 H'
    ' 2.2540 2.1340 -2.7818 0\\nM  V30 25 H 3.9159 2.7140 -2.7327 '
    '0\\nM  V30 26 H 2.6085 3.7191 -0.8320 0\\nM  V30 27 H 3.8837 2.'
    '6804 -0.2383 0\\nM  V30 28 H -2.7265 2.9469 1.3912 0\\nM  V30 2'
    '9 H -1.6892 3.1510 0.0311 0\\nM  V30 30 H -3.1016 0.6633 1.7683'
    ' 0\\nM  V30 31 H -3.1095 -0.5390 3.7954 0\\nM  V30 32 H -2.2061'
    ' -1.9892 3.6540 0\\nM  V30 33 H 1.3893 0.1581 1.7877 0\\nM  V30'
    ' 34 H 1.1369 -0.0322 -1.6909 0\\nM  V30 35 H 2.1540 -1.7468 -3.'
    '0890 0\\nM  V30 36 H 2.7694 -3.9479 -2.1315 0\\nM  V30 37 H 2.3'
    '141 -4.4175 0.2606 0\\nM  V30 END ATOM\\nM  V30 BEGIN BOND\\nM '
    ' V30 0 1 1 2\\nM  V30 1 1 1 23\\nM  V30 2 1 1 24\\nM  V30 3 1 1'
    ' 25\\nM  V30 4 1 2 3\\nM  V30 5 1 2 26\\nM  V30 6 1 2 27\\nM  V'
    '30 7 1 3 4\\nM  V30 8 2 4 5\\nM  V30 9 1 4 6\\nM  V30 10 2 6 '
    '7\\nM  V30 11 1 6 13\\nM  V30 12 1 7 8\\nM  V30 13 1 7 9\\nM  V'
    '30 14 1 8 28\\nM  V30 15 1 8 29\\nM  V30 16 1 9 10\\nM  V30 17 '
    '1 9 30\\nM  V30 18 1 10 11\\nM  V30 19 2 10 12\\nM  V30 20 1 11'
    ' 31\\nM  V30 21 1 11 32\\nM  V30 22 1 12 13\\nM  V30 23 1 12 2'
    '1\\nM  V30 24 1 13 14\\nM  V30 25 1 13 33\\nM  V30 26 2 14 15\\n'
    'M  V30 27 1 14 19\\nM  V30 28 1 15 16\\nM  V30 29 1 15 34\\nM  '
    'V30 30 2 16 17\\nM  V30 31 1 16 35\\nM  V30 32 1 17 18\\nM  V30'
    ' 33 1 17 36\\nM  V30 34 2 18 19\\nM  V30 35 1 18 37\\nM  V30 36'
    ' 1 19 20\\nM  V30 37 1 20 21\\nM  V30 38 2 21 22\\nM  V30 END'
    ' BOND\\nM  V30 END CTAB\\nM  END\\n\\n$$$$\\n",\n    "key": '
    '[\n        "amine",\n        "InChI=1S/C15H15N3O4/c1-2-21-14(19)'
    '10-9-7-5-3-4-6-8(7)22-15(20)11(9)13(17)18-12(10)16/h3-6,9,18H,'
    '2,16-17H2,1H3/t9-/m0/s1"\n    ],\n    "class": "StructUnit2"\n}')
    )
    

def test_set_bonder_centroid():
    og = mol.bonder_centroid()
    mol.set_bonder_centroid([1,2,3])
    assert np.allclose(mol.bonder_centroid(), [1,2,3], atol=1e-8)
    mol.set_bonder_centroid(og)

def test_untag_atoms():
    f = join('data', 'struct_unit', 'amine2.mol2')
    mol2 = StructUnit(f)
    mol2.untag_atoms()
    for atom in mol2.mol.GetAtoms():
        assert not atom.HasProp('fg')
