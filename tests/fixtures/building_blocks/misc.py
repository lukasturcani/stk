import pytest
import stk
from os.path import join


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.BuildingBlock('OB(O)c1ccc(B(O)O)nc1', ['boronic_acid'])


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'boronic_acid.sdf')
    )


@pytest.fixture(scope='session')
def diol2():
    return stk.BuildingBlock('Oc1cc2cc(O)c(O)nc2cc1O', ['diol'])


@pytest.fixture(scope='session')
def ring_amine():
    return stk.BuildingBlock(
        smiles='Nc1ccc2cc3cc(N)ccc3cc2c1',
        functional_groups=['ring_amine']
    )


@pytest.fixture(scope='session')
def water():
    return stk.BuildingBlock('[H]O[H]')


@pytest.fixture(scope='session')
def cycle():
    return stk.BuildingBlock('CCCC1CCCCCCCCC1')


@pytest.fixture(scope='session')
def c60():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'c60.pdb')
    )


@pytest.fixture(scope='session')
def chained_c60():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'chained_c60.mol')
    )


@pytest.fixture(scope='session')
def tmp_monodent():
    ligand = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand.func_groups = tuple(i for i in [ligand.func_groups[0]])
    return ligand


@pytest.fixture(scope='session')
def tmp_bident():
    return stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )


@pytest.fixture(scope='session')
def tmp_metal():
    from rdkit.Chem import AllChem as rdkit
    m = rdkit.MolFromSmiles('[Pd+2]')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    metal_coord_info = {
        0: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        1: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        2: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        3: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
    }
    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    return metal
