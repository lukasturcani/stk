import pytest
import stk
import os
from os.path import join
import logging
import sys
from collections import Counter, defaultdict

from .fixtures import *


logging.basicConfig(
    format='\n\n%(levelname)s:%(module)s:%(message)s',
    stream=sys.stdout
)
logging.getLogger('stk').setLevel(logging.DEBUG)


# Run tests in a directory so that that generated files are easy to
# delete.
output_dir = 'tests_output'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
os.chdir(output_dir)


def pytest_addoption(parser):
    parser.addoption('--macromodel_path', default='')
    parser.addoption('--mopac_path', default='')
    parser.addoption('--xtb_path', default='')


def pytest_generate_tests(metafunc):
    if 'macromodel_path' in metafunc.fixturenames:
        mm_path = metafunc.config.getoption('macromodel_path')
        metafunc.parametrize('macromodel_path', [mm_path])

    if 'mopac_path' in metafunc.fixturenames:
        mopac_path = metafunc.config.getoption('mopac_path')
        metafunc.parametrize('mopac_path', [mopac_path])

    if 'xtb_path' in metafunc.fixturenames:
        xtb_path = metafunc.config.getoption('xtb_path')
        metafunc.parametrize('xtb_path', [xtb_path])


@pytest.fixture(scope='session')
def make_reactor():

    def inner(building_blocks, topology_graph):
        mol = stk.ConstructedMolecule.__new__(stk.ConstructedMolecule)
        mol.topology_graph = topology_graph
        mol.atoms = []
        mol.bonds = []
        mol.construction_bonds = []
        mol.func_groups = []
        mol.building_block_counter = Counter()
        mol._position_matrix = []
        mol.building_block_vertices = defaultdict(list)
        mol.topology_graph._assign_building_blocks_to_vertices(
            mol=mol,
            building_blocks=building_blocks
        )
        vertex_clones = mol.topology_graph._clone_vertices(mol, 1)
        edge_clones = mol.topology_graph._clone_edges(vertex_clones, 1)
        mol._edge_clones = edge_clones

        mol.topology_graph._prepare(mol)
        mol.topology_graph._place_building_blocks(mol, vertex_clones)
        return stk.molecular.reactor.Reactor(mol)

    return inner


@pytest.fixture(scope='session')
def mae_path():
    return join('..', 'data', 'molecule.mae')


@pytest.fixture(scope='session')
def bb_dir():
    return join('..', 'data', 'building_block_init')
