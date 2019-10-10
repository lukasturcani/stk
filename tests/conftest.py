import pytest
import stk
import os
from os.path import join
import logging
import sys
from collections import Counter, defaultdict

from .fixtures import *


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
        mol.building_block_vertices = (
            topology_graph.assign_building_blocks_to_vertices(
                building_blocks=building_blocks
            )
        )
        vertex_clones = tuple(
            mol.topology_graph._get_vertex_clones(mol, 1)
        )
        edge_clones = tuple(
            mol.topology_graph._get_edge_clones(1)
        )
        mol._edge_clones = edge_clones

        mol.topology_graph._prepare(mol)
        mol.topology_graph._place_building_blocks(
            mol=mol,
            vertices=vertex_clones,
            edges=edge_clones
        )
        return stk.molecular.reactor.Reactor(mol)

    return inner


@pytest.fixture(scope='session')
def mae_path():
    return join('..', 'data', 'molecule.mae')


@pytest.fixture(scope='session')
def bb_dir():
    return join('..', 'data', 'building_block_init')


@pytest.fixture(scope='session')
def valid_topologies_dir():
    return join('..', 'data', 'valid_topologies')


@pytest.fixture(scope='session')
def valid_cage_dir(valid_topologies_dir):
    return join(valid_topologies_dir, 'cage_topology_tests_output')


@pytest.fixture(scope='session')
def valid_cof_dir(valid_topologies_dir):
    return join(valid_topologies_dir, 'cof_topology_tests_output')


@pytest.fixture(scope='session')
def valid_cyclic_dir(valid_topologies_dir):
    return join(valid_topologies_dir, 'cyclic_topology_tests_output')


@pytest.fixture(scope='session')
def valid_host_guest_dir(valid_topologies_dir):
    return join(
        valid_topologies_dir, 'host_guest_topology_tests_output'
    )


@pytest.fixture(scope='session')
def valid_linear_dir(valid_topologies_dir):
    return join(valid_topologies_dir, 'linear_topology_tests_output')


@pytest.fixture(scope='session')
def valid_rotaxane_dir(valid_topologies_dir):
    return join(valid_topologies_dir, 'rotaxane_topology_tests_output')
