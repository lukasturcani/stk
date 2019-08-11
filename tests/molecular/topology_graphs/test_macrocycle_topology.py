import os
from os.path import join
import stk


from ..._test_utilities import _test_dump_and_load


test_dir = 'cyclic_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_construction(tmp_macrocycle):
    repeat_units = 3
    monomer_joins = 2*repeat_units

    topology_name = tmp_macrocycle.topology_graph.__class__.__name__
    tmp_macrocycle.write(join(test_dir, f'{topology_name}.mol'))

    assert len(tmp_macrocycle.construction_bonds) == monomer_joins
    assert False


def test_topology(tmp_macrocycle):
    _test_construction(tmp_macrocycle)
    _test_dump_and_load(test_dir, tmp_macrocycle)
