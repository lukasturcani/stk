from os.path import join
import os
import numpy as np


test_dir = 'macorcycle_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_cycle_atoms(cycle_su):
    assert sorted(cycle_su.cycle_atoms()) == list(range(3, 13))


def test_write_cycle_coords(cycle_su):
    path = join(f'{test_dir}', 'cycle_coords.xyz')
    cycle_su.write_cycle_coords(path)

    with open(path, 'r') as f:
        atom_count, _, *content = f.read().split('\n')

    cycle_atoms = cycle_su.cycle_atoms()

    # Make sure the atom count was written correctly.
    assert int(atom_count) == len(cycle_atoms)

    # Make sure the correct number of atoms was written to the file.
    assert len(content) == len(cycle_atoms)

    # Make sure the coordinates were written correctly.
    cycle_atoms = set(cycle_atoms)
    for line in content:
        atom_id, *coord = line.split()
        atom_id = int(atom_id)
        coord = [float(i) for i in coord]
        assert atom_id in cycle_atoms
        # Make sure there are no duplicates.
        cycle_atoms.remove(atom_id)

        assert np.allclose(
            a=cycle_su.atom_coords(atom_id),
            b=coord,
            atol=1e-8
        )
