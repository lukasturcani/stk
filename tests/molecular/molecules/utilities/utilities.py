import numpy as np
import stk

from .building_block import is_clone_building_block
from .constructed_molecule import is_clone_constructed_molecule
from .molecule import is_clone_molecule


def normalize_ids(molecule, ids):
    if ids is None:
        return range(molecule.get_num_atoms())
    elif isinstance(ids, int):
        return (ids,)
    else:
        return ids


def has_same_structure(molecule1, molecule2):
    assert np.all(
        np.equal(
            molecule1.get_position_matrix(),
            molecule2.get_position_matrix(),
        )
    )


def get_displacement_vector(molecule, start_atom, end_atom):
    """
    Get the displacement vector between `start_atom` and `end_atom`.

    """

    position1, position2 = molecule.get_atomic_positions(
        (start_atom, end_atom)
    )
    return position2 - position1


def get_num_atom_ids(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    if atom_ids is None:
        return molecule.get_num_atoms()
    else:
        return len(tuple(atom_ids))


def is_clone(molecule1, molecule2):
    is_building_block1 = isinstance(molecule1, stk.BuildingBlock)
    is_building_block2 = isinstance(molecule2, stk.BuildingBlock)
    if is_building_block1 and is_building_block2:
        is_clone_building_block(molecule1, molecule2)
        return

    is_constructed_molecule1 = isinstance(
        molecule1,
        stk.ConstructedMolecule,
    )
    is_constructed_molecule2 = isinstance(
        molecule2,
        stk.ConstructedMolecule,
    )
    if is_constructed_molecule1 and is_constructed_molecule2:
        is_clone_constructed_molecule(molecule1, molecule2)
        return

    is_clone_molecule(molecule1, molecule2)
