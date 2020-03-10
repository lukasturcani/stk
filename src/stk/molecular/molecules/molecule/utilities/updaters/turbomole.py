from stk.utilities import periodic_table


def _with_structure_from_turbomole(self, path):
    """
    Return a clone, with its structure taken from a Turbomole file.

    Note that coordinates in ``.coord`` files are given in Bohr.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.coord`` file from which the
        structure should be updated.

    Returns
    -------
    :class:`.Molecule`
        A clone with atomic positions found in `path`.

    Raises
    ------
    :class:`RuntimeError`
        If the number of atoms in the file does not match the
        number of atoms in the molecule or if atom elements in the
        file do not agree with the atom elements in the molecule.

    """

    bohr_to_ang = 0.5291772105638411

    with open(path, 'r') as f:
        _, *content, __ = f.readlines()

    # Check the atom count is correct.
    num_atoms = len(self._atoms)
    if len(content) != num_atoms:
        raise RuntimeError(
            'The number of atoms in the coord file, '
            f'{len(content)}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Save all the coords in the file.
    new_coords = []
    for i, line in enumerate(content):
        *coords, element = line.split()
        if element.isnumeric():
            element = periodic_table[int(element)]

        if element != self._atoms[i].__class__.__name__:
            raise RuntimeError(
                f'Atom {i} element does not match file.'
            )

        new_coords.append([float(i)*bohr_to_ang for i in coords])

    # Check that the correct number of atom
    # lines was present in the file.
    if i+1 != num_atoms:
        raise RuntimeError(
            'The number of atoms lines in the coord file, '
            f'{i+1}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Update the structure.
    return self._with_position_matrix(new_coords)
