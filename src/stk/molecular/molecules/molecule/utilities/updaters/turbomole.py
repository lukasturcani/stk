import numpy as np
from stk.utilities import periodic_table


def _with_structure_from_turbomole(self, path):
    """
    Return a clone, with its structure taken from a Turbomole file.

    Note that coordinates in ``.coord`` files can be given in Bohr or
    Angstrom, which is handled. Fractional coordinates are not
    currently handled.

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

    :class:`RuntimeError`
        If the the Turbomole file has coordinates defined as fractional
        based on a unit cell.

    """

    bohr_to_ang = 0.5291772105638411

    num_atoms = len(self._atoms)
    with open(path, 'r') as f:
        content = f.readlines()

    for line_number, line in enumerate(content):
        if '$coord' in line:
            if 'angs' in line:
                coord_units = 'angstrom'
            elif 'frac' in line:
                coord_units = 'fractional'
                raise RuntimeError(
                    'Fractional coordinates are not handled currently.'
                )
            elif 'bohr' in line:
                coord_units = 'bohr'
            else:
                coord_units = 'bohr'
            coord_section = []
            coord_section = (
                content[line_number+1:line_number+1+num_atoms]
            )

    # Save all the coords in the file.
    new_coords = []
    for i, line in enumerate(coord_section):
        *coords, element = line.split()
        if element.isnumeric():
            element = periodic_table[int(element)]

        if element != self._atoms[i].__class__.__name__:
            raise RuntimeError(
                f'Atom {i} element does not match file.'
            )

        if coord_units == 'bohr':
            new_coords.append([float(i)*bohr_to_ang for i in coords])
        elif coord_units == 'fractional':
            raise RuntimeError(
                'Fractional coordinates are not handled currently.'
            )
        else:
            new_coords.append([float(i) for i in coords])

    # Check that the correct number of atom
    # lines was present in the file.
    if i+1 != num_atoms:
        raise RuntimeError(
            'The number of atoms lines in the coord file, '
            f'{i+1}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Update the structure.
    new_coords = np.array(new_coords)
    return self._with_position_matrix(new_coords)
