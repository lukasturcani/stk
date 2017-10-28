#!/usr/bin/env python3
"""
Module containing all general purpose functions shared by other modules.

LOG
---
27/07/17
    Fixed the cartesian coordinates -> fractional coordinates -> cartesian
    coordinates conversion related functions, creation of lattice array
    from unit cell parameters (triclinic system: so applicable to any)
    and conversion back to unit cell parameters. WORKS! inspiration from:
    http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm

26/07/17
    Changed the way bonds are determined. Now, rather then fixed value
    a formula and covalent radii are used as explained in the Elemental_Radii
    spreadsheet (see tables module).

TO DO LIST
----------

- In the find_windows() function, maybe change the way the EPS value for
  the DBSCAN() is estimates. Need to look how the distances change with the
  increase in size of the sampling sphere.
"""

import numpy as np
from copy import deepcopy
from multiprocessing import Pool
from scipy.optimize import brute, fmin, minimize
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import KDTree

from .tables import (
    atomic_mass, atomic_vdw_radius, opls_atom_keys, atomic_covalent_radius
    )


class _AtomKeyError(Exception):
    def __init__(self, message):
        self.message = message


class _AtomKeyConflict(Exception):
    def __init__(self, message):
        self.message = message


class _ForceFieldError(Exception):
    def __init__(self, message):
        self.message = message


class _FunctionError(Exception):
    def __init__(self, message):
        self.message = message


def is_number(number):
    """
    Return True if an object is a number - can be converted into a float.

    Parameters
    ----------
    number : any

    Returns
    -------
    bool
        True if input is a float convertable (a number), False otherwise.

    """
    try:
        float(number)
        return True
    except ValueError:
        return False


def unique(input_list):
    """
    Return a list of unique items (similar to set functionality).

    Parameters
    ----------
    input_list : list
        A list containg some items that can occur more than once.

    Returns
    -------
    list
        A list with only unique occurances of an item.

    """
    output = []
    for item in input_list:
        if item not in output:
            output.append(item)
    return output


def make_JSON_serializable(obj):
    """
    Return a dictionary with all arrays (three lvls down) changed to lists.

    It takes a dictionary (and subdictionaries up to third level down) and
    searches for numpy.ndarrays that are not json serializable. It exchanges
    these arrays into lists.

    Parameters
    ----------
    obj : dictionary
        A dictionary that is going to be dumped as json.

    Modifies
    --------
    obj : dictionary
        The input dictionary.

    Returns
    -------
    dictionary
        A processed dictionary without numpy arrays.

    """
    for key in obj.keys():
        if isinstance(obj[key], np.ndarray):
            obj[key] = obj[key].tolist()
        try:
            if obj[key].keys() is not None:
                for key2 in obj[key].keys():
                    if isinstance(obj[key][key2], np.ndarray):
                        obj[key][key2] = obj[key][key2].tolist()
                    try:
                        if obj[key][key2].keys() is not None:
                            for key3 in obj[key][key2].keys():
                                if isinstance(obj[key][key2][key3],
                                              np.ndarray):
                                    obj[key][key2][key3] = obj[key][key2][
                                        key3].tolist()
                    except AttributeError:
                        pass
        except AttributeError:
            pass
    return obj


def distance(a, b):
    """
    Return the distance between two vectors (points) a and b.

    Parameters
    ----------
    a : numpy.ndarray
        First vector.
    b : numpy.ndarray
        Second vector.

    Returns
    -------
    numpy.float64
        A distance between two vectors (points).

    """
    return (np.sum((a - b)**2))**0.5


def molecular_weight(elements):
    """
    Return molecular weight of a molecule.

    Parameters
    ----------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    Returns
    -------
    numpy.float64
        A molecular weight of a molecule.

    """
    return (np.array([atomic_mass[i.upper()] for i in elements]).sum())


def center_of_coor(coordinates):
    """
    Return the centre of coordinates.

    Parameters
    ----------
    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    Returns
    -------
    numpy.ndarray
        An 1d array with coordinates of the centre of coordinates excluding
        elements' masses.

    """
    return (np.sum(coordinates, axis=0) / coordinates.shape[0])


def center_of_mass(elements, coordinates):
    """
    Return the centre of mass (COM).

    Parameters
    ----------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    Returns
    -------
    numpy.ndarray
        An 1d array with coordinates of the centre of mass including elements'
        masses.

    """
    mass = molecular_weight(elements)
    mass_array = np.array([[atomic_mass[i.upper()]] * 3 for i in elements])
    mass_coordinates = coordinates * mass_array
    return (np.sum(mass_coordinates, axis=0) / np.array([mass, mass, mass]))


def compose_atom_list(*args):
    """
    Return an `atom list` from elements and/or atom ids and coordinates.

    An `atom list` is a special object that some pyWINDOW functions uses.
    It is a nested list of lists with each individual list containing:
        version 1 : [[element, coordinates (x, y, z)], ...]
        version 2 : [[element, atom key, coordinates (x, y, z)], ...]
    They work better for molecular re-building than two separate arrays for
    elements and coordinates do.

    Parameters (version 1)
    -----------------------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    Parameters (version 2)
    -----------------------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    atom_ids : numpy.ndarray
        An array of all forcfield dependent atom keys (type:str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    Returns
    -------
    list
        Version 1 or version 2 atom list depending on input parameters.

    Raises
    ------
    _FunctionError : Exception
        Raised when wrong number of parameters is passed to the function.

    """
    if len(args) == 2:
        atom_list = [[
            i[0],
            round(float(i[1]), 8),
            round(float(i[2]), 8),
            round(float(i[3]), 8),
        ] for i in np.concatenate(
            (args[0].reshape(-1, 1), args[1]), axis=1)]
    elif len(args) == 3:
        atom_list = [
            [
                i[0],
                i[1],
                round(float(i[2]), 8),
                round(float(i[3]), 8),
                round(float(i[4]), 8),
            ]
            for i in np.concatenate(
                (np.concatenate(
                    (args[0].reshape(-1, 1), args[1].reshape(-1, 1)
                     ), axis=1), args[2]),
                axis=1)
        ]
    else:
        raise _FunctionError(
            "The compose_atom_list() function accepts only 2 or 3 arguments.")
    return atom_list


def decompose_atom_list(atom_list):
    """
    Return elements and/or atom ids and coordinates from an `atom list`.

    Depending on input type of an atom list (version 1 or 2)
        version 1 : [[element, coordinates (x, y, z)], ...]
        version 2 : [[element, atom key, coordinates (x, y, z)], ...]
    the function reverses what pywindow.utilities.compose_atom_list() do.

    Parameters
    ----------
    atom_list : list
        A nested list of lists (version 1 or 2)

    Returns (version 1)
    -------------------
    (numpy.ndarray, numpy.ndarray)
        A touple of elements and coordinates arrays.

    Returns (version 2)
    -------------------
    (numpy.ndarray, numpy.ndarray, numpy.ndarray)
        A touple of elements, atom keys and coordinates arrays.

    """
    transpose = list(zip(*atom_list))
    if len(transpose) == 4:
        elements = np.array(transpose[0])
        array_a = np.array(transpose[1]).reshape(-1, 1)
        array_b = np.array(transpose[2]).reshape(-1, 1)
        array_c = np.array(transpose[3]).reshape(-1, 1)
        array_ab = np.concatenate((array_a, array_b), axis=1)
        coordinates = np.concatenate((array_ab, array_c), axis=1)
        return elements, coordinates
    elif len(transpose) == 5:
        elements = np.array(transpose[0])
        atom_ids = np.array(transpose[1])
        array_a = np.array(transpose[2]).reshape(-1, 1)
        array_b = np.array(transpose[3]).reshape(-1, 1)
        array_c = np.array(transpose[4]).reshape(-1, 1)
        array_ab = np.concatenate((array_a, array_b), axis=1)
        coordinates = np.concatenate((array_ab, array_c), axis=1)
        return elements, atom_ids, coordinates
    else:
        raise _FunctionError(
            "The decompose_atom_list() function accepts only list of lists "
            " with only 4 or 5 items per sublist.")


def dlf_notation(atom_key):
    """Return element for atom key using DL_F notation."""
    split = list(atom_key)
    element = ''
    number = False
    count = 0
    while number is False:
        element = "".join((element, split[count]))
        count += 1
        if is_number(split[count]) is True:
            number = True
    # In case of for example Material Studio output, integers can also be
    # in the beginning of the string. As the dlf_notation decipher function
    # is very general in use, we have to make sure these integers are deleted.
    # In standard DL_F notation the string will never start with integer so it
    # will not affect the functionality towards it.
    # EDIT2: also the '?' atoms, you can delete them manually or somewhere else
    element = "".join(i for i in element if not is_number(i))
    element = "".join(i for i in element if i != '?')
    return element


def opls_notation(atom_key):
    """Return element for OPLS forcefield atom key."""
    # warning for Ne, He, Na types overlap
    conflicts = ['ne', 'he', 'na']
    if atom_key in conflicts:
        raise _AtomKeyConflict((
            "One of the OPLS conflicting "
            "atom_keys has occured '{0}'. "
            "For how to solve this issue see the manual or "
            "MolecularSystem._atom_key_swap() doc string.").format(atom_key))
    for element in opls_atom_keys:
        if atom_key in opls_atom_keys[element]:
            return element
    # In case if atom_key was not found in the OPLS keys dictionary
    raise _AtomKeyError((
        "OPLS atom key {0} was not found in OPLS keys dictionary.").format(
            atom_key))


def decipher_atom_key(atom_key, forcefield):
    """
    Return element for deciphered atom key.

    This functions checks if the forcfield specified by user is supported
    and passes the atom key to the appropriate function for deciphering.

    Parameters
    ----------
    atom_key : str
        The atom key which is to be deciphered.

    forcefield : str
        The forcefield to which the atom key belongs to.

    Returns
    -------
    str
        A string that is the periodic table element equvalent of forcefield
        atom key.

    """
    load_funcs = {
        'DLF': dlf_notation,
        'DL_F': dlf_notation,
        'OPLS': opls_notation,
        'OPLSAA': opls_notation,
        'OPLS2005': opls_notation,
        'OPLS3': opls_notation,
    }
    if forcefield.upper() in load_funcs.keys():
        return load_funcs[forcefield.upper()](atom_key)
    else:
        raise _ForceFieldError(
            ("Unfortunetely, '{0}' forcefield is not supported by pyWINDOW."
             " For list of supported forcefields see User's Manual or "
             "MolecularSystem._decipher_atom_keys() function doc string."
             ).format(forcefield))


def shift_com(elements, coordinates, com_adjust=np.zeros(3)):
    """
    Return coordinates translated by some vector.

    Parameters
    ----------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    com_adjust : numpy.ndarray (default = [0, 0, 0])

    Returns
    -------
    numpy.ndarray
        Translated array of molecule's coordinates.

    """
    com = center_of_mass(elements, coordinates)
    com = np.array([com - com_adjust] * coordinates.shape[0])
    return coordinates - com


def max_dim(elements, coordinates):
    """
    Return the maximum diameter of a molecule.

    Parameters
    ----------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    Returns
    -------

    """
    atom_vdw_vertical = np.matrix(
        [[atomic_vdw_radius[i.upper()]] for i in elements])
    atom_vdw_horizontal = np.matrix(
        [atomic_vdw_radius[i.upper()] for i in elements])
    dist_matrix = euclidean_distances(coordinates, coordinates)
    vdw_matrix = atom_vdw_vertical + atom_vdw_horizontal
    re_dist_matrix = dist_matrix + vdw_matrix
    final_matrix = np.triu(re_dist_matrix)
    i1, i2 = np.unravel_index(final_matrix.argmax(), final_matrix.shape)
    maxdim = final_matrix[i1, i2]
    return i1, i2, maxdim


def void_diameter(elements, coordinates, com=None):
    """Return void diameter of a molecule."""
    if com is None:
        com = center_of_mass(elements, coordinates)
    atom_vdw = np.array([[atomic_vdw_radius[x.upper()]] for x in elements])
    dist_matrix = euclidean_distances(coordinates, com.reshape(1, -1))
    re_dist_matrix = dist_matrix - atom_vdw
    index = np.argmin(re_dist_matrix)
    voidd = re_dist_matrix[index][0] * 2
    return (voidd, index)


def correct_void_diameter(com, *params):
    """Return negative of a void diameter. (optimisation function)."""
    elements, coordinates = params
    return (-void_diameter(elements, coordinates, com)[0])


def opt_void_diameter(elements, coordinates, bounds=None, **kwargs):
    """Return optimised void diameter and it's COM."""
    args = elements, coordinates
    com = center_of_mass(elements, coordinates)
    if bounds is None:
        void_r = void_diameter(elements, coordinates, com=com)[0] / 2
        bounds = (
            (com[0]-void_r, com[0]+void_r),
            (com[1]-void_r, com[1]+void_r),
            (com[2]-void_r, com[2]+void_r)
        )
    minimisation = minimize(
        correct_void_diameter, x0=com, args=args, bounds=bounds)
    voidd = void_diameter(elements, coordinates, com=minimisation.x)
    return (voidd[0], voidd[1], minimisation.x)


def void_volume(void_radius):
    """Return the volume of a spherical void."""
    return (4 / 3 * np.pi * void_radius**3)


def unit_cell_to_lattice_array(cryst):
    """Return parallelpiped unit cell lattice matrix."""
    a_, b_, c_, alpha, beta, gamma = cryst
    # Convert angles from degrees to radians.
    r_alpha = np.deg2rad(alpha)
    r_beta = np.deg2rad(beta)
    r_gamma = np.deg2rad(gamma)
    # Calculate unit cell volume that is neccessary.
    volume = a_ * b_ * c_ * (
        1 - np.cos(r_alpha)**2 - np.cos(r_beta)**2 - np.cos(r_gamma)**2 + 2 *
        np.cos(r_alpha) * np.cos(r_beta) * np.cos(r_gamma))**0.5
    # Create the orthogonalisation Matrix (M^-1) - lattice matrix
    a_x = a_
    a_y = b_ * np.cos(r_gamma)
    a_z = c_ * np.cos(r_beta)
    b_x = 0
    b_y = b_ * np.sin(r_gamma)
    b_z = c_ * (
        np.cos(r_alpha) - np.cos(r_beta) * np.cos(r_gamma)) / np.sin(r_gamma)
    c_x = 0
    c_y = 0
    c_z = volume / (a_ * b_ * np.sin(r_gamma))
    lattice_array = np.array(
        [[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]])
    return lattice_array


def lattice_array_to_unit_cell(lattice_array):
    """Return crystallographic param. from unit cell lattice matrix."""
    cell_lengths = np.sqrt(np.sum(lattice_array**2, axis=0))
    gamma_r = np.arccos(lattice_array[0][1] / cell_lengths[1])
    beta_r = np.arccos(lattice_array[0][2] / cell_lengths[2])
    alpha_r = np.arccos(
        lattice_array[1][2] * np.sin(gamma_r) / cell_lengths[2]
        + np.cos(beta_r) * np.cos(gamma_r)
        )
    cell_angles = [
        np.rad2deg(alpha_r), np.rad2deg(beta_r), np.rad2deg(gamma_r)
        ]
    return np.append(cell_lengths, cell_angles)


def volume_from_lattice_array(lattice_array):
    """Return unit cell's volume from lattice matrix."""
    return np.linalg.det(lattice_array)


def volume_from_cell_parameters(cryst):
    """Return unit cell's volume from crystallographic parameters."""
    return volume_from_lattice_array(unit_cell_to_lattice_array(cryst))


def fractional_from_cartesian(coordinate, lattice_array):
    """Return a fractional coordinate from a cartesian one."""
    deorthogonalisation_M = np.matrix(np.linalg.inv(lattice_array))
    fractional = deorthogonalisation_M * coordinate.reshape(-1, 1)
    return np.array(fractional.reshape(1, -1))


def cartisian_from_fractional(coordinate, lattice_array):
    """Return cartesian coordinate from a fractional one."""
    orthogonalisation_M = np.matrix(lattice_array)
    orthogonal = orthogonalisation_M * coordinate.reshape(-1, 1)
    return np.array(orthogonal.reshape(1, -1))


def cart2frac_all(coordinates, lattice_array):
    """Convert all cartesian coordinates to fractional."""
    frac_coordinates = deepcopy(coordinates)
    for coord in range(frac_coordinates.shape[0]):
        frac_coordinates[coord] = fractional_from_cartesian(
            frac_coordinates[coord], lattice_array)
    return frac_coordinates


def frac2cart_all(frac_coordinates, lattice_array):
    """Convert all fractional coordinates to cartesian."""
    coordinates = deepcopy(frac_coordinates)
    for coord in range(coordinates.shape[0]):
        coordinates[coord] = cartisian_from_fractional(coordinates[coord],
                                                       lattice_array)
    return coordinates


def create_supercell(system, supercell=[[-1, 1], [-1, 1], [-1, 1]]):
    """Create a supercell."""
    if 'lattice' not in system.keys():
        matrix = unit_cell_to_lattice_array(system['unit_cell'])
    else:
        matrix = system['lattice']
    coordinates = deepcopy(system['coordinates'])
    multiplication_matrices = []
    for a_ in range(supercell[0][0], supercell[0][1] + 1):
        for b_ in range(supercell[1][0], supercell[1][1] + 1):
            for c_ in range(supercell[2][0], supercell[2][1] + 1):
                mult_matrix = np.array([[a_, b_, c_]])
                mult_matrix = np.repeat(
                    mult_matrix, coordinates.shape[0], axis=0)
                multiplication_matrices.append(mult_matrix)
    frac_coordinates = cart2frac_all(coordinates, matrix)
    updated_coordinates = []
    for mat in multiplication_matrices:
        updated_coor = frac_coordinates + mat
        updated_coordinates.append(updated_coor)
    supercell_frac_coordinates = np.concatenate(updated_coordinates, axis=0)
    supercell_coordinates = frac2cart_all(supercell_frac_coordinates, matrix)
    # Now for each new cell in the supercell we need to repeat the
    # elements array so that it maches
    new_elements = deepcopy(system['elements'])
    new_ids = deepcopy(system['atom_ids'])
    for i in range(len(updated_coordinates) - 1):
        new_elements = np.concatenate((new_elements, system['elements']))
        new_ids = np.concatenate((new_ids, system['atom_ids']))
    cryst = lattice_array_to_unit_cell(matrix)
    supercell_system = {
        'elements': new_elements,
        'atom_ids': new_ids,
        'coordinates': supercell_coordinates,
        'unit_cell': cryst,
        'lattice': matrix,
    }
    return supercell_system


def is_inside_polyhedron(point, polyhedron):
    if polyhedron.shape == (1, 6):
        matrix = unit_cell_to_lattice_array(polyhedron)
    if polyhedron.shape == (3, 3):
        matrix = polyhedron
    frac_coordinates = cart2frac_all(point, matrix.T)
    if point[0] <= 1.000 and point[1] <= 1.000 and point[2] <= 1.000:
        return True
    else:
        return False


def normal_vector(origin, vectors):
    """Return normal vector for two vectors with same origin."""
    return np.cross(vectors[0] - origin, vectors[1] - origin)


def discrete_molecules(system, rebuild=None, tol=0.4):
    """
    Decompose molecular system into individual discreet molecules.

    New formula for bonds: (26/07/17)
        The two atoms, x and y, are considered bonded if the distance between
        them, calculated with distance matrix, is within the ranges:
              Rcov(x) + Rcov(y) - t < R(x,y) <  Rcov(x) + Rcov(y) + t
        where Rcov is the covalent radius and the tolarenace (t) is set to
        0.4 Angstrom.
    """
    # First we check which operation mode we use.
    #    1) Non-periodic MolecularSystem.
    #    2) Periodic MolecularSystem without rebuilding.
    #    3) Periodic Molecular system with rebuilding (supercell provided).
    if rebuild is not None:
        mode = 3
    else:
        if 'unit_cell' in system.keys() or 'lattice' in system.keys():
            mode = 2
        else:
            mode = 1
    # We create a list containing all atoms, theirs periodic elements and
    # coordinates. As this process is quite complicated, we need a list
    # which we will gradually be reducing.
    try:
        elements = system['elements']
        coordinates = system['coordinates']
    except KeyError:
        raise _FunctionError(
            "The 'elements' key is missing in the 'system' dictionary "
            "attribute of the MolecularSystem object. Which means, you need to"
            " decipher the forcefield based atom keys first (see manual)."
        )
    coordinates = system['coordinates']
    args = (elements, coordinates)
    adj = 0
    # If there are forcefield 'atom ids' as well we will retain them.
    if 'atom_ids' in system.keys():
        atom_ids = system['atom_ids']
        args = (elements, atom_ids, coordinates)
        adj = 1
    atom_list = compose_atom_list(*args)
    atom_coor = decompose_atom_list(atom_list)[1 + adj]
    # Scenario 1: We load a non-periodic MolecularSystem.
    # We will not have 'unit_cell' nor 'lattice' keywords in the dictionary
    # and also we do not do any re-building.
    # Scenario 2: We load a periodic MolecularSystem. We want to only Extract
    # complete molecules that do not have been affected by the periodic
    # boundary.
    # Scenario 3: We load a periodic Molecular System. We want it to be rebuild
    # therefore, we also provide a supercell.
    # Scenarios 2 and 3 require a lattice and also their origin is at origin.
    # Scenario 1 should have the origin at the center of mass of the system.
    if mode == 2 or mode == 3:
        # Scenarios 2 or 3.
        origin = np.array([0., 0., 0.])
        if 'lattice' not in system.keys():
            matrix = unit_cell_to_lattice_array(system['unit_cell'])
        else:
            matrix = system['lattice']
        pseudo_origin_frac = np.array([0.25, 0.25, 0.25])
        pseudo_origin = cartisian_from_fractional(pseudo_origin_frac, matrix)
        # If a supercell is also provided that encloses the unit cell for the
        # reconstruction of the molecules through the periodic boundary.
        if rebuild is not None:
            selements = rebuild['elements']
            sids = rebuild['atom_ids']
            scoordinates = rebuild['coordinates']
            satom_list = compose_atom_list(selements, sids, scoordinates)
            satom_coor = decompose_atom_list(satom_list)[1 + adj]
        # There is one more step. We need to sort out for all the
        # reconstructed molecules, which are the ones that belong to the
        # unit cell. As we did the reconstruction to every chunk in the unit
        # cell we have now some molecules that belong to neighbouring cells.
        # The screening is simple. If the COM of a molecule translated to
        # fractional coordinates (so that it works for parallelpiped) is
        # within the unit cell boundaries <0, 1> then it's it. There is
        # an exception, for the trajectories, very often the unit cell
        # is centered at origin. Therefore we need to use <-0.5, 0.5>
        # boundary. We will simply decide which is the case by calculating
        # the centre of mass of the whole system.
        system_com = center_of_mass(elements, coordinates)
        if np.allclose(system_com, origin, atol=1e-00):
            boundary = np.array([-0.5, 0.5])
        else:
            boundary = np.array([0., 1.])
    else:
        # Scenario 1.
        pseudo_origin = center_of_mass(elements, coordinates)
    # Here the final discrete molecules will be stored.
    molecules = []
    # Exceptions. Usually end-point atoms that create single bonds or
    # just a separate atoms in the system.
    exceptions = ['H', 'CL', 'BR', 'F', 'HE', 'AR', 'NE', 'KR', 'XE', 'RN']
    # The upper limit for distances analysed for bonds will be assigned for
    # a given system (to save time). We take set('elements') and then find
    # the largest R(cov) in the system and set the max_dist as a double
    # of it plus the 150% tolerance (tol).
    set_of_elements = set(system['elements'])
    max_r_cov = max([atomic_covalent_radius[i.upper()] for i in set_of_elements])
    max_dist = 2 * max_r_cov + tol
    # We continue untill all items in the list have been analysed and popped.
    while atom_list:
        inside_atoms_heavy = [
            i for i in atom_list if i[0].upper() not in exceptions
        ]
        if inside_atoms_heavy:
            # Now we create an array of atom coordinates. It does seem
            # somehow counter-intuitive as this is what we started with
            # and made it into a list. But, in my opinion it's the only
            # way to do it. It's hard to control and delete items in two
            # separate arrays that we started with and we don't want
            # atoms already assigned in our array for distance matrix.
            inside_atoms_coord_heavy = decompose_atom_list(inside_atoms_heavy)[
                1 + adj]
            dist_matrix = euclidean_distances(inside_atoms_coord_heavy,
                                              pseudo_origin.reshape(1, -1))
            atom_index_x, _ = np.unravel_index(dist_matrix.argmin(),
                                               dist_matrix.shape)
            # Added this so that lone atoms (even if heavy) close to the
            # periodic boundary are not analysed, as they surely have matching
            # symmetry equivalence that bind to a bigger atom cluster inside
            # the unit_cell.
            potential_starting_point = inside_atoms_heavy[atom_index_x]
            pot_arr = np.array(potential_starting_point[1 + adj:])
            dist_matrix = euclidean_distances(
                atom_coor, pot_arr.reshape(1, -1)
                )
            idx = (dist_matrix > 0.1) * (dist_matrix < max_dist)
            if len(idx) < 1:
                pass
            else:
                working_list = [potential_starting_point]
        else:
            # Safety check.
            break
        final_molecule = []
        while working_list:
            working_list_temp = []
            try:
                atom_coor = decompose_atom_list(atom_list)[1 + adj]
            except _FunctionError:
                atom_coor = None
            for i in working_list:
                if i[0].upper() not in exceptions:
                    # It's of GREATEST importance that the i_arr variable
                    # is assigned here before entering the atom_coor loop.!
                    # Otherwise it will not be re-asigned when the satom_list
                    # still iterates, but the atom_list is already empty...
                    i_arr = np.array(i[1 + adj:])
                    if atom_coor is not None:
                        dist_matrix = euclidean_distances(
                            atom_coor, i_arr.reshape(1, -1)
                            )
                        idx = (dist_matrix > 0.1) * (dist_matrix < max_dist)
                        neighbours_indexes = np.where(idx)[0]
                        for j in neighbours_indexes:
                            j_arr = np.array(atom_coor[j])
                            r_i_j = distance(i_arr, j_arr)
                            r_cov_i_j = atomic_covalent_radius[i[0].upper()] + atomic_covalent_radius[atom_list[j][0].upper()]
                            if r_cov_i_j - tol < r_i_j < r_cov_i_j + tol:
                                working_list_temp.append(atom_list[j])
                    if rebuild is not None:
                        sdist_matrix = euclidean_distances(
                            satom_coor, i_arr.reshape(1, -1))
                        sidx = (sdist_matrix > 0.1) * (sdist_matrix < max_dist)
                        sneighbours_indexes = np.where(sidx)[0]
                        for j in sneighbours_indexes:
                            if satom_list[j] in atom_list:
                                pass
                            else:
                                j_arr = np.array(satom_coor[j])
                                r_i_j = distance(i_arr, j_arr)
                                r_cov_i_j = atomic_covalent_radius[
                                    i[0].upper()
                                    ] + atomic_covalent_radius[satom_list[j][0].upper()]
                                if r_cov_i_j - tol < r_i_j < r_cov_i_j + tol:
                                    working_list_temp.append(satom_list[j])
                    final_molecule.append(i)
                else:
                    final_molecule.append(i)
            for i in working_list:
                try:
                    atom_list.remove(i)
                except ValueError:
                    pass
            # We empty the working list as all the items were analysed
            # and moved to the final_molecule list.
            working_list = []
            # We make sure there are no duplicates in the working_list_temp.
            working_list_temp = unique(working_list_temp)
            # Now we move the entries from the temporary working list
            # to the working list for looping analysys.
            for i in working_list_temp:
                # We make sure that only new and unassigned atoms are
                # being transfered.
                if i not in final_molecule:
                    working_list.append(i)
        final_molecule_dict = {}
        final_molecule_dict['elements'] = np.array(
            [x[0] for x in final_molecule], dtype='str')
        final_molecule_dict['coordinates'] = np.array(
            [[*xyz[1 + adj:]] for xyz in final_molecule])
        if adj == 1:
            final_molecule_dict['atom_ids'] = np.array(
                [x[1] for x in final_molecule], dtype='str')
        # In general we always want the molecule so the initial bool_ is True.
        bool_ = True
        # But, for periodic only if the molecule is in the initial unit cell.
        if rebuild is not None:
            com = center_of_mass(final_molecule_dict['elements'],
                                 final_molecule_dict['coordinates'])
            com_frac = fractional_from_cartesian(com, matrix)[0]
            # If we don't round the numerical errors will come up.
            com_frac_round = np.around(com_frac, decimals=8)
            bool_ = np.all(np.logical_and(com_frac_round >= boundary[0],
                                          com_frac_round < boundary[1]),
                           axis=0)
        if bool(bool_) is True:
            molecules.append(final_molecule_dict)
    return molecules


def angle_between_vectors(x, y):
    """Calculate the angle between two vectors x and y."""
    first_step = abs(x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) / (
        np.sqrt(x[0]**2 + x[1]**2 + x[2]**2) *
        np.sqrt(y[0]**2 + y[1]**2 + y[2]**2))
    second_step = np.arccos(first_step)
    return (second_step)


def vector_analysis(vector, coordinates, elements_vdw, increment=1.0):
    """Analyse a sampling vector's path for window analysis purpose."""
    # Calculate number of chunks if vector length is divided by increment.
    chunks = int(np.linalg.norm(vector) // increment)
    # Create a single chunk.
    chunk = vector / chunks
    # Calculate set of points on vector's path every increment.
    vector_pathway = np.array([chunk * i for i in range(chunks + 1)])
    analysed_vector = np.array([
        np.amin(
            euclidean_distances(coordinates, i.reshape(1, -1)) - elements_vdw)
        for i in vector_pathway
    ])
    if all(i > 0 for i in analysed_vector):
        pos = np.argmin(analysed_vector)
        # As first argument we need to give the distance from the origin.
        dist = np.linalg.norm(chunk * pos)
        return np.array(
            [dist, analysed_vector[pos] * 2, *chunk * pos, *vector])


def optimise_xy(xy, *args):
    """Return negative void diameter for x and y coordinates optimisation."""
    z, elements, coordinates = args
    window_com = np.array([xy[0], xy[1], z])
    return -void_diameter(elements, coordinates, com=window_com)[0]


def optimise_z(z, *args):
    """Return void diameter for coordinates optimisation in z direction."""
    x, y, elements, coordinates = args
    window_com = np.array([x, y, z])
    return void_diameter(elements, coordinates, com=window_com)[0]


def window_analysis(window,
                    elements,
                    coordinates,
                    elements_vdw,
                    increment2=0.1,
                    z_bounds=[None, None],
                    lb_z=True,
                    **kwargs):
    """
    Return window diameter and window's centre.

    Parameters
    ----------
    widnow: list

    elements: numpy.array

    coordinates: numpy.array

    elements_vdw: numpy.array

    step: float

    """
    # Copy the coordinates as we will manipulate them.
    coordinates = deepcopy(coordinates)
    # Find the vector with the largest window sampling diameter from the pool.
    vector_ = window[window.argmax(axis=0)[1]][5:8]
    vector_analysed = vector_analysis(
        vector_, coordinates, elements_vdw, increment=increment2)
    # A safety check, if the refined analysis give None we end the function.
    if vector_analysed is not None:
        pass
    else:
        return None
    vector = vector_analysed[5:8]
    # Unit vectors.
    vec_a = [1, 0, 0]
    vec_b = [0, 1, 0]
    vec_c = [0, 0, 1]
    # Angles needed for rotation (in radians) to rotate and translate the
    # molecule for the vector to become the Z-axis.
    angle_1 = angle_between_vectors(np.array([vector[0], vector[1], 0]), vec_a)
    angle_2 = angle_between_vectors(vector, vec_c)
    # Depending in which cartesian coordinate system area the vector is
    # We need a rotation into a different direction and by different value.
    if vector[0] >= 0 and vector[1] >= 0 and vector[2] >= 0:
        angle_1 = -angle_1
        angle_2 = -angle_2
    if vector[0] < 0 and vector[1] >= 0 and vector[2] >= 0:
        angle_1 = np.pi * 2 + angle_1
        angle_2 = angle_2
    if vector[0] >= 0 and vector[1] < 0 and vector[2] >= 0:
        angle_1 = angle_1
        angle_2 = -angle_2
    if vector[0] < 0 and vector[1] < 0 and vector[2] >= 0:
        angle_1 = np.pi * 2 - angle_1
    if vector[0] >= 0 and vector[1] >= 0 and vector[2] < 0:
        angle_1 = -angle_1
        angle_2 = np.pi + angle_2
    if vector[0] < 0 and vector[1] >= 0 and vector[2] < 0:
        angle_2 = np.pi - angle_2
    if vector[0] >= 0 and vector[1] < 0 and vector[2] < 0:
        angle_2 = angle_2 + np.pi
    if vector[0] < 0 and vector[1] < 0 and vector[2] < 0:
        angle_1 = -angle_1
        angle_2 = np.pi - angle_2
    # Rotation matrix for rotation around Z-axis with angle_1.
    rotation_around_z = np.array([[np.cos(angle_1), -np.sin(angle_1), 0],
                                  [np.sin(angle_1), np.cos(angle_1), 0],
                                  [0, 0, 1]])
    # Rotate the whole molecule around with rotation_around_z.
    coordinates = np.array([np.dot(rotation_around_z, i) for i in coordinates])
    # Rotation matrix for rotation around Y-axis with angle_2
    rotation_around_y = np.array([[np.cos(angle_2), 0, np.sin(angle_2)],
                                  [0, 1, 0],
                                  [-np.sin(angle_2), 0, np.cos(angle_2)]])
    # Rotate the whole molecule around with rotation_around_y.
    coordinates = np.array([np.dot(rotation_around_y, i) for i in coordinates])
    # Third step is translation. We are now at [0, 0, -z].
    # We shift the molecule so that center of the window is at the origin.
    # The `z` is from original vector analysis. It is the point on the vector
    # where the largest sampling sphere was (vector_analysed[0]).
    new_z = vector_analysed[0]
    # Translate the whole molecule to shift window's center to origin.
    coordinates = coordinates - np.array([[0, 0, new_z]] *
                                         coordinates.shape[0])
    # !!!Here the window center (xy and z) optimisation take place!!!
    window_com = np.array([0, 0, 0], dtype=float)
    # The lb_z parameter is 'lower bound equal to z' which means,
    # that we set the lower bound for the z optimisation to be equal
    # to the -new_z as in some cages it's the COM - void that is the
    # limiting diameter. But, no lower than new_z because we don't want to
    # move into the other direction.
    if lb_z:
        z_bounds[0] = -new_z
    window_diameter, _ = void_diameter(elements, coordinates, com=window_com)
    # SciPy minimisation on z coordinate.
    z_args = (window_com[0], window_com[1], elements, coordinates)
    z_optimisation = minimize(
        optimise_z, x0=window_com[2], args=z_args, bounds=[z_bounds])
    # Substitute the z coordinate for a minimised one.
    window_com[2] = z_optimisation.x[0]
    # SciPy brute optimisation on x and y coordinates in window plane.
    xy_args = (window_com[2], elements, coordinates)
    xy_bounds = ((-window_diameter / 2, window_diameter / 2),
                 (-window_diameter / 2, window_diameter / 2))
    xy_optimisation = brute(
        optimise_xy, xy_bounds, args=xy_args, full_output=True, finish=fmin)
    # Substitute the x and y coordinates for the optimised ones.
    window_com[0] = xy_optimisation[0][0]
    window_com[1] = xy_optimisation[0][1]
    # Additional SciPy minimisation on z coordinate. Added on 18 May 2017.
    # We can argue which aproach is best. Whether z opt and then xy opt
    # or like now z opt -> xy opt -> additional z opt etc. I have also tested
    # a loop of optimisations until some convergence and optimisation of
    # xyz coordinates at the same time by optimising these two optimisations.
    # In the end. I think this approach is best for cages.
    z_args = (window_com[0], window_com[1], elements, coordinates)
    # The z_bounds should be passed in kwargs.
    z_optimisation = minimize(
        optimise_z, x0=window_com[2], args=z_args, bounds=[z_bounds])
    # Substitute the z coordinate for a minimised one.
    window_com[2] = z_optimisation.x[0]
    # Calculate the new window diameter.
    window_diameter, _ = void_diameter(elements, coordinates, com=window_com)
    # To get the window true centre of mass we need to revere the rotation and
    # translation operations on the window com.
    # Reverse the translation by substracting the new_z.
    window_com[2] = window_com[2] + new_z
    angle_2_1 = -angle_2
    reverse_around_y = np.array([[np.cos(angle_2_1), 0, np.sin(angle_2_1)],
                                 [0, 1, 0],
                                 [-np.sin(angle_2_1), 0, np.cos(angle_2_1)]])
    # Reversing the second rotation around Y-axis.
    window_com = np.dot(reverse_around_y, window_com)
    angle_1_1 = -angle_1
    reverse_around_z = np.array([[np.cos(angle_1_1), -np.sin(angle_1_1), 0],
                                 [np.sin(angle_1_1), np.cos(angle_1_1), 0],
                                 [0, 0, 1]])
    # Reversing the first rotation around Z-axis.
    window_com = np.dot(reverse_around_z, window_com)
    return (window_diameter, window_com)


def find_windows(elements,
                 coordinates,
                 ncpus=1,
                 mol_size=None,
                 adjust=1,
                 void_opt=True,
                 increment=1.0,
                 output='all',
                 **kwargs):
    """Return windows diameters and center of masses for a molecule."""
    # Copy the coordinates as will perform many opertaions on them
    coordinates = deepcopy(coordinates)
    # Center of our cartesian system is always at origin
    origin = np.array([0, 0, 0])
    # Initial center of mass to reverse translation at the end
    initial_com = center_of_mass(elements, coordinates)
    # Shift the cage to the origin using either the standard center of mass
    # or if void_opt flag is True, the optimised void center as center of mass
    if void_opt is True:
        # Normally the void is calculated from the COM of a molecule.
        # So, esentially the molecule's COM is the void's center.
        # To shift the molecule so that the center of the optimised void
        # is at the origin of the system and not the center of the not
        # optimised one, we need to adjust the shift. We also have to update
        # the initial com.
        com_adjust = initial_com - opt_void_diameter(elements, coordinates, **
                                                     kwargs)[2]
        initial_com = initial_com - com_adjust
        coordinates = shift_com(elements, coordinates, com_adjust=com_adjust)
    else:
        # Otherwise, we just shift the cage to the origin.
        coordinates = shift_com(elements, coordinates)
    # We create an array of vdw radii of elements.
    elements_vdw = np.array([[atomic_vdw_radius[x.upper()]] for x in elements])
    # We calculate maximum diameter of a molecule to determine the radius
    # of a sampling sphere neccessary to enclose the whole molecule.
    shpere_radius = max_dim(elements, coordinates)[2] / 2
    sphere_surface_area = 4 * np.pi * shpere_radius**2
    # Here we determine the number of sampling points necessary for a fine
    # sampling. Smaller molecules require more finner density of sampling
    # points on the sampling sphere's surface, whereas largen require less.
    # This formula was created so that larger molecule do not take much longer
    # to analyse, as number_sampling_points*length_of_sampling_vectors
    # results in quadratic increase of sampling time. The 250 factor was
    # specificly determined to produce close to 1 sampling point /Angstrom^2
    # for a sphere of radius ~ 24 Angstrom. We can adjust how fine is the
    # sampling by changing the adjust factor.
    number_of_points = int(np.log10(sphere_surface_area) * 250 * adjust)
    points_per_1A_surface = number_of_points / sphere_surface_area
    # Here I use code by Alexandre Devert for spreading points on a sphere:
    # http://blog.marmakoide.org/?p=1
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(number_of_points)
    z = np.linspace(1 - 1.0 / number_of_points, 1.0 / number_of_points - 1.0,
                    number_of_points)
    radius = np.sqrt(1 - z * z)
    points = np.zeros((number_of_points, 3))
    points[:, 0] = radius * np.cos(theta) * shpere_radius
    points[:, 1] = radius * np.sin(theta) * shpere_radius
    points[:, 2] = z * shpere_radius
    # Here we will compute the eps parameter for the sklearn.cluster.DBSCAN
    # (3-dimensional spatial clustering algorithm) which is the mean distance
    # to the closest point of all points.
    values = []
    tree = KDTree(points)
    for i in points:
        dist, ind = tree.query(i.reshape(1, -1), k=10)
        values.append(dist[0][1])
        values.append(dist[0][2])
        values.append(dist[0][3])
    mean_closest_distance = np.mean(values)
    # The best eps is parametrized when adding the mean distance and it's root.
    eps = mean_closest_distance + mean_closest_distance**0.5
    # Here we either run the sampling points vectors analysis in serial
    # or parallel. The vectors that go through molecular voids return
    # as analysed list with the increment at vector's path with largest
    # included sphere, coordinates for this narrow channel point. vectors
    # that find molecule on theirs path are return as NoneType object.
    if ncpus == 1:
        results = [
            vector_analysis(
                point, coordinates, elements_vdw, increment=increment)
            for point in points
        ]
        results = [x for x in results if x is not None]
        dataset = np.array([x[5:8] for x in results])
    # Parralel analysis on user's defined number of CPUs.
    if ncpus > 1:
        pool = Pool(processes=ncpus)
        parallel = [
            pool.apply_async(
                vector_analysis,
                args=(
                    point,
                    coordinates,
                    elements_vdw, ),
                kwds={'increment': increment}) for point in points
        ]
        results = [p.get() for p in parallel if p.get() is not None]
        pool.terminate()
        # Dataset is an array of sampling points coordinates.
        dataset = np.array([x[5:8] for x in results])
    # If not a single vector was returned from the analysis it mean that
    # no molecular channels (what we call windows here) connects the
    # molecule's interior with the surroungsings (exterior space).
    # The number of windows in that case equals zero and zero is returned.
    # Otherwise we continue our search for windows.
    if len(results) == 0:
        return None
    else:
        # Perfomr DBSCAN to cluster the sampling points vectors.
        # the n_jobs will be developed later.
        # db = DBSCAN(eps=eps, n_jobs=_ncpus).fit(dataset)
        db = DBSCAN(eps=eps).fit(dataset)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = set(db.labels_)
        # Assing cluster label to a sampling point.
        clusters = [[i, j] for i, j in zip(results, db.labels_)]
        clustered_results = {label: [] for label in labels}
        # Create a dictionary of clusters with points listed.
        [clustered_results[i[1]].append(i[0]) for i in clusters]
        # No for the sampling point vector in each cluster that had
        # the widest channel's 'neck' is assumed to pass the closest
        # to the window's center and therefore will be passed to
        # window analysis function.
        # We also pass user defined settings for window analysis.
        # Again either in serlia or in parallel.
        # Noisy points get a cluster label -1, therefore we have to exclude it.
        if ncpus == 1:
            window_results = [
                window_analysis(
                    np.array(clustered_results[cluster]), elements,
                    coordinates, elements_vdw, **kwargs)
                for cluster in clustered_results if cluster != -1
            ]
        elif ncpus > 1:
            pool = Pool(processes=ncpus)
            parallel = [
                pool.apply_async(
                    window_analysis,
                    args=(np.array(clustered_results[cluster]), elements,
                          coordinates, elements_vdw),
                    kwds=kwargs) for cluster in clustered_results
                if cluster != -1
            ]
            window_results = [p.get() for p in parallel if p.get() is not None]
            pool.terminate()
    # The function returns two numpy arrays, one with windows diameters
    # in Angstrom, second with corresponding windows center's coordinates
    windows = np.array([result[0] for result in window_results
                        if result is not None])
    windows_coms = np.array(
        [np.add(result[1], initial_com) for result in window_results
         if result is not None])
    # Safety measures, if one of the windows is None or negative a warning
    # should be raised.
    for result in window_results:
        if result is None:
            msg_ = " ".join(
                ['Warning. One of the analysed windows has',
                 'returned as None. See manual.']
            )
            # print(msg_)
        elif result[0] < 0:
            msg_ = " ".join(
                ['Warning. One of the analysed windows has a vdW',
                 'corrected diameter smaller than 0. See manual.']
            )
            # print(msg_)
    if output == 'all':
        return (windows, windows_coms)
    if output == 'windows':
        return windows


def vector_analysis_reversed(vector, coordinates, elements_vdw, shpere_radius,
                             increment=0.1):
    """Analyse a sampling vector's path for avarge diamatere of a molecule."""
    # Calculate number of chunks if vector length is divided by increment.
    chunks = int(np.linalg.norm(vector) // increment)
    # Create a single chunk.
    chunk = vector / chunks
    # Calculate set of points on vector's path every increment.
    vector_pathway = [chunk * i for i in range(chunks + 1)]
    reversed_vector_pathway = np.array(vector_pathway[::-1])
    analysed_vector = np.array([
        np.amin(
            euclidean_distances(coordinates, i.reshape(1, -1)) - elements_vdw)
        for i in reversed_vector_pathway
    ])
    if all(i > 0 for i in analysed_vector):
        return None
    else:
        count = -1
        for i in analysed_vector:
            if i > 0:
                pass
            else:
                break
            count += 1
        dist_origin = np.linalg.norm(reversed_vector_pathway[count])
        dist_origin_corrected = dist_origin - analysed_vector[count]
        return [dist_origin_corrected, reversed_vector_pathway[count]]


def find_avarage_diameter(elements, coordinates, adjust=1, increment=0.1,
                          **kwargs):
    """Return avarage diameter for a molecule."""
    # Copy the coordinates as will perform many opertaions on them
    coordinates = deepcopy(coordinates)
    # Center of our cartesian system is always at origin
    origin = np.array([0, 0, 0])
    # Initial center of mass to reverse translation at the end
    initial_com = center_of_mass(elements, coordinates)
    # We just shift the cage to the origin.
    coordinates = shift_com(elements, coordinates)
    # We create an array of vdw radii of elements.
    elements_vdw = np.array([[atomic_vdw_radius[x.upper()]] for x in elements])
    # We calculate maximum diameter of a molecule to determine the radius
    # of a sampling sphere neccessary to enclose the whole molecule.
    shpere_radius = max_dim(elements, coordinates)[2]
    sphere_surface_area = 4 * np.pi * shpere_radius**2
    # Here we determine the number of sampling points necessary for a fine
    # sampling. Smaller molecules require more finner density of sampling
    # points on the sampling sphere's surface, whereas largen require less.
    # This formula was created so that larger molecule do not take much longer
    # to analyse, as number_sampling_points*length_of_sampling_vectors
    # results in quadratic increase of sampling time. The 250 factor was
    # specificly determined to produce close to 1 sampling point /Angstrom^2
    # for a sphere of radius ~ 24 Angstrom. We can adjust how fine is the
    # sampling by changing the adjust factor.
    number_of_points = int(np.log10(sphere_surface_area) * 250 * adjust)
    points_per_1A_surface = number_of_points / sphere_surface_area
    # Here I use code by Alexandre Devert for spreading points on a sphere:
    # http://blog.marmakoide.org/?p=1
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(number_of_points)
    z = np.linspace(1 - 1.0 / number_of_points, 1.0 / number_of_points - 1.0,
                    number_of_points)
    radius = np.sqrt(1 - z * z)
    points = np.zeros((number_of_points, 3))
    points[:, 0] = radius * np.cos(theta) * shpere_radius
    points[:, 1] = radius * np.sin(theta) * shpere_radius
    points[:, 2] = z * shpere_radius
    # Here we will compute the eps parameter for the sklearn.cluster.DBSCAN
    # (3-dimensional spatial clustering algorithm) which is the mean distance
    # to the closest point of all points.
    values = []
    tree = KDTree(points)
    for i in points:
        dist, ind = tree.query(i.reshape(1, -1), k=10)
        values.append(dist[0][1])
        values.append(dist[0][2])
        values.append(dist[0][3])
    mean_closest_distance = np.mean(values)
    # The best eps is parametrized when adding the mean distance and it's root.
    eps = mean_closest_distance + mean_closest_distance**0.5
    # Here we either run the sampling points vectors analysis in serial
    # or parallel. The vectors that go through molecular voids return
    # as analysed list with the increment at vector's path with largest
    # included sphere, coordinates for this narrow channel point. vectors
    # that find molecule on theirs path are return as NoneType object.
    results = [
        vector_analysis_reversed(
            point, coordinates, elements_vdw, shpere_radius,
            increment=increment)
        for point in points
    ]
    results_cleaned = [x[0] for x in results if x is not None]
    avarage_molecule_diameter = np.mean(results_cleaned)
    points_density = []
    for i in np.arange(1, int(shpere_radius)+1, increment):
        surface = 4 * np.pi * i**2
        density = number_of_points/surface
        points_density.append([i, density])
    normalised = []
    for i in points_density:
        normalised.append([i[0], points_density[0][1]/i[1]])
    weighted_avarage = [[] for x in range(len(normalised))]
    for i in results_cleaned:
        for j in range(len(normalised)-1):
            if normalised[j][0] < i <= normalised[j+1][0]:
                weighted_avarage[j].append(i)
    average_1 = 0
    average_2 = 0
    for i, j in zip(weighted_avarage, normalised):
        if i:
            average_1 += np.mean(i) * j[1]
            average_2 += j[1]
    average = average_1 / average_2
    return average * 2


def vector_analysis_pore_shape(vector, coordinates, elements_vdw,
                               increment=1.0):
    """Analyse a sampling vector's path for pore shape determination."""
    # Calculate number of chunks if vector length is divided by increment.
    chunks = int(np.linalg.norm(vector) // increment)
    # Create a single chunk.
    chunk = vector / chunks
    # Calculate set of points on vector's path every increment.
    vector_pathway = np.array([chunk * i for i in range(chunks + 1)])
    analysed_vector = np.array([
        np.amin(
            euclidean_distances(coordinates, i.reshape(1, -1)) - elements_vdw)
        for i in vector_pathway
    ])
    if all(i > 0 for i in analysed_vector):
        return None
    else:
        count = -1
        for i in analysed_vector:
            if i > 0:
                pass
            else:
                break
            count += 1
        return vector_pathway[count]


def calculate_pore_shape(elements, coordinates, adjust=1, increment=0.1,
                         **kwargs):
    """Return avarage diameter for a molecule."""
    # Copy the coordinates as will perform many opertaions on them
    coordinates = deepcopy(coordinates)
    # Center of our cartesian system is always at origin
    origin = np.array([0, 0, 0])
    # Initial center of mass to reverse translation at the end
    initial_com = center_of_mass(elements, coordinates)
    # We just shift the cage to the origin.
    coordinates = shift_com(elements, coordinates)
    # We create an array of vdw radii of elements.
    elements_vdw = np.array([[atomic_vdw_radius[x.upper()]] for x in elements])
    # We calculate maximum diameter of a molecule to determine the radius
    # of a sampling sphere neccessary to enclose the whole molecule.
    shpere_radius = max_dim(elements, coordinates)[2]/2
    sphere_surface_area = 4 * np.pi * shpere_radius**2
    # Here we determine the number of sampling points necessary for a fine
    # sampling. Smaller molecules require more finner density of sampling
    # points on the sampling sphere's surface, whereas largen require less.
    # This formula was created so that larger molecule do not take much longer
    # to analyse, as number_sampling_points*length_of_sampling_vectors
    # results in quadratic increase of sampling time. The 250 factor was
    # specificly determined to produce close to 1 sampling point /Angstrom^2
    # for a sphere of radius ~ 24 Angstrom. We can adjust how fine is the
    # sampling by changing the adjust factor.
    number_of_points = int(np.log10(sphere_surface_area) * 250 * adjust)
    points_per_1A_surface = number_of_points / sphere_surface_area
    # Here I use code by Alexandre Devert for spreading points on a sphere:
    # http://blog.marmakoide.org/?p=1
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(number_of_points)
    z = np.linspace(1 - 1.0 / number_of_points, 1.0 / number_of_points - 1.0,
                    number_of_points)
    radius = np.sqrt(1 - z * z)
    points = np.zeros((number_of_points, 3))
    points[:, 0] = radius * np.cos(theta) * shpere_radius
    points[:, 1] = radius * np.sin(theta) * shpere_radius
    points[:, 2] = z * shpere_radius
    # Here we will compute the eps parameter for the sklearn.cluster.DBSCAN
    # (3-dimensional spatial clustering algorithm) which is the mean distance
    # to the closest point of all points.
    values = []
    tree = KDTree(points)
    for i in points:
        dist, ind = tree.query(i.reshape(1, -1), k=10)
        values.append(dist[0][1])
        values.append(dist[0][2])
        values.append(dist[0][3])
    mean_closest_distance = np.mean(values)
    # The best eps is parametrized when adding the mean distance and it's root.
    eps = mean_closest_distance + mean_closest_distance**0.5
    # Here we either run the sampling points vectors analysis in serial
    # or parallel. The vectors that go through molecular voids return
    # as analysed list with the increment at vector's path with largest
    # included sphere, coordinates for this narrow channel point. vectors
    # that find molecule on theirs path are return as NoneType object.
    results = [
        vector_analysis_pore_shape(
            point, coordinates, elements_vdw, increment=increment)
        for point in points
    ]
    results_cleaned = [x for x in results if x is not None]
    ele = np.array(['He'] * len(results_cleaned))
    coor = np.array(results_cleaned)
    return ele, coor
