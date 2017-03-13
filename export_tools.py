"""
Defines functions only useful when using MMEA as a library.

The functions here are not for use when running the GA.

"""

import os
import glob
import networkx as nx
import rdkit.Chem as chem
from functools import partial
from multiprocessing import Pool
import itertools as it

from .ga import Population
from .molecular import (StructUnit, functional_groups)
from .molecular.optimization import *


def fg_prune(ifolder, fg, fg_num, ext):
    """
    Deletes molecules without a given functional group from a folder.

    Parameters
    ----------
    ifolder : str
        The full path of the folder holding .mol files. Any
        molecules without functional group `fg` `fg_num` amount of times
        are deleted.

    fg : str
        The name of the functional group which the copied molecules must
        possess. The name must correspond to one of the name of a
        functional group defined within
        ``FGInfo.functional_groups_list``.

    fg_num : int
        The number of functional groups of type `fg` which the molecule
        must have in order to not be deleted.

    ext : str
        The extension of the molecular structure files in `ifolder`.
        All other files are skipped.

    Returns
    -------
    None : NoneType

    """

    # Go through all the files in `ifolder` if they do not have the
    # file extension `ext` go to the next file.
    for file_name in os.listdir(ifolder):
        if not file_name.endswith(ext):
            continue

        path = os.path.join(ifolder, file_name)
        # Make a StructUnit object and substitute the functional group
        # of type `fg`.
        mol = StructUnit(path)
        mol.untag_atoms()
        mol.func_grp = next((x for x in functional_groups if
                             x.name == fg), None)

        mol.tag_atoms()

        # Check that the correct number is present.
        if len(mol.functional_group_atoms()) != fg_num:
            print('Deleting {}.'.format(path))
            os.remove(path)


def fg_distance_prune(folder, fg, ext):
    """
    Deletes molecules with functional groups seperated by 1 atom.

    Parameters
    ----------
    folder : str
        The full path of the folder which holdes the molecules. The
        files are removed from this folder.

    fg : str
        The name of the functional group.

    ext : str
        The file extension of the structure files in `folder`. All
        other files are skipped.

    """

    # Go through all the files in `ifolder` if they do not have the
    # file extension `ext` go to the next file.
    for file_name in os.listdir(folder):
        path = os.path.join(folder, file_name)
        if not file_name.endswith(ext):
            continue

        # Make a StructUnit object and substitute the functional group
        # of type `fg`.
        mol = StructUnit(path, fg)

        # Make a mathematical graph of the molecule. Useful for finding
        # the separation between nodes (atoms).
        g = mol.graph()
        # Each bonder atom can act as either a start or end node on a
        # graph. Find the seperations between such nodes. If the
        # separation is 3 this means there is only one nodes between
        # the start and end nodes. As a result the functional groups
        # are separated by 1 atom and should be deleted.
        for start, end in it.combinations(mol.bonder_ids, 2):
            if nx.shortest_path_length(g, start, end) < 3:
                print('Removing {}.'.format(path))
                os.remove(path)
                break


def optimize_folder(folder, macromodel_path,
                    timeout=False, temp=300,
                    confs=500, eq_time=50, sim_time=500):
    """
    Optimizes all molecules found in `folder`.

    Parameters
    ----------
    folder : str
        The full path of the folder which holds .mol file for
        optimization.

    macromodel_path : str
        The full path to the Schrodinger home directory.

    timeout : float (default=False)
        The seconds before the optimization is cancalled. If ``False``
        it will run until completion.

    temp : float (default=300)
        The temperature (K) of the MD conformer search.

    confs : float (default=500)
        The number of conformers sampled during the MD conformer search.

    eq_time : float (default=50)
        The equilibration time before the MD conformer search.

    sim_time : float (default=500)
        The simulation time of the MD conformer search.

    Returns
    -------
    None : NoneType

    """

    # First make a list holding all macromolecule objects to be
    # optimized. Because the objects need to be initialized from a .mol
    # file a StructUnit instance not a MacroMolecule instance is used.
    # minimal = True, because it is all that is needed to run
    # optimziations, plus you want to avoid doing functional group
    # substitutions.

    # Get the names of all the .mol files.
    names = [os.path.join(folder, file_name) for file_name in
             os.listdir(folder) if file_name.endswith(".mol")]
    # Make the StructUnit instances from the .mol files.
    macro_mols = [StructUnit(file_path) for file_path in  names]

    # Run the opt.
    md_opt = partial(macromodel_md_opt, macromodel_path=macromodel_path,
        timeout=timeout, temp=temp, confs=confs,
        eq_time=eq_time, sim_time=sim_time)

    with Pool() as p:
        p.map(md_opt, macro_mols)


def redump_pop(*folders, ofolder, cls):
    """
    Collects all the MacroMolecule .dmp files and creates a dump file.

    The function creates a population of all the .dmp files in `folders`
    and creates a dump file of that population. The dump file gets
    placed in `ofolder`.

    Parameters
    ----------
    *folders : str
        The full paths of the folders holding the MacroMolecule .dmp
        files.

    ofolder : str
        The folder in which the dump file should be placed.

    cls : MacroMolecule child class
        The class to which the .dmp files belong.

    Modifies
    --------
    ofolder/pop_dump
        Creates a file holding the population dump at this location.

    Returns
    -------
    None : NoneType

    """

    # Collect the paths of all the .dmp files.
    cage_paths = []
    for folder in folders:
        s = os.path.join(folder, '*.dmp')
        cage_paths.extend(glob.glob(s))

    # Make a population holding all the MacroMolecule objects and dump
    # it.
    pop = Population(*(cls.load(x) for x in cage_paths))
    pop.dump(os.path.join(ofolder, 'pop_dump'))


def substruct_prune(folder, ext, substruct):
    """
    Deletes molecules which contain the substructure `substruct`.

    Parameters
    ----------
    folder : str
        The full path of the folder from which the files are checked for
        the substructure and deleted.

    ext : str
        The extension of the structure files in `folder`. All other
        files are skipped.

    substruct : str
        The SMARTS string of the substructure, which if present in a
        molecule causes it to be deleted from `folder`.

    """

    # Create a rdkit molecule of the substructure.
    substruct_mol = chem.MolFromSmarts(substruct)
    # Go through all the files in `ifolder` if they do not have the
    # file extension `ext` go to the next file.
    for file_name in os.listdir(folder):
        if not file_name.endswith(ext):
            continue

        # Make a molecule from the file.
        path = os.path.join(folder, file_name)
        mol = StructUnit(path, minimal=True)
        # Check for substruct and delete as appropriate.
        if mol.mol.HasSubstructMatch(substruct_mol):
            print('Removing {}.'.format(path))
            os.remove(path)
