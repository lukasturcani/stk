"""
Defines classes which describe molecules.

There are a couple of major classes defined here. The most important
ones are :class:`StructUnit` and :class:`MacroMolecule`. Here is an
overview of their role and interaction. This is followed by a
step-by-step guide to macromolecular assembly.

:class:`StructUnit` represents the monomers that make up
macromolecules. These are commonly refered to as ``building blocks`` in
the documentation. :class:`StructUnit` holds information concerning
only a single building block molecule. For example, the number of atoms
and bonds a building block may have. It also has information about the
functoinal groups present in the building block molecule (see
:class:`.FGInfo`). The class also allows manipulation of the lone
building block molecule, such as rotations and translations.

The :class:`StructUnit` should be inherited as necessary. For example,
:class:`StructUnit2` adds manipulations relavant to molecules with 2
functional groups. :class:`StructUnit3` adds manipulations relavant to
molecules with 3 or more functional groups. If you have a monomer which
needs specific information or manipulations, give it its own class.

The :class:`MacroMolecule` class represents an assembled macromolecule.
It requires at least 2 basic pieces of information: which monomers are
used to assemble the macromolecule and what the topology/structure
of the macromolecule is.

:attr:`MacroMolecule.building_blocks` holds a :class:`list` of
:class:`StructUnit` instances. These represent the monomers which
make up the macromolecule. Only one :class:`StructUnit` instance per
monomer type is held. So if 4 of one type of monomer and 2 of another
type of monomer form a macromolecule, only 2 :class:`StructUnit`
instances are in :attr:`MacroMolecule.building_blocks`.

:attr:`MacroMolecule.topology` holds a :class:`Topology` instance. This
instance is responsible for assembling the macromolecule from the
building blocks. The building should happen in
:meth:`MacroMolecule.__init__` via :meth:`.Topology.build`. The
:meth:`~.Topology.build` method places the assembled macromoleclue in
:attr:`MacroMolecule.mol` as an ``rdkit`` molecule.

:class:`StructUnit` labels atoms in the functional groups of the
building blocks as either ``'bonder'`` or ``'del'`` (see its
documentation). This tells the :class:`.Topology` instance which atoms
form bonds and which are removed during assembly.

.. _`macromolecular assembly`:

A more detailed description of macromolecular assembly.
-------------------------------------------------------

This is a step-by-step guide of how macromolecular assembly is carried
out and what the classes do.

First you create :class:`StructUnit` instances of the building blocks
which make up the macromolecule:

.. code-block:: python

    bb = StructUnit('/path/to/struct/file.mol2', 'amine')

The :class:`StructUnit` instances are initialized using paths to
molecular structure files. (Initializing a :class:`StructUnit`
automatically completes steps 1 to 4.)

    1. Place an ``rdkit`` instance of the molecule into
       :attr:`StructUnit.mol`, i.e.

       .. code-block:: python

           bb.mol  # <rdkit.Chem.rdchem.Mol at 0x7f961a8f1f80>

    2. Scan the path of the structure file for the name of a functional
       group. (Alternatively the name of a functional group can be
       supplied to the initializer).

       .. code-block:: python

           bb.func_group.name  # 'amine'

Which functional groups are recognized by ``stk``?

The module :mod:`.functional_groups` defines the class :class:`.FGInfo`
and a :class:`tuple` of instances of this class called
:data:`functional_groups`. If you put an :class:`.FGInfo` instance into
:data:`functional_groups`, the functional group will be recognized.

    3. Place the :class:`.FGInfo` instance of the functional group into
       :attr:`StructUnit.func_grp`.

    4. Using :class:`.FGInfo`, tag atoms in the building block as
       either ``'bonder'`` or ``'del'``. ``'bonder'`` signifies that
       the atoms form a bond during macromolecular assembly, while
       ``'del'`` means they are deleted. Also tag all atoms in the
       functional group with the tag ``'fg'`` which holds the name
       of the functional group.

    5. Give the :class:`StructUnit` and :class:`.Topology` instances to
       the macromolecule's initializer.

       .. code-block:: python

           macro_mol = MacroMolecule([bb1, bb2], Topology())

       Normally, :class:`MacroMolecule` and :class:`.Topology` will
       not be used directly. Instead, classes derived from these
       will be used. For example,

           .. code-block:: python

               macro_mol = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

    6. Run :meth:`.Topology.build` inside
       :meth:`MacroMolecule.__init__`.

    7. The details of :meth:`.Topology.build` will vary depending on
       the :class:`.Topology` class used. However, the basic structure
       is the same (steps 8 - 10).

    8. Use :meth:`.Topology.place_mols` to combine the ``rdkit``
       molecules of all the building blocks into a single ``rdkit``
       instance. The combined ``rdkit`` instance is placed into
       ``macro_mol.mol``. :meth:`.Topology.place_mols` also gives
       each functional group a unique id. Every atom in the
       assembled molecule that belongs to a functional group is
       given the tag ``'fg_id'``, whose value holds a unique id for
       each functional group.

    9. Use :meth:`.Topology.prepare` to run any additional operations
       before joining up the building blocks and deleting extra
       atoms, this method may do nothing.

    10. Use :meth:`.Topology.bonded_fgs` to get the ids of functional
        groups in the molecule to react using :func:`.react`.

    11. Run :meth:`.Topology.cleanup` to perform any final operations
        on the assembled molecule. Can be nothing.

After all this you should have a ``rdkit`` instance of the
macromolecule which should be placed into :attr:`MacroMolecule.mol`.

.. _`adding macromolecules`:

Extending stk: Adding new macromolecules.
-----------------------------------------

To add new macromolecules create a new class which inherits
:class:`MacroMolecule`.

If you're adding a new class of macromolecules, it quite likely you
want to add a new :class:`.Topology` class. See the
:mod:`.topologies.base` for guidance on adding these. The topology
class does the assembly of the macromolecule from the building blocks.

"""

import tempfile
import warnings
import logging
import json
import os
import numpy as np
import networkx as nx
import itertools as it
import math
import rdkit.Geometry.rdGeometry as rdkit_geo
import rdkit.Chem.AllChem as rdkit
from rdkit.Chem import rdMolTransforms

from rdkit import DataStructs
from glob import glob
from functools import total_ordering, partial
from scipy.spatial.distance import euclidean
from scipy.optimize import minimize
from sklearn.metrics.pairwise import euclidean_distances

from collections import Counter, defaultdict
from inspect import signature

from . import topologies
from .functional_groups import functional_groups, react
from .energy import Energy
import pywindow
from ..utilities import (flatten,
                         periodic_table,
                         normalize_vector,
                         rotation_matrix,
                         vector_theta,
                         mol_from_mae_file,
                         rotation_matrix_arbitrary_axis,
                         atom_vdw_radii,
                         Cell)


logger = logging.getLogger(__name__)
# Toggles caching when making molecules.
CACHE_SETTINGS = {'ON': True}


class Cached(type):
    """
    A metaclass for creating classes which create cached instances.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = dict()

    def __call__(self, *args, **kwargs):
        sig = signature(self.__init__)
        sig = sig.bind_partial(self, *args, **kwargs)
        sig.apply_defaults()
        sig = sig.arguments
        key = self.gen_key(sig['building_blocks'], sig['topology'])

        if key in self.cache and CACHE_SETTINGS['ON']:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            if CACHE_SETTINGS['ON']:
                self.cache[key] = obj
            return obj


class CachedStructUnit(type):
    """
    A metaclass for making :class:`StructUnit` create cached instances.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = dict()

    def __call__(self, *args, **kwargs):
        # Get the arguments given to the initializer as a dictionary
        # mapping argument name to argument value.
        sig = signature(self.__init__)
        sig = sig.bind_partial(self, *args, **kwargs)
        sig.apply_defaults()
        sig = sig.arguments

        _, ext = os.path.splitext(sig['file'])

        # Ensure a valid file type was provided.
        if ext not in self.init_funcs:
            raise TypeError(('Unable to initialize'
                             ' from "{}" files.').format(ext))

        mol = self.init_funcs[ext](sig['file'])

        # Get the name of the functional group provided to the
        # initializer or get it from the path.
        if sig['functional_group']:
            fg = sig['functional_group']
        else:
            fg = next((x.name for x in functional_groups if
                       x.name in sig['file']), None)

        key = self.gen_key(mol, fg)
        if key in self.cache and CACHE_SETTINGS['ON']:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            if CACHE_SETTINGS['ON']:
                self.cache[key] = obj
            return obj


class Molecule:
    """
    The most basic class representing molecules.

    This class defines the operations which any class
    describing molecules should inherit or may find useful. Examples of
    such are :class:`StructUnit` and :class:`MacroMolecule`. This class
    should not be used directly.

    Attributes
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        An ``rdkit`` molecule instance representing the molecule.

    inchi : :class:`str`
        The InChI of the molecule.

    energy : :class:`.Energy`
        Handles all things energy.

    atom_props : :class:`dict`
        Maps atom id to a :class:`dict` holding the ``rdkit`` atom
        properties of the atom. For example

        .. code-block:: python

            atom_props = {0: {'bonder': 1,
                              'fg': 'amine',
                              'fg_id': 1},
                          5: {'fg': 'amine',
                              'fg_id': 2,
                              'del': 1}}

    bonder_ids : :class:`list`
        Holds the id of bonder atoms in a nested :class:`list`

        .. code-block:: python

            bonder_ids = [[1, 3], [5, 6], [11, 14]]

        This means atoms with ids ``1`` and ``3`` are bonder atoms with
        ``fg_id`` of ``0``, atoms ``5`` and ``6`` are bonder atoms with
        ``fg_id`` of ``1`` and ``11`` and ``14`` are bonder atoms with
        ``fg_id`` of ``2``.

    optimized : :class:`bool`
        Indicates whether a :class:`Molecule` has been passed through
        an optimization function or not.

    name : :class:`str`
        A name which can be optionally given to the molecule for easy
        identification.

    note : :class:`str`
        A note or comment about the molecule. Purely optional but can
        be useful for labelling and debugging.

    """

    def __init__(self, name="", note=""):
        self.optimized = False
        self.energy = Energy(self)
        self.name = name
        self.note = note

        self.save_atom_props()

    def all_atom_coords(self, conformer=-1):
        """
        Yields the coordinates of atoms in :attr:`mol`.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Yields
        ------
        :class:`tuple`
            The yielded :class:`tuple` has the form

            .. code-block:: python

                (32, numpy.array([12, 34, 3]))

            Where the first element is the atom id and the second
            element is an array holding the coordinates of the atom.

        """

        # Get the conformer from the rdkit instance.
        conf = self.mol.GetConformer(conformer)

        # Go through all the atoms and ask the conformer to return
        # the position of each atom. This is done by supplying the
        # conformers `GetAtomPosition` method with the atom's id.
        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_position = conf.GetAtomPosition(atom_id)
            yield atom_id, np.array([*atom_position])

    def atom_coords(self, atom_id, conformer=-1):
        """
        Return coordinates of an atom.

        Parameters
        ----------
        atom_id : :class:`int`
            The id of the atom whose coordinates are desired.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`numpy.array`
            An array holding the x, y and z coordinates of the atom.

        """

        conf = self.mol.GetConformer(conformer)
        atom_position = conf.GetAtomPosition(atom_id)
        return np.array([*atom_position])

    def atom_distance(self, atom1_id, atom2_id, conformer=-1):
        """
        Return the distance between 2 atoms.

        Parameters
        ----------
        atom1_id : :class:`int`
            The id of the first atom.

        atom2_id : :class:`int`
            The id of the second atom.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`scipy.double`
            The distance between the first and second atoms.

        """

        # Get the atomic positions of each atom and use the scipy
        # function to calculate their distance in Euclidean space.
        atom1_coords = self.atom_coords(atom1_id, conformer)
        atom2_coords = self.atom_coords(atom2_id, conformer)
        return euclidean(atom1_coords, atom2_coords)

    def atom_symbol(self, atom_id):
        """
        Returns the symbol of the atom with id `atom_id`.

        Parameters
        ----------
        atom_id : :class:`int`
            The id number of the atom.

        Returns
        -------
        :class:`str`
            The atomic symbol of the atom.

        """

        return self.mol.GetAtomWithIdx(atom_id).GetSymbol()

    def bonder_centroids(self, conformer=-1):
        """
        Calculates the centriod of bonder atoms in each fg.

        The centroids are yielded in order, with the centroid of
        the functional group with ``fg_id`` of ``0`` first.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`numpy.array`
            The bonder centroid of a functional group.

        """

        for atoms in self.bonder_ids:
            yield (sum(self.atom_coords(a, conformer) for a in atoms)
                   / len(atoms))

    def bonder_centroid(self, conformer=-1):
        """
        Returns the centroid of the bonder atoms.

        The calculation has two stages. First, the centroids of
        bonder atoms within the same functional groups are found.
        Second, the centroid of the centroids found in stage 1 is
        calculated and returned.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.array`
            An array holding the midpoint of the bonder atoms.

        """

        return sum(self.bonder_centroids(conformer)) / len(self.bonder_ids)

    def _cavity_size(self, origin, conformer):
        """
        Calculates diameter of the molecule from `origin`.

        The cavity is measured by finding the atom nearest to
        `origin`, correcting for van der Waals diameter and multiplying
        by -2.

        This function should not be used directly. Use
        :meth:`cavity_size` instead, which finds the optimal value of
        `origin` to use.

        Parameters
        ----------
        origin : :class:`numpy.array`
            Holds the x, y and z coordinate of the position from which
            the cavity is measured.

        conformer : :class:`int`
            The id of the conformer to be used.

        Returns
        -------
        :class:`float`
            The (negative) diameter of the molecules cavity.

        """

        atom_vdw = np.array([atom_vdw_radii[x.GetSymbol()] for x
                            in self.mol.GetAtoms()])

        distances = euclidean_distances(self.position_matrix(conformer).T,
                                        np.matrix(origin))
        distances = distances.flatten() - atom_vdw
        return -2*min(distances)

    def cavity_size(self, conformer=-1):
        """
        Calculates the diameter of the molecule's cavity.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`float`
            The diameter of the molecule's cavity in Angstroms.

        """

        # This function uses _cavity_size() to calculate the cavity
        # size. _cavity_size() finds the closest atom to `origin` to
        # get its value of the cavity.

        # What this function does is finds the value of `origin` which
        # causes _cavity_size() to calculate the largest possible
        # cavity.
        ref = self.center_of_mass(conformer)
        icavity = 0.5*self._cavity_size(ref, conformer)
        bounds = [(coord+icavity, coord-icavity) for coord in ref]
        cavity_origin = minimize(lambda x: self._cavity_size(x, conformer),
                                 x0=ref,
                                 bounds=bounds).x
        cavity = -self._cavity_size(cavity_origin, conformer)
        return 0 if cavity < 0 else cavity

    def center_of_mass(self, conformer=-1):
        """
        Returns the centre of mass of the molecule.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`numpy.array`
            An array holding the coordinates of the center of mass.

        References
        ----------
        https://en.wikipedia.org/wiki/Center_of_mass

        """

        center = np.array([0., 0., 0.])
        total_mass = 0.
        for atom_id, coord in self.all_atom_coords(conformer):
            mass = self.mol.GetAtomWithIdx(atom_id).GetMass()
            total_mass += mass
            center += mass*coord
        return np.divide(center, total_mass)

    def centroid(self, conformer=-1):
        """
        Returns the centroid of the molecule.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`numpy.array`
            A numpy array holding the position of the centroid.

        """

        centroid = sum(x for _, x in self.all_atom_coords(conformer))
        return np.divide(centroid, self.mol.GetNumAtoms())

    def dihedral_strain(self,
                        dihedral_SMARTS='',
                        target=180,
                        conformer=-1):
        """
        Returns the difference between the average dihedral and target.

        The differences is a returned as a percent.

        Parameters
        ----------
        dihedral_SMARTS : :class:`str`
            The SMARTS code for the dihedral of interest.

        target : :class:`float`
            Float representing the target value for the dihedral angle.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`float`
            The percent difference between the average dihedral in the
            molecule and the target value.

        """

        match = rdkit.MolFromSmarts(dihedral_SMARTS)
        atoms_dihedral = self.mol.GetSubstructMatches(match)

        dihedral_info = []
        if len(atoms_dihedral) > 0 and len(atoms_dihedral[0]) != 0:
            for atoms_group in atoms_dihedral:
                # Calculate the dihedral angle
                dihedral_value = rdMolTransforms.GetDihedralDeg(
                                    self.mol.GetConformer(conformer),
                                    atoms_group[0],
                                    atoms_group[1],
                                    atoms_group[2],
                                    atoms_group[3])
                # check that the dihedral is calculated in the right direction
                if abs(dihedral_value) > 90:
                    dihedral_value = abs(dihedral_value)
                else:
                    dihedral_value = 180 - abs(dihedral_value)

                dihedral_info.append(dihedral_value)

            # Calculate the average dihedral value
            avg_dihedral = np.mean([abs(x) for x in dihedral_info])
            # Calculate the relative diff with the target dihedral value
            diff = (abs(target - avg_dihedral) / target) * 100
        else:
            # If the molecule does not contain the bond, give 1% strain.
            diff = 1

        return diff

    def dump(self, path):
        """
        Writes a JSON :class:`dict` of the molecule to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file to which the JSON dict should be
            written.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            json.dump(self.json(), f, indent=4)

    def fg_centroid(self, fg_id, conformer=-1):
        """
        The centroid of bonder atoms in a functional group.

        Parameters
        ----------
        fg_id : :class:`int`
            The id of the functional group.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.array`
            The coordinates of a bonder centroid.

        Raises
        ------
        :class:`RuntimeError`
            If `fg_id` is not found on any atoms.

        """

        s = np.array([0., 0., 0.])
        c = 0
        for a in self.mol.GetAtoms():
            if (a.HasProp('fg_id') and
                a.GetIntProp('fg_id') == fg_id and
               a.HasProp('bonder')):
                s += self.atom_coords(a.GetIdx())
                c += 1
        if not c:
            raise RuntimeError(f'No fg_id of {fg_id}.')
        return s / c

    def fg_distance(self, fg1, fg2, conformer=-1):
        """
        The distance between the bonder centroids of two fgs.

        Parameters
        ----------
        fg1 : :class:`int`
            The ``fg_id`` of the first fg.

        fg2 : :class:`int`
            The ``fg_id`` of the second fg.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The distance between `fg1` and `fg2`.

        """

        c1 = self.fg_centroid(fg1, conformer)
        c2 = self.fg_centroid(fg2, conformer)
        return euclidean(c1, c2)

    @classmethod
    def from_dict(self, json_dict, optimized=True, load_names=True):
        """
        Creates a :class:`Molecule` from a JSON :class:`dict`.

        The :class:`Molecule` returned has the class specified in
        `json_dict`, not :class:`Molecule`.

        Parameters
        ----------
        json_dict : :class:`dict`
            A :class:`dict` holding the JSON representation of a
            molecule.

        optimized : :class:`bool`, optional
            The value passed to :attr:`Molecule.optimized`

        load_names : :class:`bool`, optional
            If ``True`` then the ``name`` key stored in `json_dict`
            is loaded.

        Returns
        -------
        :class:`Molecule`
            The molecule represented by `json_dict`.

        """
        # Get the class of the object.
        c = globals()[json_dict['class']]
        json_dict['optimized'] = optimized
        json_dict['load_names'] = load_names
        return c._json_init(json_dict)

    def graph(self):
        """
        Returns a mathematical graph representing the molecule.

        Returns
        -------
        :class:`networkx.Graph`
            A graph where the nodes are the ids of the atoms and the
            edges are the bonds.

        """

        # Create a graph instance and add the atom ids as nodes. Use
        # the atom ids from each end of a bond to define edges. Do this
        # for all bonds to account for all edges.

        graph = nx.Graph()

        for atom in self.mol.GetAtoms():
            graph.add_node(atom.GetIdx())

        for bond in self.mol.GetBonds():
            graph.add_edge(bond.GetBeginAtomIdx(),
                           bond.GetEndAtomIdx())

        return graph

    @property
    def inchi(self):
        """
        Returns the InChI of the molecule.

        Returns
        -------
        :class:`str`
            The InChI of the molecule.

        """

        self.update_stereochemistry()
        return rdkit.MolToInchi(self.mol)

    @classmethod
    def load(cls, path, optimized=True, load_names=True):
        """
        Creates a :class:`Molecule` from a JSON file.

        The returned :class:`Molecule` has the class specified in the
        JSON file, not :class:`Molecule`.

        Parameters
        ----------
        path : :class:`str`
            The full path holding a JSON representation to a molecule.

        optimized : :class:`bool`, optional
            The value passed to :attr:`Molecule.optimized`.

        load_names : :class:`bool`, optional
            If ``True`` then the ``name`` key stored in the JSON file
            is loaded.

        Returns
        -------
        :class:`Molecule`
            The molecule held in `path`.

        """

        with open(path, 'r') as f:
            json_dict = json.load(f)

        return cls.from_dict(json_dict, optimized, load_names)

    def max_diameter(self, conformer=-1):
        """
        Returns the largest distance between 2 atoms in the molecule.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`tuple` of form (float, int, int)
            A :class:`tuple` of the form

            .. code-block:: python

                max_diameter = (312.3, 4, 54)

            Where the first element is the largest inter-atomic
            distance in the molecule. The next 2 elements are the ids
            of the involved atoms.

        """

        coords = self.position_matrix(conformer).T
        dist = euclidean_distances(coords, coords)
        vdw = np.matrix([atom_vdw_radii[self.atom_symbol(i)] for
                         i in range(self.mol.GetNumAtoms())])
        dist = dist + vdw + vdw.T
        maxid1, maxid2 = np.unravel_index(dist.argmax(), dist.shape)
        return dist[maxid1, maxid2], maxid1, maxid2

    def mdl_mol_block(self, conformer=-1):
        """
        Returns a V3000 mol block of the molecule.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`str`
            The V3000 mol block representing the molecule.

        """

        main_string = ("\n"
                       "     RDKit          3D\n"
                       "\n"
                       "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                       "M  V30 BEGIN CTAB\n"
                       "M  V30 COUNTS {0} {1} 0 0 0\n"
                       "M  V30 BEGIN ATOM\n"
                       "!!!ATOM!!!BLOCK!!!HERE!!!\n"
                       "M  V30 END ATOM\n"
                       "M  V30 BEGIN BOND\n"
                       "!!!BOND!!!BLOCK!!!HERE!!!\n"
                       "M  V30 END BOND\n"
                       "M  V30 END CTAB\n"
                       "M  END\n"
                       "\n"
                       "$$$$\n")

        # id atomic_symbol x y z
        atom_line = "M  V30 {} {} {:.4f} {:.4f} {:.4f} 0{}\n"
        atom_block = ""

        # id bond_order atom1 atom2
        bond_line = "M  V30 {} {} {} {}\n"
        bond_block = ""

        main_string = main_string.format(self.mol.GetNumAtoms(),
                                         self.mol.GetNumBonds())
        # Kekulize the mol, which means that each aromatic bond is
        # converted to a single or double. This is necessary because
        # .mol V3000 only supports integer bonds. However, this fails
        # sometimes on big molecules.
        try:
            rdkit.Kekulize(self.mol)
        except ValueError:
            pass

        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_sym = periodic_table[atom.GetAtomicNum()]
            charge = atom.GetFormalCharge()
            charge = '' if charge == 0 else f' CHG={charge}'
            x, y, z = self.atom_coords(atom_id, conformer)
            atom_block += atom_line.format(atom_id+1,
                                           atom_sym,
                                           x, y, z,
                                           charge)

        for bond in self.mol.GetBonds():
            bond_id = bond.GetIdx()
            atom1_id = bond.GetBeginAtomIdx() + 1
            atom2_id = bond.GetEndAtomIdx() + 1
            bond_order = int(bond.GetBondTypeAsDouble())
            # Ensure that no information was lost when converting
            # double to int.
            assert bond_order == bond.GetBondTypeAsDouble()
            bond_block += bond_line.format(bond_id, bond_order,
                                           atom1_id, atom2_id)

        main_string = main_string.replace(
                            "!!!ATOM!!!BLOCK!!!HERE!!!\n", atom_block)

        return main_string.replace(
                            "!!!BOND!!!BLOCK!!!HERE!!!\n", bond_block)

    def position_matrix(self, conformer=-1):
        """
        Returns the position of all atoms as a matrix.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`numpy.matrix`
            The matrix has a shape ``[3, n]``. Each column holds the x,
            y and z coordinates of an atom. The index of the column
            corresponds to the id of the atom in the molecule.

        """

        pos_array = np.array([])
        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            pos_vect = np.array([*self.atom_coords(atom_id, conformer)])
            pos_array = np.append(pos_array, pos_vect)

        return np.matrix(pos_array.reshape(-1, 3).T)

    def same(self, other):
        """
        Check if `other` has the same molecular structure.

        Parameters
        ----------
        other : :class:`Molecule`
            The :class:`Molecule` instance you are checking has
            the same structure.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the structures match.

        """

        return self.inchi == other.inchi

    def rotate(self, theta, axis, conformer=-1):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the molecular centroid.

        Parameters
        ----------
        theta : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.array`
            The axis about which the rotation happens.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Save the original position.
        og_position = self.centroid(conformer)
        # Move the centroid of the molecule to the origin, so that the
        # rotation occurs about this point.
        self.set_position([0, 0, 0], conformer)
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        new_pos_mat = np.dot(rot_mat, self.position_matrix(conformer))
        # Apply the rotation.
        self.set_position_from_matrix(new_pos_mat, conformer)
        # Return the centroid of the molecule to the origin position.
        self.set_position(og_position, conformer)

    def save_atom_props(self):
        """
        Updates :attr:`atom_props` and :attr:`bonder_ids`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.bonder_ids = bonder_ids = []
        self.atom_props = atom_props = defaultdict(dict)

        for atom in self.mol.GetAtoms():
            atomid = atom.GetIdx()
            for name, value in atom.GetPropsAsDict(False, False).items():
                atom_props[atomid][name] = value

                if name == 'bonder':
                    fg_id = atom.GetIntProp('fg_id')

                    # Make sure bonder_ids does not raise an index error.
                    if len(bonder_ids) < fg_id + 1:
                        diff = fg_id + 1 - len(bonder_ids)
                        bonder_ids.extend([] for i in range(diff))

                    bonder_ids[fg_id].append(atomid)

    def set_orientation(self, start, end, conformer=-1):
        """
        Rotates the molecule by a rotation from `start` to `end`.

        Given two direction vectors, `start` and `end`, this method
        applies the rotation required transform `start` to `end` on
        the molecule. The rotation occurs about the centroid of the
        molecule.

        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you
        can associate a geometric feature of the molecule with a
        vector, then the molecule can be rotated so that this vector is
        aligned with `end`. The defined vector can be virtually
        anything. This means that any geometric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Notes
        -----
        The difference between this method and
        :meth:`StructUnit._set_orientation2` is about which point the
        rotation occurs: centroid of the entire molecule versus
        centroid of the bonder atoms, respectively.

        Parameters
        ----------
        start : :class:`numpy.array`
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : :class:`numpy.array`
            This array holds the vector, onto which `start` is rotated.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)

        # Record the position of the molecule then translate the
        # centroid to the origin. This is so that the rotation occurs
        # about this point.
        og_center = self.centroid(conformer)
        self.set_position([0, 0, 0], conformer)

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix(conformer))

        # Set the positions of the molecule.
        self.set_position_from_matrix(new_pos_mat, conformer)
        self.set_position(og_center, conformer)

        return self.mol

    def set_position(self, position, conformer=-1):
        """
        Sets the centroid of the molecule to `position`.

        Parameters
        ----------
        position : :class:`numpy.array`
            This array holds the position on which the centroid of the
            molecule should be placed.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule with the centroid placed at
            `position`. This is the same instance as that in
            :attr:`Molecule.mol`.

        """

        conf_id = self.mol.GetConformer(conformer).GetId()

        # Get the original centroid.
        centroid = self.centroid(conf_id)
        # Find out how much it needs to shift to reach `position`.
        shift = position - centroid
        # Apply the shift and get the resulting rdkit conformer object.
        new_conf = self.shift(shift, conf_id).GetConformer()
        new_conf.SetId(conf_id)

        # Replace the old rkdit conformer with one where the centroid
        # is at `position`.
        self.mol.RemoveConformer(conf_id)
        self.mol.AddConformer(new_conf)

        return self.mol

    def set_position_from_matrix(self, pos_mat, conformer=-1):
        """
        Set atomic positions of the molecule to those in `pos_mat`.

        Parameters
        ----------
        pos_mat : :class:`numpy.array`
            The matrix holds the coordinates on which the atoms of the
            molecule should be placed.

            The shape of the matrix is ``[3, n]``. Each column of
            `pos_mat` represents the coordinates of a single atom. The
            1st column sets the coordinates of the atom with id of 0.
            The next column sets the coordinates of the atom with id 1,
            and so on.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        conf = self.mol.GetConformer(conformer)
        for i, coord_mat in enumerate(pos_mat.T):
            coord = rdkit_geo.Point3D(coord_mat.item(0),
                                      coord_mat.item(1),
                                      coord_mat.item(2))
            conf.SetAtomPosition(i, coord)

    def shift(self, shift, conformer=-1):
        """
        Shifts the coordinates of all atoms.

        This does not modify the molecule. A modified copy is returned.

        Parameters
        ----------
        shift : :class:`numpy.array`
            A numpy array holding the value of the shift along each
            axis.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            A copy of the molecule where the coordinates have been
            shifted by `shift`.

        """

        # The function does not modify the existing conformer, as a
        # result a new instance is created and used for modification.
        conf = rdkit.Conformer(self.mol.GetConformer(conformer))

        # For each atom, get the atomic positions from the conformer
        # and shift them. Create a new geometry instance from these new
        # coordinate values. The geometry instance is used by rdkit to
        # store the coordinates of atoms. Finally, set the conformers
        # atomic position to the values stored in this newly generated
        # geometry instance.
        for atom in self.mol.GetAtoms():

            # Remember the id of the atom you are currently using. It
            # is used to change the position of the correct atom at the
            # end of the loop.
            atom_id = atom.GetIdx()

            # `atom_position` in an instance holding in the x, y and z
            # coordinates of an atom in its 'x', 'y' and 'z'
            # attributes.
            atom_position = np.array(conf.GetAtomPosition(atom_id))

            # Inducing the shift.
            new_atom_position = atom_position + shift

            # Creating a new geometry instance.
            new_coords = rdkit_geo.Point3D(*new_atom_position)

            # Changes the position of the atom in the conformer to the
            # values stored in the new geometry instance.
            conf.SetAtomPosition(atom_id, new_coords)

        # Create a new copy of the rdkit molecule instance representing
        # the molecule - the original instance is not to be modified.
        new_mol = rdkit.Mol(self.mol)

        # The new rdkit molecule was copied from the one held in the
        # `mol` attribute, as result it has a copy of its conformer. To
        # prevent the rdkit molecule from holding multiple conformers
        # the `RemoveAllConformers` method is run first. The shifted
        # conformer is then given to the rdkit molecule, which is
        # returned.
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(conf)
        return new_mol

    def update_from_mae(self, path, conformer=-1):
        """
        Updates molecular structure to match an ``.mae`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mae`` file from which the structure
            should be updated.

        conformer : :class:`int`, optional
            The conformer to be updated.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()

        mol = mol_from_mae_file(path)
        conf = rdkit.Conformer(mol.GetConformer())
        conf.SetId(conformer)
        self.mol.RemoveConformer(conformer)
        self.mol.AddConformer(conf)

    def update_from_mol(self, path, conformer=-1):
        """
        Updates molecular structure to match an ``.mol`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mol`` file from which the structure
            should be updated.

        conformer : :class:`int`, optional
            The conformer to be updated.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()

        mol = rdkit.MolFromMolFile(path, sanitize=False, removeHs=False)
        conf = rdkit.Conformer(mol.GetConformer())
        conf.SetId(conformer)
        self.mol.RemoveConformer(conformer)
        self.mol.AddConformer(conf)

    def update_stereochemistry(self, conformer=-1):
        """
        Updates stereochemistry tags in :attr:`Molecule.mol`.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        for atom in self.mol.GetAtoms():
            atom.UpdatePropertyCache()
        rdkit.AssignAtomChiralTagsFromStructure(self.mol, conformer)
        rdkit.AssignStereochemistry(self.mol, True, True, True)

    def write(self, path, conformer=-1):
        """
        Writes a molecular structure file of the molecule.

        This bypasses the need to the writining functions in ``rdkit``.
        These have issues with macromolecules due to poor ring finding
        and sanitization issues.

        Parameters
        ----------
        path : :class:`str`
            The `path` to which the molecule should be written.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        write_funcs = {'.mol': self._write_mdl_mol_file,
                       '.sdf': self._write_mdl_mol_file,
                       '.pdb': self._write_pdb_file}

        _, ext = os.path.splitext(path)
        write_func = write_funcs[ext]
        write_func(path, conformer)

    def _write_mdl_mol_file(self, path, conformer=-1):
        """
        Writes a V3000 ``.mol`` file of the molecule

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            f.write(self.mdl_mol_block(conformer))

    def _write_pdb_file(self, path, conformer=-1):
        """
        Writes a ``.pdb`` file of the molecule

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # First write the file using rdkit.
        rdkit.MolToPDBFile(self.mol, path, conformer)

        # Edit the file because rkdit does poor atom labelling.
        new_content = ''
        with open(path, 'r') as pdb:
            for line in pdb:
                if 'HETATM' in line:
                    words = line.split()
                    lbl_word = words[2]
                    rpl_word = words[-1]
                    rpl_word += " "*(len(lbl_word)-len(rpl_word))
                    line = line.replace(lbl_word, rpl_word)

                new_content += line

        with open(path, 'w') as pdb:
            pdb.write(new_content)


class StructUnit(Molecule, metaclass=CachedStructUnit):
    """
    Represents the building blocks of macromolecules.

    The goal of this class is to conveniently store information about,
    and perform operations on, single instances of macromolecular
    building blocks.

    Attributes
    ----------
    init_funcs : :class:`dict`
        This dictionary holds the various functions which can be used
        to initialize ``rdkit`` molecules and pairs them with the
        appropriate file extension.

    file : :class:`str`
        The full path to the molecular structure file holding the
        molecule. The supported file formats are the keys in
        :attr:`init_funcs`. As long as a file with one of these
        extensions is provided, the correct initialization function
        will be used.

    func_grp : :class:`.FGInfo`
        The :class:`.FGInfo` instance holding information about the
        functional group which will react when the building block
        assembles to form macromolecules.

    """

    init_funcs = {'.mol': partial(rdkit.MolFromMolFile,
                                  sanitize=False, removeHs=False),

                  '.sdf': partial(rdkit.MolFromMolFile,
                                  sanitize=False, removeHs=False),

                  '.mol2': partial(rdkit.MolFromMol2File,
                                   sanitize=False, removeHs=False),

                  '.mae': mol_from_mae_file,

                  '.pdb': partial(rdkit.MolFromPDBFile,
                                  sanitize=False, removeHs=False)}

    def __init__(self, file, functional_group=None, name="", note=""):
        """
        Initializes a :class:`StructUnit` instance.

        Parameters
        ----------
        file : :class:`str`
            The full path of the molecular structure file holding the
            building block.

        functional_group : :class:`str`, optional
            The name of the functional group which is to have atoms
            tagged. If ``None``, a functional group name found in the
            path `file`  is used. If no functional group is provided
            to this parameter and the name of one is not present in
            `file`, no tagging is done.

        name : :class:`str`, optional
            A name which can be optionally given to the molcule for
            easy identification.

        note : :class:`str`, optional
            A note or comment about the molecule.

        """

        self.file = file
        _, ext = os.path.splitext(file)

        if ext not in self.init_funcs:
            raise TypeError(f'Unable to initialize from "{ext}" files.')

        self.mol = self.init_funcs[ext](file)
        # Update the property cache of each atom. This updates things
        # like valence.
        for atom in self.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # Define a generator which yields an ``FGInfo`` instance from
        # `functional_groups`. The yielded ``FGInfo``instance
        # represents the functional group of the molecule which will
        # undergo bond formation. The generator determines the
        # functional group of the molecule from the path of of the
        # structure file.

        # The database of precursors should be organized so that any
        # given structure file has the name of its functional group in
        # its path. Each file should have the name of only one
        # functional group in its path. If this is not the case, the
        # generator will return the functional group which appears
        # first in `functional_groups`.

        # Assign the FGInfo instance from `functional_groups` which
        # describes the functional group provided in `functional_group`
        # or is found in the path name.
        if functional_group:
            self.func_grp = next((x for x in functional_groups if
                                  x.name == functional_group), None)
        else:
            self.func_grp = next((x for x in functional_groups if
                                  x.name in file), None)

        # Calling this function labels the atoms in the rdkit molecule
        # as either atoms which form a bond during reactions or atoms
        # which get removed.
        if self.func_grp:
            self.tag_atoms()

        super().__init__(name, note)

    @classmethod
    def init_random(cls, db, fg=None, name="", note=""):
        """
        Picks a random file from `db` to initialize from.

        Parameters
        ----------
        db : :class:`str`
            A path to a database of molecular files.

        fg : :class:`str`, optional
            The name of a functional group which the molecules in `db`
            have. By default it is assumed the name is present in the
            path of the files.

        name : :class:`str`, optional
            The name to be given to the created molecule.

        note : :class:`str`, optional
            A note to be given to the created molecule.

        Returns
        -------
        :class:`StructUnit`
            A random molecule from `db`.

        None : :class:`NoneType`
            If no files in `db` could be initialized from.

        """

        files = glob(os.path.join(db, '*'))
        np.random.shuffle(files)

        for molfile in files:
            try:
                return cls(molfile, fg, name, note)

            except Exception:
                logger.error(
                    'Could not initialize {} from {}.'.format(
                                                cls.__name__, molfile))

    def all_bonder_distances(self, conformer=-1):
        """
        Yield distances between all pairs of bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`tuple`
            A :class:`tuple` of the form ``(3, 54, 12.54)``.
            The first two elements are the ``fg_id`` of the involved
            functional groups and the third element is the distance
            between them.

        """

        centroids = it.combinations(self.bonder_centroids(conformer), 2)
        ids = it.combinations(range(len(self.bonder_ids)), 2)
        for (id1, id2), (c1, c2) in zip(ids, centroids):
                yield id1, id2, euclidean(c1, c2)

    def bonder_direction_vectors(self, conformer=-1):
        """
        Yields the direction vectors between bonder atoms.

        First the centroids of bonder atoms within the same functional
        group are found. Then direction vectors between all pairs of
        these centroids are yielded. Each pair is only yielded once.

        The yielded vector is normalized.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`tuple`
            They yielded tuple has the form

            .. code-block:: python

                (3, 54, np.array([12.2, 43.3, 9.78]))

            The first two elements of the tuple represent the
            ``fg_id`` of the start and end fgs of the vector,
            respectively. The array is the direction vector running
            between the functional group positions.

        """

        centroids = it.combinations(self.bonder_centroids(conformer), 2)
        ids = it.combinations(range(len(self.bonder_ids)), 2)
        for (id1, id2), (c1, c2) in zip(ids, centroids):
            yield id2, id1, normalize_vector(c1-c2)

    def bonder_position_matrix(self, conformer=-1):
        """
        Returns a matrix holding the positions of bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.matrix`
            The matrix has the shape ``[3, n]``. Each column holds the
            x, y and z coordinates of a bonder centroid. The index of
            the column corresponds to the ``fg_id`` of the bonder
            centroid.

        """

        return np.matrix(list(self.bonder_centroids(conformer))).T

    def centroid_centroid_dir_vector(self, conformer=-1):
        """
        Returns the direction vector between the 2 molecular centroids.

        The first molecular centroid is the centroid of the entire
        molecule. The second molecular centroid is given by
        :meth:`bonder_centroid`.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.array`
            The normalized direction vector running from the centroid
            of the bonder atoms to the molecular centroid.

        """

        # If the bonder centroid and centroid are in the same position,
        # the centroid - centroid vector should be orthogonal to the
        # bonder direction vector.
        if np.allclose(self.centroid(conformer),
                       self.bonder_centroid(conformer), atol=1e-5):
            *_, bvec = next(self.bonder_direction_vectors(conformer))
            # Construct a secondary vector by finding the minimum
            # component of bvec and setting it to 0.
            vec2 = list(bvec)
            minc = min(vec2)
            vec2[vec2.index(min(vec2))] = 0 if abs(minc) >= 1e-5 else 1
            # Get a vector orthogonal to bvec and vec2.
            a = normalize_vector(np.cross(bvec, vec2))
            return a

        else:
            return normalize_vector(self.centroid(conformer) -
                                    self.bonder_centroid(conformer))

    def core(self):
        """
        Return the molecule with no H or functional group atoms.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The "core" of the molecule.

        """

        emol = rdkit.EditableMol(self.mol)
        for atom in reversed(self.mol.GetAtoms()):
            atomid = atom.GetIdx()
            if not self.is_core_atom(atomid):
                emol.RemoveAtom(atomid)
        return emol.GetMol()

    def functional_group_atoms(self):
        """
        Returns a container of atom ids of atoms in functional groups.

        Returns
        -------
        :class:`tuple`
            The form of the returned tuple is:

            .. code-block:: python

                ((1, 2, 3), (4, 5, 6), (7, 8, 9))

            This means that all atoms with ids ``1`` to ``9`` are in a
            functional group and that the atoms ``1``, ``2`` and ``3``
            all form one functional group together. So do ``4``, ``5``
            and ``6`` and so on.

        """

        # Generate a ``rdkit.Chem.rdchem.Mol`` instance which
        # represents the functional group of the molecule.
        func_grp_mol = rdkit.MolFromSmarts(self.func_grp.fg_smarts)

        # Do a substructure search on the the molecule in `mol` to find
        # which atoms match the functional group. Return the atom ids
        # of those atoms.
        return self.mol.GetSubstructMatches(func_grp_mol)

    def is_core_atom(self, atomid):
        """
        Returns ``True`` if atom is not H or part of a fg.

        Parameters
        ----------
        atomid : :class:`int`
            The id of the atom being queried.

        Returns
        -------
        :class:`bool`
            Indicates whether the atom with `atomid` is part of the
            core.

        """

        atom = self.mol.GetAtomWithIdx(atomid)
        if atom.GetAtomicNum() == 1:
            return False
        if atomid not in self.atom_props:
            return True

        return 'fg' not in self.atom_props[atomid]

    def json(self):
        """
        Returns a JSON representation of the molecule.

        The representation has the following form:

        .. code-block:: python

            {
                'class' : 'StructUnit',
                'mol_block' : '''A string holding the V3000 mol
                                 block of the molecule.''',
                'note' : 'This molecule is nice.',
                'name' : 'benzene'
            }

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        return {

            'class': self.__class__.__name__,
            'func_grp': (self.func_grp if
                         self.func_grp is None else
                         self.func_grp.name),
            'mol_block': self.mdl_mol_block(),
            'note': self.note,
            'name': self.name

        }

    @classmethod
    def _json_init(cls, json_dict):
        """
        Completes a JSON initialization.

        This function is not to be used. Use :meth:`Molecule.load`
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : :class:`dict`
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        with tempfile.NamedTemporaryFile('r+t', suffix='.mol') as f:
            f.write(json_dict['mol_block'])
            f.seek(0)
            obj = cls(f.name, json_dict['func_grp'],
                      (json_dict['name'] if
                       json_dict['load_names'] else ""),
                      json_dict['note'])
            obj.optimized = json_dict['optimized']
            return obj

    @staticmethod
    def gen_key(rdkit_mol, functional_group):
        """
        Generates the key used when caching the molecule.

        Parameters
        ----------
        rdkit_mol : :class:`rdkit.Chem.rdchem.Mol`
            An ``rdkit`` instance of the molecule.

        functional_group : :class:`str`
            The name of the functional group being used to make
            macromolecules.

        Returns
        -------
        :class:`tuple`
            The key used for caching the molecule. Has the form

            .. code-block:: python

                ('amine', 'InChIString')

        """

        return functional_group, rdkit.MolToInchi(rdkit_mol)

    def minimize_theta(self, v1, v2, axis, centroid, conformer=-1):
        """
        Rotates the molecule to minimize angle between `v1` and `v2`.

        The rotation is done about the vector `axis`.

        Parameters
        ----------
        v1 : :class:`numpy.array`
            The vector which is rotated.

        v2 : :class:`numpy.array`
            The vector which is stationary.

        axis : :class:`numpy.array`
            The vector about which the rotation happens.

        centroid : :class:`numpy.array`
            The position vector at the center of the rotation.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # If the vector being rotated is not finite exit. This is
        # probably due to a planar molecule.
        if not all(np.isfinite(x) for x in v1):
            return

        # Save the initial position and change the origin to
        # `centroid`.
        iposition = self.centroid(conformer)
        self.mol = self.shift(-centroid, conformer)

        # 1. First transform the problem.
        # 2. The rotation axis is set equal to the z-axis.
        # 3. Apply this transformation to all vectors in the problem.
        # 4. Take only the x and y components of `v1` and `v2`.
        # 5. Work out the angle between them.
        # 6. Apply that rotation along the original rotation axis.

        rotmat = rotation_matrix(axis, [0, 0, 1])
        tstart = np.dot(rotmat, v1)
        tstart = np.array([tstart[0], tstart[1], 0])

        # If the `tstart` vector is 0 after these transformations it
        # means that it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], atol=1e-8):
            self.set_position(iposition, conformer)
            return

        tend = np.dot(rotmat, v2)
        tend = np.array([tend[0], tend[1], 0])
        angle = vector_theta(tstart, tend)

        # Check in which direction the rotation should go.
        # This is done by applying the rotation in each direction and
        # seeing which one leads to a smaller theta.
        r1 = rotation_matrix_arbitrary_axis(angle, [0, 0, 1])
        t1 = vector_theta(np.dot(r1, tstart), tend)
        r2 = rotation_matrix_arbitrary_axis(-angle, [0, 0, 1])
        t2 = vector_theta(np.dot(r2, tstart), tend)

        if t2 < t1:
            angle *= -1

        rotmat = rotation_matrix_arbitrary_axis(angle, axis)
        posmat = np.dot(rotmat, self.position_matrix(conformer))
        self.set_position_from_matrix(posmat, conformer)
        self.set_position(iposition, conformer)

    @classmethod
    def rdkit_init(cls, mol, functional_group=None, name="", note=""):
        """
        Uses an ``rdkit`` molecule for initialization.

        Parameters
        ----------
        mol : :class:`rdkit.Chem.rdchem.Mol`
            An ``rdkit`` molecule used for initialization.

        Returns
        -------
        :class:`StructUnit`
            A :class:`StructUnit` of `mol`.

        """

        with tempfile.NamedTemporaryFile('r+t', suffix='.mol') as f:
            f.write(rdkit.MolToMolBlock(mol, forceV3000=True))
            f.seek(0)
            return cls(f.name, functional_group, name, note)

    def rotate2(self, theta, axis, conformer=-1):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the centroid of the bonder atoms.

        Parameters
        ----------
        theta : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.array`
            The axis about which rotation happens.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Save the origin position of the bonder atom centroid.
        og_position = self.bonder_centroid(conformer)
        # Change the position of the centroid of the bonder atoms to
        # the origin so that the rotation occurs about this point.
        self.set_bonder_centroid([0, 0, 0], conformer)
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation on the original atomic coordinates to get
        # the new ones.
        new_pos_mat = np.dot(rot_mat, self.position_matrix(conformer))
        # Set the atomic positions to the new coordinates.
        self.set_position_from_matrix(new_pos_mat, conformer)
        # Return the centroid to its original position.
        self.set_bonder_centroid(og_position, conformer)

    def set_bonder_centroid(self, position, conformer=-1):
        """
        Move the molecule so that the bonder centroid is on `position`.

        Parameters
        ----------
        position : :class:`numpy.array`
            An array holding the desired the position. It holds
            the x, y and z coordinates, respectively.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule after it has been shifted. The same
            instance as in :attr:`~Molecule.mol`.

        """

        conf_id = self.mol.GetConformer(conformer).GetId()

        center = self.bonder_centroid(conf_id)
        shift = position - center
        new_conf = self.shift(shift, conf_id).GetConformer()
        new_conf.SetId(conf_id)

        # Make sure the rkdit molecule has only one conformer.
        self.mol.RemoveConformer(conf_id)
        self.mol.AddConformer(new_conf)

        return self.mol

    def _set_orientation2(self, start, end, conformer):
        """
        Rotates the molecule by a rotation from `start` to `end`.

        Given two direction vectors, `start` and `end`, this method
        applies the rotation required transform `start` to `end` on
        the molecule. The rotation occurs about the centroid of the
        molecule.

        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you
        can associate a gemotric feature of the molecule with a vector,
        then the molecule can be roatated so that this vector is
        aligned with `end`. The defined vector can be virtually
        anything. This means that any geometric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Notes
        -----
        The difference between this method and
        :meth:`~Molecule.set_orientation` is about which point the
        rotation occurs: centroid of bonder atoms versus centroid of
        entire molecule, respectively.

        Parameters
        ----------
        start : :class:`numpy.array`
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : :class:`numpy.array`
            This array holds the vector, onto which `start` is rotated.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)

        # Record the position of the molecule then translate the bonder
        # atom centroid to the origin. This is so that the rotation
        # occurs about this point.
        og_center = self.bonder_centroid(conformer)
        self.set_bonder_centroid(np.array([0, 0, 0]), conformer)

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix(conformer))

        # Set the positions in the rdkit molecule.
        self.set_position_from_matrix(new_pos_mat, conformer)
        self.set_bonder_centroid(og_center, conformer)

        return self.mol

    def similar_molecules(self, mols):
        """
        Returns molecules from `mols` ordered by similarity.

        The most similar molecule is at index 0.

        This method uses the Morgan fingerprints of radius 4 to
        evaluate how similar the molecules in `mols` are.

        Parameters
        ----------
        mols : :class:`iterable` of :class:`rdkit.Chem.rdchem.Mol`
            A group of molecules to which similarity is compared.

        Returns
        -------
        :class:`list`
            A :class:`list` of the form,

            .. code-block:: python

                returned_list = [(8.9, mol1), (7.3, mol2), (3.4, mol3)]

            where the :class:`float` is the similarity of a given
            molecule in `mols` while the ```mol`` is corresponding
            ``rdkit`` molecule. Most similar molecule yielded first.

        """

        # First get the fingerprint of `self`.
        rdkit.GetSSSR(self.mol)
        self.mol.UpdatePropertyCache(strict=False)
        fp = rdkit.GetMorganFingerprint(self.mol, 4)

        # For every structure file in the database create a rdkit
        # molecule. Place these in a list.
        similarities = []
        for mol in mols:
            rdkit.GetSSSR(mol)
            mol.UpdatePropertyCache(strict=False)
            mol_fp = rdkit.GetMorganFingerprint(mol, 4)
            similarity = DataStructs.DiceSimilarity(fp, mol_fp)
            similarities.append((similarity, mol))

        return sorted(similarities, reverse=True, key=lambda x: x[0])

    @classmethod
    def smiles_init(cls,
                    smiles,
                    functional_group=None,
                    note="",
                    name=""):
        """
        Initialize from a SMILES string.

        Parameters
        ----------
        smiles : :class:`str`
            A SMILES string of the molecule.

        functional_group : :class:`str`, optional
            The name of the functional group which is to have atoms
            tagged. If no functional group is provided to this
            parameter, no tagging is done.

        note : :class:`str`, optional
            A note or comment about the molecule.

        name : :class:`str`, optional
            A name which can be optionally given to the molcule for
            easy identification.

        Returns
        -------
        :class:`StructUnit`
            A :class:`StructUnit` instance of the molecule in `smarts`.

        """

        mol = rdkit.MolFromSmiles(smiles)
        rdkit.SanitizeMol(mol)
        mol = rdkit.AddHs(mol)
        key = cls.gen_key(mol, functional_group)
        if key in cls.cache and CACHE_SETTINGS['ON']:
            return cls.cache[key]

        rdkit.EmbedMolecule(mol, rdkit.ETKDG())
        obj = cls.__new__(cls)
        obj.file = smiles
        obj.key = key
        obj.mol = mol
        obj.func_grp = next((x for x in functional_groups if
                            x.name == functional_group), None)
        if obj.func_grp:
            obj.tag_atoms()

        Molecule.__init__(obj, note, name)

        cls.cache[key] = obj
        return obj

    def tag_atoms(self):
        """
        Adds atom properties to atoms.

        All atoms which form the functional group of the molecule have
        the property ``'fg'`` added. Its value is set to the name of
        the functional group. In addition each such atom is given
        the property ``'fg_id'`` which is unique to each functional
        group.

        The atoms which form bonds during assembly have the property
        called ``'bonder'`` added and set to ``1``. Atoms which are
        deleted during reactions have the property ``'del'`` set to
        ``1``.

        Returns
        -------
        None : :class:`NoneType`

        """

        for fg_id, fg in enumerate(self.functional_group_atoms()):
            for atom_id in fg:
                atom = self.mol.GetAtomWithIdx(atom_id)
                atom.SetProp('fg', self.func_grp.name)
                atom.SetIntProp('fg_id', fg_id)

        # Give all atoms which form bonds during reactions the tag
        # 'bonder' and set its value to '1'. Add their ids to
        # `bonder_ids`.
        bond_mol = rdkit.MolFromSmarts(self.func_grp.bonder_smarts)
        bond_atoms = self.mol.GetSubstructMatches(bond_mol)
        for atom_id in flatten(bond_atoms):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetIntProp('bonder', 1)

        # Give all atoms which form bonds during reactions the tag
        # 'del' and set its value to '1'.
        del_mol = rdkit.MolFromSmarts(self.func_grp.del_smarts)
        del_atoms = self.mol.GetSubstructMatches(del_mol)
        for atom_id in flatten(del_atoms):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetIntProp('del', 1)

    def untag_atoms(self):
        """
        Removes the tags added by :meth:`tag_atoms`.

        Returns
        -------
        None : :class:`NoneType`

        """

        for atom in self.mol.GetAtoms():
            atom.ClearProp('fg')
            atom.ClearProp('bonder')
            atom.ClearProp('del')
            atom.ClearProp('fg_id')

    def __str__(self):
        return "{} {}".format(self.__class__.__name__, list(self.key))

    def __repr__(self):
        return str(self)


class StructUnit2(StructUnit):
    """
    Represents building blocks with 2 functional groups.

    """

    def set_orientation2(self, end, conformer=-1):
        """
        Rotate the molecule so that bonder atoms lie on `end`.

        The molecule is rotated about the centroid of the bonder atoms.
        It is rotated so that the direction vector running between the
        2 bonder centroids is aligned with the vector `end`.

        Parameters
        ----------
        end : :class:`numpy.array`
            The vector with which the molecule's bonder atoms should be
            aligned.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        *_, start = next(self.bonder_direction_vectors(conformer))
        return self._set_orientation2(start, end, conformer)

    def minimize_theta2(self, vector, axis, conformer=-1):
        """
        Rotates molecule about `axis` to minimze theta with `vector`.

        The molecule is rotated about `axis` passing through the bonder
        centroid. It is rotated so that the vector between the bonder
        and molecular centroids lies on the same plane as `vector`.

        Parameters
        ----------
        vector : :class:`numpy.array`
            The vector to which the distance should be minimized.

        axis : :class:`numpy.array`
            The direction vector along which the rotation happens.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.minimize_theta(self.centroid_centroid_dir_vector(conformer),
                            vector,
                            axis,
                            self.bonder_centroid(conformer),
                            conformer)


class StructUnit3(StructUnit):
    """
    Represents building blocks with 3 functional groups.

    """

    def bonder_plane(self, conformer=-1):
        """
        Returns coeffs of the plane formed by the bonder centroids.

        A plane is defined by the scalar plane equation::

            ax + by + cz = d.

        This method returns the ``a``, ``b``, ``c`` and ``d``
        coefficients of this equation for the plane formed by the
        bonder centroids. The coefficents ``a``, ``b`` and ``c``
        describe the normal vector to the plane. The coefficent ``d``
        is found by substituting these coefficients along with the
        ``x``, ``y`` and ``z`` variables in the scalar equation and
        solving for ``d``. The variables ``x``, ``y`` and ``z`` are
        substituted by the coordinates of some point on the plane. For
        example, the position of one of the bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.array`
            This array has the form ``[a, b, c, d]`` and represents the
            scalar equation of the plane formed by the bonder
            centroids.

        References
        ----------
        https://tinyurl.com/okpqv6

        """

        centroid = next(self.bonder_centroids(conformer))
        d = -np.sum(self.bonder_plane_normal(conformer) * centroid)
        return np.append(self.bonder_plane_normal(conformer), d)

    def bonder_plane_normal(self, conformer=-1):
        """
        Returns the normal to the plane formed by bonder centroids.

        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the molecule.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.array`
            A unit vector which describes the normal to the plane of
            the bonder centroids.

        """

        if sum(1 for _ in self.bonder_direction_vectors(conformer)) < 2:
            raise ValueError(("StructUnit3 molecule "
                             "has fewer than 3 functional groups."))

        vgen = (v for *_, v in self.bonder_direction_vectors(conformer))
        v1, v2 = it.islice(vgen, 2)

        normal_v = normalize_vector(np.cross(v1, v2))

        theta = vector_theta(normal_v,
                             self.centroid_centroid_dir_vector(conformer))

        if theta > np.pi/2:
            normal_v *= -1

        return normal_v

    def minimize_theta2(self, fg_id, vector, axis, conformer=-1):
        """
        Rotate molecule to minimize angle between `fg_id` and `vector`.

        The rotation is done about `axis` and is centered on the
        bonder centroid, as given by
        :meth:`~.StructUnit.bonder_centroid`.

        Parameters
        ----------
        fg_id : :class:`int`
            The id of functional group which is to have angle
            minimized.

        vector : :class:`numpy.array`
            A vector with which the angle is minimized.

        axis : :class:`numpy.array`
            The vector about which the rotation happens.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        centroid = self.bonder_centroid(conformer)

        fg_centroids = self.bonder_centroids(conformer)
        for i in range(fg_id+1):
            fg_centroid = next(fg_centroids)
        v1 = fg_centroid - centroid
        self.minimize_theta(v1, vector, axis, centroid, conformer)

    def set_orientation2(self, end, conformer=-1):
        """
        Rotates the molecule so the plane normal is aligned with `end`.

        Here "plane normal" referes to the normal of the plane formed
        by the bonder centroids. The molecule is rotated about
        :meth:`~Molecule.bonder_centroid`. The rotation results in the
        normal of the plane being aligned with `end`.

        Parameters
        ----------
        end : :class:`numpy.array`
            The vector with which the normal of plane of bonder
            centroids shoould be aligned.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        start = self.bonder_plane_normal(conformer)
        return self._set_orientation2(start, end, conformer)


@total_ordering
class MacroMolecule(Molecule, metaclass=Cached):
    """
    A representing assembled macromolecules.

    Because of the computational cost associated with macromolecule
    assembly, instances of this class are cached. This means that
    providing the same arguments to the initializer will not build a
    different instance with the same attribute values. It will yield
    the original instance, retrieved from memory.

    This class is not intended to be used directly but should be
    inherited by subclasses representing specific macromolecules. The
    :class:`Cage` and :class:`Polymer` classes are examples of this.
    Any information or methods that apply generally to all
    macromolecule should be defined within this class while specific
    non-general data should be included in the derived classes.

    Note the equality operator ``==`` compares :attr:`fitness`. This
    means two different macromolecules compare equal if they happen to
    have the same fitness. The operator is not to be used to check if
    one macromolecule is the same structurally as another. To do this,
    use :meth:`~Molecule.same`.

    Attributes
    ----------
    building_blocks : :class:`list` of :class:`StructUnit`
        This attribute holds :class:`StructUnit` instances which
        represent the monomers forming the macromolecule. Only one
        :class:`StructUnit` instance is needed per building block, even
        if multiples of that molecule join up to form the
        macromolecule.

    bb_counter : :class:`collections.Counter`
        A counter keeping track of how much of each building block is
        used to form the macromolecule.

    topology : :class:`.Topology`
        Defines the shape of macromolecule and assembles it.

    fitness : :class:`float`
        The fitness value of the macromolecule, used by the GA.

    unscaled_fitness : :class:`dict`
        The dictionary holds the name of a fitness function as the
        key and the value it calculated for unscaled_fitness as the
        value. For example,

        .. code-block:: python

            unscaled_fitness = {'fitness_func1': 12.3,
                                'fitness_func2': 49.2}

        where ``'fitness_func1'`` and ``'fitness_func2'`` are the
        names of the fintess functions applied on the molecule.

    progress_params : :class:`dict`
        Holds the fitness parameters which the GA should track to make
        progress plots. The key is the name of a fitness function. Has
        the form

        .. code-block:: python

            unscaled_fitness = {'fitness_func1': [8, 49]
                                'fitness_func2': [78, 4.2, 32.3]}

    bonds_made : :class:`int`
        The number of bonds made during assembly.

    """

    def __init__(self,
                 building_blocks,
                 topology,
                 name="",
                 note="",
                 bb_conformers=None):
        """
        Initialize a :class:`MacroMolecule` instance.

        Parameters
        ---------
        building_blocks : :class:`list` of :class:`StructUnit`
            The :class:`StructUnit` instances of building blocks
            forming the macromolecule.

        topology : :class:`.Topology`
            Defines the shape of macromolecule and assembles it.

        name : :class:`str`, optional
            A name which can be given to the molcule for easy
            identification.

        note : :class:`str`, optional
            A note or comment about the molecule.

        bb_conformers : :class:`list` of :class:`int`, optional
            The ids of the building block conformers to be used. Must
            be equal in length to `building_blocks` and orders must
            correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        """

        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(building_blocks))]

        self.fitness = None
        self.unscaled_fitness = {}
        self.progress_params = {}
        self.building_blocks = building_blocks
        self.bb_counter = Counter()
        self.topology = topology

        try:
            # Ask the ``Topology`` instance to assemble/build the
            # macromolecule. This creates the `mol` attribute.
            topology.build(self, bb_conformers)

        except Exception as ex:
            self.mol = rdkit.Mol()
            errormsg = ('Build failure.\n'
                        '\n'
                        'topology\n'
                        '--------\n'
                        '{}\n'
                        '\n'
                        'building blocks\n'
                        '---------------\n').format(topology)

            bb_blocks = []
            for i, bb in enumerate(building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append(
                    ('{0.__class__.__name__} {0.func_grp.name}\n'
                     '{1}').format(bb, bb.mdl_mol_block(bb_conf)))

            errormsg += '\n'.join(bb_blocks)

            logger.error(errormsg, exc_info=True)

        super().__init__(name, note)

    def add_conformer(self, bb_conformers):
        """
        Assembles a new conformer.

        Parameters
        ----------
        bb_conformers : :class:`list` of :class:`int`
            The ids of the building block conformers to be used. Must
            be equal in length to :attr:`building_blocks` and orders
            must correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        Returns
        -------
        :class:`int`
            The id of the new conformer.

        """

        # Save the original rdkit molecule.
        original_mol = self.mol
        # Build a new molecule.
        try:
            # Ask the ``Topology`` instance to assemble/build the
            # macromolecule. This creates the `mol` attribute.
            self.topology.build(self, bb_conformers)

        except Exception as ex:
            self.mol = rdkit.Mol()
            errormsg = ('Build failure.\n'
                        '\n'
                        'topology\n'
                        '--------\n'
                        '{}\n'
                        '\n'
                        'building blocks\n'
                        '---------------\n').format(self.topology)

            bb_blocks = []
            for i, bb in enumerate(self.building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append(
                    ('{0.__class__.__name__} {0.func_grp.name}\n'
                     '{1}').format(bb, bb.mdl_mol_block(bb_conf)))

            errormsg += '\n'.join(bb_blocks)

            logger.error(errormsg, exc_info=True)

        # Get the new conformer.
        new_conf = rdkit.Conformer(self.mol.GetConformer())
        # Add it to the original molecule.
        new_id = original_mol.AddConformer(new_conf, True)
        self.mol = original_mol
        return new_id

    def building_block_cores(self, bb):
        """
        Yields the "cores" of the building block molecules.

        The structure of the yielded cores has the geometry found in
        the macromolecule.

        Parameters
        ----------
        bb : :class:`int`
            The index of a building block molecule within
            :attr:`building_blocks`. The cores of this molecule are
            yielded.

        Yields
        ------
        :class:`rdkit.Chem.rdchem.Mol`
            The core of a building block molecule, as found in the
            macromolecule.

        """

        done = set()
        for atom_id, props in self.atom_props.items():
            # Ignore fragments which do not correspond to the molecule
            # `bb`.
            if props['bb_index'] != bb or props['mol_index'] in done:
                continue

            done.add(props['mol_index'])
            # For each fragment make a new core. To preserve bond
            # information, start with the macromolecule.
            coremol = rdkit.EditableMol(self.mol)
            # Remove any atoms not part of the core - atoms not in the
            # fragment itself, hydrogens, atoms in the functional
            # group.
            for atom in reversed(self.mol.GetAtoms()):
                atomid = atom.GetIdx()
                bb_index = self.atom_props[atomid]['bb_index']
                mol_index = self.atom_props[atomid]['mol_index']
                fg = self.atom_props[atomid].get('fg', None)
                if (bb_index != props['bb_index'] or
                    mol_index != props['mol_index'] or
                    atom.GetAtomicNum() == 1 or
                   fg is not None):
                    coremol.RemoveAtom(atomid)
            yield coremol.GetMol()

    def json(self):
        """
        Returns a JSON representation of the molecule.

        The representation has the form

        .. code-block:: python

            {
                'class' : 'Polymer',
                'mol_block' : '''A string holding the V3000 mol
                                 block of the molecule.'''
                'building_blocks' : {bb1.json(), bb2.json()}
                'topology' : 'Copolymer(repeating_unit="AB")'
                'unscaled_fitness' : {'fitness_func1' : fitness1,
                                      'fitness_func2' : fitness2},
                'note' : 'A nice molecule.',
                'name' : 'Poly-Benzene'
            }

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        return {
            'bb_counter': [(key.json(), val) for key, val in
                           self.bb_counter.items()],
            'bonds_made': self.bonds_made,
            'class': self.__class__.__name__,
            'mol_block': self.mdl_mol_block(),
            'building_blocks': [x.json() for x in
                                self.building_blocks],
            'topology': repr(self.topology),
            'unscaled_fitness': repr(self.unscaled_fitness),
            'progress_params': self.progress_params,
            'note': self.note,
            'name': self.name,
            'atom_props': self.atom_props

        }

    @classmethod
    def _json_init(cls, json_dict):
        """
        Completes a JSON initialization.

        This function is not to be used. Use :meth:`Molecule.load`
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : :class:`dict`
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        bbs = [Molecule.from_dict(x) for x in
               json_dict['building_blocks']]

        topology = eval(json_dict['topology'],  topologies.__dict__)

        key = cls.gen_key(bbs, topology)
        if key in cls.cache and CACHE_SETTINGS['ON']:
            return cls.cache[key]

        obj = cls.__new__(cls)
        obj.mol = rdkit.MolFromMolBlock(json_dict['mol_block'],
                                        sanitize=False,
                                        removeHs=False)
        obj.topology = topology
        obj.unscaled_fitness = eval(json_dict['unscaled_fitness'],
                                    np.__dict__)
        obj.fitness = None
        obj.progress_params = json_dict['progress_params']
        obj.bb_counter = Counter({Molecule.from_dict(key): val for
                                  key, val in json_dict['bb_counter']})
        obj.bonds_made = json_dict['bonds_made']
        obj.energy = Energy(obj)
        obj.optimized = json_dict['optimized']
        obj.note = json_dict['note']
        obj.name = json_dict['name'] if json_dict['load_names'] else ""
        obj.key = key
        obj.building_blocks = bbs
        obj.atom_props = {int(key): value for key, value in
                          json_dict['atom_props'].items()}

        # Remake bonder_ids
        obj.bonder_ids = bonder_ids = []
        for atom_id, props in obj.atom_props.items():
            if 'bonder' in props:
                fg_id = props['fg_id']
                if len(bonder_ids) < fg_id + 1:
                    diff = fg_id + 1 - len(bonder_ids)
                    bonder_ids.extend([] for i in range(diff))
                bonder_ids[fg_id].append(atom_id)

        if CACHE_SETTINGS['ON']:
            cls.cache[key] = obj

        return obj

    @staticmethod
    def gen_key(building_blocks, topology):
        """
        Generates the key used for caching the molecule.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`StructUnit`
            The building blocks used to make the macromolecule.

        topology : :class:`.Topology`
            The topology used to make the macromolecule.

        Returns
        -------
        :class:`tuple`
            The key used for caching the macromolecule.

        """

        return (frozenset(x.key for x in building_blocks),
                repr(topology))

    def bb_distortion(self, bb_conformers=None, conformer=-1):
        """
        Rmsd difference of building blocks before and after assembly.

        The function looks at each building block in the macromolecule
        and calculates the rmsd between the "free" verson and the one
        present in the macromolecule. The mean of these rmsds is
        returned.

        Atoms which form the functional group of the building blocks
        and hydrogens are excluded from the calculation.

        Parameters
        ----------
        bb_conformers : :class:`list` of :class:`int`
            The ids of building block conformers to use. 1 id for each
            building block, in an order corresponding to
            :attr:`building_blocks`. If ``None``, all conformer ids
            default to ``-1``.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The mean rmsd of the macromole's building blocks to their
            "free" counterparts.

        """

        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(self.building_blocks))]

        # Go through each of the building blocks. For each building
        # block get the core. Get the corrospending cores in the
        # macromolecules and add the rmsd to the sum. Increment the
        # count to calculate the mean later.
        rmsd = 0
        n = 0
        for i, bb in enumerate(self.building_blocks):
            free = bb.core()
            am = [(x, x) for x in range(free.GetNumAtoms())]
            for frag in self.building_block_cores(i):
                rmsd += rdkit.AlignMol(free,
                                       frag,
                                       bb_conformers[i],
                                       conformer,
                                       atomMap=am)
                n += 1
        return rmsd / n

    def update_cache(self):
        """
        Update attributes of cached molecule.

        When an instance of :class:`MacroMolecule` is first created it
        is cached. Using ``multiprocessing`` to perform optimizations
        or calculate fitness returns modified copies of the cached
        molecules. In order to ensure that the cached molecules have
        their attributes updated to the values of the copies, this
        method must be run on the copies.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.key in self.__class__.cache:
            self.__class__.cache[self.key].__dict__ = dict(vars(self))

    def __eq__(self, other):
        return self.fitness == other.fitness

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __str__(self):
        return "{}(building_blocks={}, topology={!r})".format(
                        self.__class__.__name__,
                        [str(x) for x in self.building_blocks],
                        self.topology)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return id(self)

    """
    The following methods are inteded for convenience while
    debugging or testing and should not be used during typical
    execution of the program.

    """

    @classmethod
    def testing_init(cls, bb1, bb2, topology):

        key = frozenset({bb1, bb2}), repr(topology)

        if key in MacroMolecule.cache.keys():
            return MacroMolecule.cache[key]
        else:
            macro_mol = cls.__new__(cls)
            macro_mol.building_blocks = [bb1, bb2]
            macro_mol.topology = topology
            MacroMolecule.cache[key] = macro_mol
            return macro_mol


class Cage(MacroMolecule):
    """
    Used to represent molecular cages.

    """

    def window_difference(self, conformer=-1):
        """
        The total difference in all window sizes.

        Every combination of windows is considered and all the
        size differences are summed and returned. Only
        differences between windows of the same type are
        considered.

        Consider a triangular-based prism cage topology. Such a
        cage will have triangular windows and square windows. You
        only want to compare the triangulars with other
        triangular windows and squares only with other squares.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The total difference of window size when considering
            every combination of windows of the same type.

        None : :class:`NoneType`
            If not all windows were found.

        """

        windows = self.windows(conformer)

        if windows is None or len(windows) < self.topology.n_windows:
            return None

        windows = np.array(windows)

        # Cluster the windows into groups so that only size
        # differences between windows of the same type are taken
        # into account. To do this, first sort the windows by
        # size. If two windows types are present split the
        # windows at the two groups at the point where the window
        # sizes have the biggest difference. If there are three
        # types split it at the two biggest differences and so
        # on.

        diffs = list(abs(np.ediff1d(windows)))
        sorted_diffs = sorted(diffs, reverse=True)

        # Get indices of where the list should be split.
        split = []
        for x in range(self.topology.n_window_types-1):
            i = diffs.index(sorted_diffs[x]) + 1
            split.append(i)

        # Get the sub-lists.
        og = list(windows)
        clusters = []
        for i in sorted(split, reverse=True):
            clusters.append(og[i:])
            og = og[:i]

        if self.topology.n_window_types == 1:
            clusters.append(og)

        # After this sum the differences in each group and then
        # sum the group totals.
        diff_sums = []
        for cluster in clusters:
            diff_sum = sum(abs(w1 - w2) for w1, w2 in
                           it.combinations(cluster, 2))

            diff_num = sum(1 for _ in it.combinations(cluster, 2))

            diff_sums.append(diff_sum / diff_num)

        return sum(diff_sums)

    def windows(self, conformer=-1):
        """
        Returns window sizes found by ``pyWindow``.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        None : :class:`NoneType`
            If the function for finding windows and their sizes
            found fewer than the required number of windows or
            if it failed for some other reason.

        :class:`list` of :class:`float`
            Each :class:`float` represents the size of a
            window in the cage. If the window finding function
            found more than the expected number of windows, only
            the largest ``n`` windows are returned. Where ``n`` is the
            number of expected windows for that cage topology.

        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Load an RDKit molecule object to pyWINDOW.

            # As pyWindow doesnt support multiple conformers, first
            # make an rdkit molecule holding only the desired conformer.
            new_mol = rdkit.Mol(self.mol)
            new_mol.RemoveAllConformers()
            new_mol.AddConformer(self.mol.GetConformer(conformer))

            pw_molecule = pywindow.molecular.Molecule.load_rdkit_mol(new_mol)
            # Find windows and get a single array with windows' sizes.
            all_windows = pw_molecule.calculate_windows(output='windows')

        # If pyWindow failed, return ``None``.
        if all_windows is None:
            return None

        return sorted(all_windows[0], reverse=True)[:self.topology.n_windows]


class Polymer(MacroMolecule):
    """
    Used to represent polymers.

    """

    pass


class Periodic(MacroMolecule):
    """
    Used to represent periodic structures.

    This class is essentially the same as :class:`MacroMolecule`,
    with additional methods and attributes relevant to periodic
    materials being added.

    Attributes
    ----------
    deleters : :class:`dict`
        A :class:`dict` of the form

        .. code-block:: python

            {
                12: [[coord, elem, bond_type],
                     [coord, elem, bond_type],
                     [coord, elem, bond_type]],

                45: [[coord, elem, bond_type],
                     [coord, elem, bond_type],
                     [coord, elem, bond_type]]
            }

        The key is an :class:`int` which represents a bonder atom.
        The value holds the position, element and bond type of
        every deleter atom removed from that bonder. The position is a
        :class:`numpy.array`, the element is an :class:`int` and the
        bond type is an :class:`rdkit.Chem.rdchem.BondType`.

    periodic_bonds : :class:`list` of :class:`.PeriodicBond`
        When periodic topologies are assembled, periodic bonds
        do not get added to the ``rdkit`` molecule in the
        :attr:`~.MacroMolecule.mol` attribute. Instead,
        :meth:`~.PeriodicLattice.join_mols` adds
        :class:`.PeriodicBond` instances representing the bonds into
        this list.

    cell_dimensions : :class:`list` of :class:`numpy.array`
        The dimensions of the unit cell. The first array is the vector
        ``a`` the second is ``b`` and the third is ``c``. This should
        be added during the build process by the periodic topology.

    """

    def __init__(self, building_blocks, topology, name="", note=""):
        self.periodic_bonds = []
        self._ids_updated = False
        super().__init__(building_blocks, topology, name, note)

    def island(self, dimensions):
        """
        Build a terminated supercell.

        Terminated means that the periodic bonds are replaced with
        bonds to terminating atoms.

        Parameters
        ----------
        dimensions : :class:`list` of :class:`int`
            A 3 member :class:`list`, holding the number of unit cells
            in the x, y and z directions used to make the supercell.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            An ``rdkit`` molecule of the island.

        """

        cells, island = self._place_island(dimensions)
        return self._join_island(cells, island)

    def _join_island(self, cells, island):
        """
        Adds bonds between unit cells of `island`.

        Notes
        -----
        For internal use by :meth:`island`.

        Parameters
        ----------
        cells : nested :class:`list` of :class:`Cell`

        island : :class:`rdkit.Chem.rdchem.Mol`
            The island molecule holding unit cells placed side by
            side like in a supercell but with no bonds running between
            them.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The island with bonds added.

        """

        # `self.periodic_bonds` holds objects of the
        # ``PeriodicBond`` class. Each ``PeriodicBond`` object has the
        # ids of two fgs in the unit cell which are connected
        # by a bond running across the periodic boundary. The
        # `direction` attribute descibes the axes along which the
        # bond is periodic. For example, if `direction1` is [1, 0, 0]
        # it means that the fg in `periodic_bond.fg1` has a
        # perdiodic bond connecting it to `periodic_bond.fg2` going
        # in the positive direction along the x-axis.

        # When iterating through all the unit cells composing the
        # island, you can use the `direction` vector to get index of
        # the unit cell which holds fg connected the present cell.
        # Then just form bonds between the correct fgs by mapping
        # the fg ids in the unit cells to the ids of the equivalent
        # fgs in the original unit cell  and checking the
        # `periodic_bond` to see which fg ids are connected.
        for cell in flatten(cells):
            for periodic_bond in self.periodic_bonds:

                # Get the indices of the cell which holds the atom
                # bonded to the equivalent atom of
                # `periodic_bond.atom1` in the present `cell`.
                x, y, z = cell.id + periodic_bond.direction
                if (x < 0 or y < 0 or z < 0 or
                    x >= len(cells) or
                    y >= len(cells[0]) or
                   z >= len(cells[0][0])):
                    continue

                # ccel as in "connected cell".
                ccell = cells[x][y][z]

                # `fg1` is the id of a fg, found in `cell`
                # and equivalent to `periodic_bond.fg1`, having a
                # bond added.
                fg1 = cell.fgs[periodic_bond.fg1]
                # `fg2` is the id of a fg, found in
                # `ccell` and equivalent to `periodic_bond.fg2`,
                # having a bond added.
                fg2 = ccell.fgs[periodic_bond.fg2]

                island, _ = react(island, True, fg1, fg2)

        return island

    def _place_island(self, dimensions):
        """
        Places unit cells side by side to form an island.

        Notes
        -----
        For internal use by :meth:`island`.

        Parameters
        ----------
        dimensions : :class:`list` of :class:`int`
            The number of unit cells in the x, y and z directions to be
            placed side by side.

        Returns
        -------
        :class:`tuple`
            The first member of the tuple is a :class:`list` holding
            :class:`Cell` objects, one for each unit cell placed. The
            :class:`Cell` object is placed within nested lists so that
            a cell with the coordinates x, y, z can be accessed from
            the list using ``[x][y][z]``. For example,

                >>> cells, island, bonder_map = \
periodic._place_island([4, 4, 4])
                >>> cells[2][1][3]
                <Cell at 0x7fa0155d54e0>

            where the returned :class:`Cell` object represents the
            3rd unit cell along the x axis, the second along the y axis
            and the fourth along the z axis.

            The second member is an ``rdkit`` molecule of the island
            being built. The third member is a :class:`dict` mapping
            the ids of fgs in the island back to the id of the
            equivalent fg in the original unit cell.

        """

        a, b, c = self.cell_dimensions
        cells = np.full(dimensions, None, object).tolist()
        island = rdkit.Mol()

        xdim, ydim, zdim = (range(d) for d in dimensions)
        nfgs = 1 + max(atom.GetIntProp('fg_id') for atom in
                       self.mol.GetAtoms() if atom.HasProp('fg_id'))
        # For each dimension place a unit cell.
        for i, (x, y, z) in enumerate(it.product(xdim, ydim, zdim)):
            unit_cell = self.shift(x*a + y*b + z*c)

            # Update fg_ids.
            for atom in unit_cell.GetAtoms():
                if not atom.HasProp('fg_id'):
                    continue
                atom.SetIntProp('fg_id', atom.GetIntProp('fg_id')+i*nfgs)

            island = rdkit.CombineMols(island, unit_cell)
            # `bonders` maps a bonder id in the original unit
            # cell to the one currently being added to the
            # island.
            fgs = {fg: i*len(self.bonder_ids) + fg for
                   fg in range(len(self.bonder_ids))}
            cells[x][y][z] = Cell((x, y, z), fgs)

        return cells, island

    def write_gulp_input(self, path, keywords,
                         cell_fix=[0, 0, 0, 0, 0, 0], atom_fix=None):
        """
        Writes a GULP input file of the unit cell.

        Parameters
        ----------
        path : :class:`str`
            The `path` of the file to which the molecule should be
            written.

        keywords : :class:`list` of :class:`str`
            The keywords to be placed on the first line of the input
            file.

        cell_fix : :class:`list` of :class:`int`, optional
            A 6 member list holding the fix parameters for the unit
            cell.

        atom_fix : :class:`numpy.array` of :class:`int`, optional
            An n by 3 array where n is the number of atoms in the
            unit cell. Each row has the fix parameters for a given
            atom.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atom_fix is None:
            atom_fix = np.ones([self.mol.GetNumAtoms(), 3])

        with open(path, 'w') as f:
            f.write(' '.join(keywords) + '\n\n')
            f.write('name {}\n\n'.format(self.name))
            # Write the cell parameters.
            f.write('cell\n')
            # The sizes of cell vectors a, b and c are written first.
            for vector in self.cell_dimensions:
                f.write(str(np.round(np.linalg.norm(vector), 6)) + ' ')
            # Then angles alpha, beta and gamma.
            a, b, c = self.cell_dimensions
            angle1 = round(math.degrees(vector_theta(a, c)), 6)
            angle2 = round(math.degrees(vector_theta(b, c)), 6)
            angle3 = round(math.degrees(vector_theta(a, b)), 6)
            f.write(str(angle1) + ' ')
            f.write(str(angle2) + ' ')
            f.write(str(angle3))
            # Finally the fix parameters for the cell.
            for fix in cell_fix:
                f.write(' ' + str(fix))
            f.write('\n')
            # Add atom coordinates.
            f.write('cart\n')
            for (id_, coords), fix in zip(self.all_atom_coords(),
                                          atom_fix):
                # Don't write deleter atoms.
                if self.atom_props.get(id_, {}).get('del', False):
                    continue

                x, y, z = [round(x, 4) for x in coords]
                fx, fy, fz = [int(x) for x in fix]
                f.write('{} core {} {} {} {} {} {}\n'.format(
                         self.atom_symbol(id_), x, y, z, fx, fy, fz))
            f.write('\n')
            # Add bonds.
            for bond in self.mol.GetBonds():
                a1 = bond.GetBeginAtomIdx() + 1
                a2 = bond.GetEndAtomIdx() + 1
                f.write('connect {} {} 0 0 0\n'.format(a1, a2))

            # Add periodic bonds.
            for bond in self.periodic_bonds:
                a1 = next(a for a, props in self.atom_props.items() if
                          props.get('fg_id', None) == bond.fg1 and
                          props.get('bonder', False)) + 1
                a2 = next(a for a, props in self.atom_props.items() if
                          props.get('fg_id', None) == bond.fg2 and
                          props.get('bonder', False)) + 1
                dx, dy, dz = bond.direction
                f.write(f'connect {a1} {a2} {dx:+} {dy:+} {dz:+}\n')
