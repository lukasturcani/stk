"""
Defines classes which describe molecules.

There are a couple of major classes defined here. The most important
ones are ``StructUnit`` and ``MacroMolecule``. Here is an overview of
their role/interaction. This is followed by a step-by-step guide to
macromolecular assembly.

The StructUnit class represents the monomers that make up
macromolecules. These are commonly refered to as ``building blocks`` in
the documentation. The class holds information concerning only a single
building block molecule. Such as the number of atoms and bonds it may
have. It also has information about the functoinal groups present on
the building block molecule (FGInfo class defined in
/mmea/classes/fg_info.py). The class also allows manipulation of the
lone building block molecule, such as rotations and translations.

The StructUnit class should be inherited as necessary. For example
StructUnit2 adds manipulations relavant to molecules with 2 functional
groups. The StructUnit3 class adds manipulations relavant to 3 or more
functional groups. If you have a monomer which needs specific
information or manipulations, give it its own class.

The MacroMolecule class represents an assembled macromolecule. It
requires at least 2 basic pieces of information: which monomers are
used to assemble the macromolecule and what is the topology/structure
of the macromolecule.

The MacroMolecule holds in its `building_blocks` attribute a list of
StructUnit (or derived classes) instances. These represent the
monomers which form the macromolecule. Only one StructUnit instance per
monomer type is held. So if 4 of one type of monomer and 2 of another
type of monomer form a macromolecule, only 2 StructUnit instances are
in `building_blocks`.

In its `topology` attribute a MacroMolecule holds a ``Topology` child
class instance. This instance is responsible for assembling the
macromolecule from the building blocks. The building should happen in
the ``__init__()`` method of the MacroMolecule via the Topology
instance's ``build()`` method. The ``build()`` method should place the
assembled macromoleclue in the `mol` attribute of the MacroMolecule
instance.

The StructUnit class labels atoms in the functional groups of the
building blocks as either ``bonder`` or ``del`` (see its
documentation). This tells the Topology intance which atoms form bonds
and which are removed during assembly.

No need to worry about the metaclasses.

-------------------------------------------------------
A more detailed description of macromolecular assembly.
-------------------------------------------------------
This is a step-by-step guide of how macromolecular assembly is carried
out and what the classes do.

First you create StructUnit instances of the building blocks which make
up the macromolecule:

    bb = StructUnit('/path/to/struct/file.mol2')

The StructUnit instances are initialized using paths of molecular
structure files. (Initializing a StructUnit automatically completes
steps 1 to 4.)

    1) Place an rdkit instance of the molecule into the `mol` attribute
       of the StructUnit.
    2) Scan the path of the structure file for the name of a functional
       group. (Alternatively the name of a functional group can be
       supplied to the initializer).

What functional groups are recognized by MMEA?
The module ``/mmea/molecular/fg_info.py`` defines a class called
FGInfo. Instances of this class are held in the list
`functional_groups` (also in that module). If you put an FGInfo
instance in that list, the functional group will be recognized.

    3) Place the FGInfo instance of the functional group in the
       `func_grp` attribute of the StructUnit.

    4) Using the FGInfo instance, tag atoms in the building block as
       either 'bonder' or 'del'. 'bonder' signifies that the atoms form
       a bond during macromolecular assembly, while 'del' means they
       are deleted.

    5) Give the StructUnit instances and the class representing the
       topology to the macromolecule's initializer.

    6) Make an instance of the topology class in the MacroMolecule
      initializer.

    7) Run the ``build()`` of the instance generated in 6) during
       MacroMolecule initialization.

    8) The ``build()`` will vary depending on the Topology class.
       However the basic structure is the same (steps 9 - 11).

    9) Combine the StructUnit rdkit molecules into a single rdkit
       molecule. Make sure that the building blocks are arranged in the
       shape of the macromolecule. All the manipulations available via
       the StructUnit class are useful here to make sure all the
       building blocks are oriented correctly when forming the
       macromolecule.

    10) Create bonds between all the disjoined building block
        molecules. This is where the tagging done by StructUnit is
        needed.

    11) Delete all the atoms tagged for deletion.

After all this you should have a rdkit instance of the macromolecule
which should be placed in the `mol` attribute of the MacroMolecule
intance.

Extending MMEA: Adding new macromolecules.
------------------------------------------
To add new macromolecules create a new class which inherits
``MacroMolecule``. The initializer of ``MacroMolecule`` can be used or
a custom one may be used. However, all the same attributes must be
defined.

If you're adding a new class of macromolecules, it quite likely you
want to add a new ``Topology`` class. See the docstring of the
/mmea/molecular_tools/topologies/base.py module for guidance on adding
these. The topology class does the assembly of the macromolecule from
the building blocks. The assembled macromolecule instance is then held
by the macromolecular class.

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

from collections import Counter, defaultdict, ChainMap
from inspect import signature

from . import topologies
from .fg_info import functional_groups
from .energy import Energy
from ..addons.pyWINDOW import pywindow as pywindow
from ..convenience_tools import (flatten, periodic_table,
                                 normalize_vector, rotation_matrix,
                                 vector_theta, mol_from_mae_file,
                                 rotation_matrix_arbitrary_axis,
                                 atom_vdw_radii, bond_dict, Cell)


logger = logging.getLogger(__name__)


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

        if key in self.cache:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            self.cache[key] = obj
            return obj


class CachedStructUnit(type):
    """
    A metaclass for making StructUnit create cached instances.

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
        if key in self.cache:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            self.cache[key] = obj
            return obj


class Molecule:
    """
    The most basic class representing molecules.

    This class defines the operations which any class
    describing molecules should inherit or may find useful. Examples of
    such classes are ``StructUnit`` and ``MacroMolecule``. This calls
    is unlikely to be useful as in and of itself. It lacks an
    `__init__()` method because it does not need one. Child classes
    should define it.

    It is assumed that child classes will have some basic attributes.

    Attributes
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A rdkit molecule instance representing the molecule.

    inchi : str
        The InChI of the molecule.

    energy : Energy
        An instance of the ``Energy`` class. It handles all things
        energy.

    optimized : bool
        Indicates whether a Molecule has been passed through an
        optimization function or not.

    name : str (default = "")
        A name which can be optionally given to the molecule for easy
        identification.

    note : str (default = "")
        A note or comment about the molecule. Purely optional but can
        be useful for labelling and debugging.

    """

    def __init__(self, name="", note=""):
        self.optimized = False
        self.energy = Energy(self)
        self.bonder_ids = []
        self.name = name
        self.note = note

    def all_atom_coords(self):
        """
        Yields the coordinates of atoms in `mol`.

        Yields
        ------
        tuple of (int, numpy.array)
            The ``int`` represents the atom id of the atom whose
            coordinates are being yielded.

            The array represents the position.

        """

        # Get the conformer from the rdkit instance.
        conformer = self.mol.GetConformer()

        # Go through all the atoms and ask the conformer to return
        # the position of each atom. This is done by supplying the
        # conformers `GetAtomPosition` method with the atom's id.
        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_position = conformer.GetAtomPosition(atom_id)
            yield atom_id, np.array([*atom_position])

    def atom_coords(self, atom_id):
        """
        Return coordinates of an atom.

        Parameters
        ----------
        atom_id : int
            The id of the atom whose coordinates are desired.

        Returns
        -------
        numpy.array
            An array holding the x, y and z coordinates of the atom.

        """

        conf = self.mol.GetConformer()
        atom_position = conf.GetAtomPosition(atom_id)
        return np.array([*atom_position])

    def atom_distance(self, atom1_id, atom2_id):
        """
        Return the distance between 2 atoms.

        Parameters
        ----------
        atom1_id : int
            The id of the first atom.

        atom2_id : int
            The id of the second atom.

        Returns
        -------
        scipy.double
            The distance between the first and second atoms.

        """

        # Get the atomic positions of each atom and use the scipy
        # function to calculate their distance in Euclidean space.
        atom1_coords = self.atom_coords(atom1_id)
        atom2_coords = self.atom_coords(atom2_id)
        return euclidean(atom1_coords, atom2_coords)

    def atom_symbol(self, atom_id):
        """
        Returns the symbol of the atom with id `atom_id`.

        Parameters
        ----------
        atom_id : int
            The id number of the atom.

        Returns
        -------
        str
            The atomic symbol of the atom.

        """

        return self.mol.GetAtomWithIdx(atom_id).GetSymbol()

    def _cavity_size(self, origin):
        """
        Calculates diameter of the molecule's from `origin`.

        The cavity is measured by finding the atom nearest to
        `origin`, correcting for van der Waals diameter and multiplying
        by -2.

        This function should not be used. Use cavity_size() instead.
        cavity_size() finds the optimal value of `origin` to use.

        Parameters
        ----------
        origin : numpy.array
            Holds the x, y and z coordinate of the position from which
            the cavity is measured.

        Returns
        -------
        float
            The (negative) diameter of the molecules cavity.

        """

        atom_vdw = np.array([atom_vdw_radii[x.GetSymbol()] for x
                            in self.mol.GetAtoms()])

        distances = euclidean_distances(self.position_matrix().T,
                                        np.matrix(origin))
        distances = distances.flatten() - atom_vdw
        return -2*min(distances)

    def cavity_size(self):
        """
        Calculates the diameter of the molecule's cavity.

        Returns
        -------
        The diameter of the molecule's cavity in Angstroms.

        """

        # This function uses _cavity_size() to calculate the cavity
        # size. _cavity_size() finds the closest atom to `origin` to
        # get its value of the cavity.

        # What this function does is finds the value of `origin` which
        # causes _cavity_size() to calculate the largest possible
        # cavity.
        ref = self.center_of_mass()
        icavity = 0.5*self._cavity_size(ref)
        bounds = [(coord+icavity, coord-icavity) for coord in ref]
        cavity_origin = minimize(self._cavity_size,
                                 x0=ref,
                                 bounds=bounds).x
        cavity = -self._cavity_size(cavity_origin)
        return 0 if cavity < 0 else cavity

    def center_of_mass(self):
        """
        Returns the centre of mass of the molecule.

        Returns
        -------
        numpy.array
            An array holding the coordinates of the center of mass.

        References
        ----------
        https://en.wikipedia.org/wiki/Center_of_mass

        """

        center = np.array([0., 0., 0.])
        total_mass = 0.
        for atom_id, coord in self.all_atom_coords():
            mass = self.mol.GetAtomWithIdx(atom_id).GetMass()
            total_mass += mass
            center += mass*coord
        return np.divide(center, total_mass)

    def centroid(self):
        """
        Returns the centroid of the molecule.

        Returns
        -------
        numpy.array
            A numpy array holding the position of the centroid.

        """

        centroid = sum(x for _, x in self.all_atom_coords())
        return np.divide(centroid, self.mol.GetNumAtoms())

    # @classmethod
    def dihedral_strain(self, dihedral_SMARTS, target):
        """
        Calculates the relative % difference between all the average dihedral
        angle values within the molecule and the target value.

        Parameters
        ----------
        dihedral_SMARTS : str
            The SMARTS code for the dihedral of interest.

        target : float
            Float representing the target value for the dihedral angle.

        Returns
        -------
        diff : float
            Float representing the % relative difference between the average
            dihedral value in the molecule and the target value.
        """

        # Sanitize the molecule
        rdkit.SanitizeMol(self.mol)

        match = rdkit.MolFromSmarts(dihedral_SMARTS)

        atoms_dihedral = self.mol.GetSubstructMatches(match)

        dihedral_info = []
        if len(atoms_dihedral) > 0:
            for atoms_group in atoms_dihedral:
                # Calculate the dihedral angle
                dihedral_value = rdMolTransforms.GetDihedralDeg(
                                    self.mol.GetConformer(), atoms_group[0],
                                    atoms_group[1], atoms_group[2],
                                    atoms_group[3])
                dihedral_info.append(dihedral_value)

            # Calculate the average dihedral value
            avg_dihedral = np.mean([abs(x) for x in dihedral_info])
            # Calculate the relative diff with the target dihedral value
            diff = (abs(target - avg_dihedral) / target) * 100
        else:
            diff = 1

        return diff

    def dump(self, path):
        """
        Writes a JSON dict of the molecule to a file.

        Parameters
        ----------
        path : str
            The full path of the file to which the JSON dict of the
            molecule can be written.

        Returns
        -------
        None : NoneType

        """

        with open(path, 'w') as f:
            json.dump(self.json(), f, indent=4)

    @classmethod
    def fromdict(self, json_dict, optimized=True, load_names=True):
        """
        Creates a Molecule from a JSON dict.

        The Molecule returned has the class specified in the JSON
        dictionary, not ``Molecule``.

        Parameters
        ----------
        json_dict : dict
            A dict holding the JSON representation of a molecule.

        optimized : bool
            The value passed to the loaded molecules `optimized`
            attribute.

        load_names : bool (default = True)
            If ``True`` then the `name` attribute stored in the JSON
            objects is loaded. If ``False`` then it's not.

        Returns
        -------
        Molecule
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
        networkx.Graph
            A graph where the nodes are the ids of the atoms in the
            rdkit molecule and the edges are the bonds.

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
        Returns the InChi of the molecule.

        """

        self.update_stereochemistry()
        return rdkit.MolToInchi(self.mol)

    @classmethod
    def load(cls, path, optimized=True, load_names=True):
        """
        Creates a Molecule from a JSON file.

        The Molecule returned has the class specified in the JSON
        dictionary, not ``Molecule``.

        Parameters
        ----------
        path : str
            The full path holding a JSON representation to a molecule.

        optimized : bool (default = True)
            The value passed to the loaded molecules `optimized`
            attribute.

        load_names : bool (default = True)
            If ``True`` then the `name` attribute stored in the JSON
            objects is loaded. If ``False`` then it's not.

        Returns
        -------
        Molecule
            The molecule held in `path`.

        """

        with open(path, 'r') as f:
            json_dict = json.load(f)

        return cls.fromdict(json_dict, optimized, load_names)

    def max_diameter(self):
        """
        Returns the largest distance between 2 atoms in the molecule.

        Returns
        -------
        tuple of form (float, int, int)
            The float represents the largest inter-atomic distance on
            the molecule. The ints are the ids of the atoms.

        """

        maxid1, maxid2 = max((x for x in
                              it.combinations(
                                range(self.mol.GetNumAtoms()), 2)),
                             key=lambda x: self.atom_distance(*x))

        maxd = self.atom_distance(maxid1, maxid2)
        maxd += (atom_vdw_radii[self.atom_symbol(maxid1)] +
                 atom_vdw_radii[self.atom_symbol(maxid2)])

        return maxd, maxid1, maxid2

    def mdl_mol_block(self):
        """
        Returns a V3000 mol block of the molecule.

        Returns
        -------
        str
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
        atom_line = "M  V30 {0} {1} {2:.4f} {3:.4f} {4:.4f} 0\n"
        atom_block = ""

        # id bond_order atom1 atom2
        bond_line = "M  V30 {0} {1} {2} {3}\n"
        bond_block = ""

        main_string = main_string.format(self.mol.GetNumAtoms(),
                                         self.mol.GetNumBonds())
        # Kekulize the mol, which means that each aromatic bond is
        # converted to a single or double. This is necessary because
        # .mol V3000 only supports integer bonds.
        rdkit.Kekulize(self.mol)
        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_sym = periodic_table[atom.GetAtomicNum()]
            x, y, z = self.atom_coords(atom_id)
            atom_block += atom_line.format(atom_id+1,
                                           atom_sym, x, y, z)

        for bond in self.mol.GetBonds():
            bond_id = bond.GetIdx()
            atom1_id = bond.GetBeginAtomIdx() + 1
            atom2_id = bond.GetEndAtomIdx() + 1
            bond_order = int(bond.GetBondTypeAsDouble())
            bond_block += bond_line.format(bond_id, bond_order,
                                           atom1_id, atom2_id)

        main_string = main_string.replace(
                            "!!!ATOM!!!BLOCK!!!HERE!!!\n", atom_block)

        return main_string.replace(
                            "!!!BOND!!!BLOCK!!!HERE!!!\n", bond_block)

    def position_matrix(self):
        """
        Returns the position of all atoms in the as a matrix.

        Returns
        -------
        numpy.matrix
            The matrix is 3 x n. Each column holds the x, y and z
            coordinates of an atom. The index of the column corresponds
            to the id of the atom in the molecule.

        """

        pos_array = np.array([])
        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            pos_vect = np.array([*self.atom_coords(atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        return np.matrix(pos_array.reshape(-1, 3).T)

    def same(self, other):
        """
        Check if `other` has the same molecular structure.

        Parameters
        ----------
        other : MacroMolecule
            The ``MacroMolecule`` instance you are checking has
            the same structure.

        Returns
        -------
        bool
            Returns ``True`` if the building blocks and topology
            of the macromolecules are the same.

        """

        return self.inchi == other.inchi

    def rotate(self, theta, axis):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the molecular centroid.

        Parameters
        ----------
        theta : float
            The size of the rotation in radians.

        axis : numpy.array
            The axis about which the rotation happens.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The positions of all the atoms are changed due to the
            rotation.

        Returns
        -------
        None : NoneType

        """

        # Save the original position.
        og_position = self.centroid()
        # Move the centroid of the molecule to the origin, so that the
        # rotation occurs about this point.
        self.set_position([0, 0, 0])
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())
        # Apply the rotation.
        self.set_position_from_matrix(new_pos_mat)
        # Return the centroid of the molecule to the origin position.
        self.set_position(og_position)

    def save_bonders(self):
        """
        Adds atoms tagged with 'bonder' to `bonder_ids`.

        Modifies
        --------
        bonder_ids : list of ints
            Updates this set with the ids of atoms tagged 'bonder'.

        Returns
        -------
        None : NoneType

        """

        # Clear the set in case the method is run twice.
        self.bonder_ids = []
        for atom in self.mol.GetAtoms():
            if atom.HasProp('bonder'):
                self.bonder_ids.append(atom.GetIdx())

    def set_orientation(self, start, end):
        """
        Rotates the molecule by a rotation from `start` to `end`.

        Note: The difference between this method and
        :meth:`StructUnit._set_orientation2()` is about which point the
        rotation occurs: centroid of entire molecule versus centroid of
        bonder atoms, respectively.

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
        anything. This means that any geomteric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Parameters
        ----------
        start : numpy.array
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : numpy.array
            This array holds the vector, onto which `start` is rotated.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about its centroid.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule in `mol`.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)

        # Record the position of the molecule then translate the
        # centroid to the origin. This is so that the rotation occurs
        # about this point.
        og_center = self.centroid()
        self.set_position([0, 0, 0])

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())

        # Set the positions of the molecule.
        self.set_position_from_matrix(new_pos_mat)
        self.set_position(og_center)

        return self.mol

    def set_position(self, position):
        """
        Sets the centroid of the molecule to `position`.

        Parameters
        ----------
        position : numpy.array
            This array holds the position on which the centroid of the
            molecule should be placed.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed so that its
            centroid falls on `position`.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule with the centroid placed at `position`.
            This is the same instance as that in `mol`.

        """

        # Get the original centroid.
        centroid = self.centroid()
        # Find out how much it needs to shift to reach `position`.
        shift = position - centroid
        # Apply the shift and get the resulting rdkit conformer object.
        new_conf = self.shift(shift).GetConformer()

        # Replace the old rkdit conformer with one where the centroid
        # is at `position`.
        self.mol.RemoveAllConformers()
        self.mol.AddConformer(new_conf)

        return self.mol

    def set_position_from_matrix(self, pos_mat):
        """
        Set atomic positions of the molecule to those in `pos_mat`.

        Parameters
        ----------
        pos_mat : numpy.array
            The matrix holds the coordinates on which the atoms of the
            molecule should be placed.

            The dimensions are 3 x n. Each column of `pos_mat`
            represents the coordinates of a single atom. The 1st column
            sets the coordinate of the atom with id of 0. The next
            column sets the coordinate of the atom with id 1, and so
            on.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The coordinates of atoms in this molecule are set to the
            coordinates in `pos_mat`.

        Returns
        -------
        None : NoneType

        """

        conf = self.mol.GetConformer()
        for i, coord_mat in enumerate(pos_mat.T):
            coord = rdkit_geo.Point3D(coord_mat.item(0),
                                      coord_mat.item(1),
                                      coord_mat.item(2))
            conf.SetAtomPosition(i, coord)

    def shift(self, shift):
        """
        Shifts the coordinates of all atoms.

        This does not modify the molecule. A modified copy is returned.

        Parameters
        ----------
        shift : numpy.array
            A numpy array holding the value of the shift along each
            axis.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            A copy of the molecule where the coordinates have been
            shifted by `shift`.

        """

        # The function does not modify the existing conformer, as a
        # result a new instance is created and used for modification.
        conformer = rdkit.Conformer(self.mol.GetConformer())

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
            atom_position = np.array(
                                    conformer.GetAtomPosition(atom_id))

            # Inducing the shift.
            new_atom_position = atom_position + shift

            # Creating a new geometry instance.
            new_coords = rdkit_geo.Point3D(*new_atom_position)

            # Changes the position of the atom in the conformer to the
            # values stored in the new geometry instance.
            conformer.SetAtomPosition(atom_id, new_coords)

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
        new_mol.AddConformer(conformer)
        return new_mol

    def update_from_mae(self, mae_path):
        """
        Updates molecular structure to match an .mae file.

        Parameters
        ----------
        mae_path : str
            The full path of the .mae file from which the structure
            should be updated.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The rdkit molecule held in this attribute is changed so
            that it matches the moleclue held in the .mae file.

        Returns
        -------
        None : NoneType

        """

        self.mol = mol_from_mae_file(mae_path)

    def update_stereochemistry(self):
        """
        Updates stereochemistry tags on `mol` attribute.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The '_CIPCode' property on each atom of the rkdit molecule
            is updated to reflect the R or S configuration of the atom
            in the current geometry.

        Returns
        -------
        None : NoneType

        """

        for atom in self.mol.GetAtoms():
            atom.UpdatePropertyCache()
        rdkit.AssignAtomChiralTagsFromStructure(self.mol)
        rdkit.AssignStereochemistry(self.mol, True, True, True)

    def write(self, path):
        """
        Writes a molecular structure file of the molecule.

        This bypasses the need to use rdkit's writing functions, which
        have issues with macromolecules due to poor ring finding and
        sanitization issues.

        Parameters
        ----------
        path : str
            The `path` to which the molecule should be written.

        Returns
        -------
        None : NoneType

        """

        write_funcs = {'.mol': self._write_mdl_mol_file,
                       '.pdb': self._write_pdb_file}

        _, ext = os.path.splitext(path)
        write_func = write_funcs[ext]
        write_func(path)

    def _write_mdl_mol_file(self, path):
        """
        Writes a V3000 .mol file of the molecule

        This function should not be used directly, only via the
        ``write()`` method.

        Parameters
        ----------
        path : str
            The full path to which to the file should be written.

        Modifies
        --------
        A file is created/ changed at `path`.

        Returns
        -------
        None : NoneType

        """

        with open(path, 'w') as f:
            f.write(self.mdl_mol_block())

    def _write_pdb_file(self, path):
        """
        Writes a .pdb file of the molecule

        This function should not be used directly, only via the
        ``write()`` method.

        Parameters
        ----------
        path : str
            The full path to which to the file should be written.

        Modifies
        --------
        A file is created/ changed at `path`.

        Returns
        -------
        None : NoneType

        """

        # First write the file using rdkit.
        rdkit.MolToPDBFile(self.mol, path)

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
    Represents the building blocks of macromolecules examined by MMEA.

    ``Building blocks`` in this case refers to the smallest molecular
    units of the assembled molecules examined by MMEA.

    The goal of this class is to conveniently store information about,
    and perform operations on, single instances of macromolecular
    building blocks.

    An important part of this class is labelling atoms in functional
    groups.

    This class should only deal with issues that concern a single
    building block in and of itself.

    Class attributes
    ----------------
    init_funcs : dict of {str : function}
        This dictionary holds the various functions which can be used
        to initialize rdkit molecules and pairs them with the
        appropriate file extension.

    Attributes
    ----------
    file : str
        The full path of the molecular structure file holding the
        molecule. The supported file formats are the keys in the
        `init_funcs` dictionary. As long as a file of one of these
        types is provided, MMEA will automatically use the correct
        initialization function.

    mol : rdkit.Chem.rdchem.Mol
        The rdkit instance of the molecule held in `file`.

    inchi : str
        The InChI of the molecule.

    func_grp : FGInfo
        The ``FGInfo`` instance holding information about the
        functional group which will react when the building block
        assembles to form macromolecules.

    bonder_ids : list of ints
        A list holding the atom ids of the atoms which form bonds
        during macromolecular assembly.

    energy : Energy
        This attribute handles information about the energy of the
        instance. See the documentation of ``Energy`` to see what is
        available.

    optimized : bool (default = False)
        A flag to monitor whether an optimization has been performed on
        the molecule.

    key : MacroMolKey
        The key used for caching the molecule.

    name : str (default = "")
        A name which can be optionally given to the molcule for easy
        identification.

    note : str (default = "")
        A note or comment about the molecule. Purely optional but can
        be useful for labelling and debugging.

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
        Initializes a ``StructUnit`` instance.

        Parameters
        ----------
        file : str
            The full path of the molecular structure file holding the
            building block.

        functional_group : str (default = None)
            The name of the functional group which is to have atoms
            tagged. If ``None``, a functional group name found in the
            path `file`  is used. If no functional group is provided
            to this parameter and the name of one is not present in
            `file`, no tagging is done.

        note : str (default = "")
            A note or comment about the molecule.

        name : str (default = "")
            A name which can be optionally given to the molcule for
            easy identification.

        """

        self.file = file
        _, ext = os.path.splitext(file)

        if ext not in self.init_funcs:
            raise TypeError(
                   'Unable to initialize from "{}" files.'.format(ext))

        self.mol = self.init_funcs[ext](file)
        # Update the property cache of each atom. This updates things
        # like valence.
        for atom in self.mol.GetAtoms():
            atom.UpdatePropertyCache()

        super().__init__(name, note)

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

    @classmethod
    def init_random(cls, db, fg=None, name="", note=""):
        """
        Picks a random file from `db` to initialize from.

        Parameters
        ----------
        db : str
            A path to a database of molecular files.

        fg : str, optional
            The name of a functional group which the molecules in `db`
            have. By default it is assumed the name is present in the
            path of the files.

        name : str, optional
            The name to be given to the created molecule.

        note : str, optional
            A note to be given to the created molecule.

        Returns
        -------
        StructUnit
            A random molecule from `db`.

        None : NoneType
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

    def all_bonder_distances(self):
        """
        Yield distances between all pairs of bonder atoms.

        All distances are only yielded once. This means that if the
        distance between atoms with ids ``1`` and ``2``is yielded as
        ``(12.4, 1, 2)``, no tuple of the form ``(12.4, 2, 1)`` will be
        yielded.

        Yields
        ------
        tuple of form (int, int, scipy.double)
            The ints represnt the atoms ids and the double is their
            distance.

        """

        # Iterate through each pair of atoms - do not allow
        # recombinations.
        for atom1, atom2 in it.combinations(self.bonder_ids, 2):
                yield (atom1, atom2, self.atom_distance(atom1, atom2))

    def bonder_centroid(self):
        """
        Returns the centroid of the bonder atoms.

        Returns
        -------
        numpy.array
            A numpy array holding the midpoint of the bonder atoms.

        """

        centroid = sum(self.atom_coords(x) for x in self.bonder_ids)
        return np.divide(centroid, len(self.bonder_ids))

    def bonder_direction_vectors(self):
        """
        Yields the direction vectors between all pairs of bonder atoms.

        The yielded vector is normalized. If a pair (1,2) is yielded,
        the pair (2,1) will not be.

        Yields
        ------
        tuple of (int, int, numpy.array)
            The ints in the tuple represent the ids of the start and
            end atoms, respectively. The array is the direciont vector
            running between the atomic positions.

        """

        for atom1_id, atom2_id in it.combinations(self.bonder_ids, 2):
            p1 = self.atom_coords(atom1_id)
            p2 = self.atom_coords(atom2_id)

            yield atom2_id, atom1_id, normalize_vector(p1-p2)

    def bonder_position_matrix(self):
        """
        Returns a matrix holding the positions of bonder atoms.

        Returns
        -------
        numpy.matrix
            The matrix is 3 x n. Each column holds the x, y and z
            coordinates of a bonder atom. The index of the column
            corresponds to the index of the atom id in `bonder_ids`.

        """

        pos_array = np.array([])

        for atom_id in self.bonder_ids:
            pos_vect = np.array([*self.atom_coords(atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        return np.matrix(pos_array.reshape(-1, 3).T)

    def centroid_centroid_dir_vector(self):
        """
        Returns the direction vector between the 2 molecular centroids.

        The first molecule centroid is the centroid of the entire
        molecule. The second molecular centroid is the centroid of the
        bonder atoms.

        Returns
        -------
        numpy.array
            The normalized direction vector running from the centroid
            of the bonder atoms to the molecular centroid.

        """

        # If the bonder centroid and centroid are in the same position,
        # the centroid - centroid vector should be orthogonal to the
        # bonder direction vector.
        if np.allclose(self.centroid(),
                       self.bonder_centroid(), atol=1e-5):
            *_, bvec = next(self.bonder_direction_vectors())
            # Construct a secondary vector by finding the minimum
            # component of bvec and setting it to 0.
            vec2 = list(bvec)
            minc = min(vec2)
            vec2[vec2.index(min(vec2))] = 0 if abs(minc) >= 1e-5 else 1
            # Get a vector orthogonal to bvec and vec2.
            a = normalize_vector(np.cross(bvec, vec2))
            return a

        else:
            return normalize_vector(self.centroid() -
                                    self.bonder_centroid())

    def core(self):
        """
        Return the molecule with no H or functional group atoms.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The "core" of the molecule.

        """

        emol = rdkit.EditableMol(self.mol)
        for atom in reversed(self.mol.GetAtoms()):
            atomid = atom.GetIdx()
            if not self.is_core_atom(atomid):
                emol.RemoveAtom(atomid)
                continue
        return emol.GetMol()

    def functional_group_atoms(self):
        """
        Returns a container of atom ids of atoms in functional groups.

        Returns
        -------
        tuple of tuples of ints
            The form of the returned tuple is:
            ((1,2,3), (4,5,6), (7,8,9)). This means that all atoms with
            ids 1 to 9 are in a functional group and that the atoms 1,
            2 and 3 all form one functional group together. So do 4, 5
            and 5 and so on.

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
        Returns ``True`` if atom is not H or part of the fg.

        Parameters
        ----------
        atomid : int
            The id of the atom being queried.

        Returns
        -------
        bool
            Indicates whether the atom with `atomid` is part of the
            core.

        """

        # Tags are removed by multiprocessing - make sure to reapply
        # them.
        self.tag_atoms()
        atom = self.mol.GetAtomWithIdx(atomid)
        if atom.HasProp('fg') or atom.GetAtomicNum() == 1:
            return False
        return True

    def json(self):
        """
        Returns a JSON representation of the object.

        The representation has the following form:

            {
                'class' : 'StructUnit',
                'mol_block' : 'A string holding the V3000 mol
                               block of the molecule.',
                'note' : 'This molecule is nice.',
                'name' : 'benzene'
            }

        Returns
        -------
        dict
            A dict which represents the molecule.

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

        This function is not to be used. Use ``Molecule.load()``
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : dict
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : NoneType

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
        Generates the key for caching the molecule.

        Parameters
        ----------
        rdkit_mol : rdkit.Chem.rdchem.Mol
            An rdkit instance of the molecule.

        functional_group : str
            The name of the functional group being used to make
            macromolecules.

        Returns
        -------
        tuple
            The key used for caching the molecule.

        """

        return functional_group, rdkit.MolToInchi(rdkit_mol)

    def minimize_theta(self, v1, v2, axis, centroid):
        """
        Rotate mol to minimize angle between `v1` and `v2`.

        The rotation is done about the vector `axis`.

        Parameters
        ----------
        v1 : numpy.array
            The vector which is rotated.

        v2 : numpy.array
            The vector which is stationary.

        axis : numpy.array
            The vector about which the rotation happens.

        centroid : numpy.array
            The position vector at the center of the rotation.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to the
            rotation.

        Returns
        -------
        None : NoneType

        """

        # If the vector being rotated is not finite exit. This is
        # probably due to a planar molecule.
        if not all(np.isfinite(x) for x in v1):
            return

        # Save the initial position and change the origin to
        # `centroid`.
        iposition = self.centroid()
        self.mol = self.shift(-centroid)

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
            self.set_position(iposition)
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
        posmat = np.dot(rotmat, self.position_matrix())
        self.set_position_from_matrix(posmat)
        self.set_position(iposition)

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

    def rotate2(self, theta, axis):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the centroid of the bonder atoms.

        Parameters
        ----------
        theta : float
            The size of the rotation in radians.

        axis : numpy.array
            The axis about which rotation happens.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The atoms in this molecule are rotated.

        Returns
        -------
        None : NoneType

        """

        # Save the origin position of the bonder atom centroid.
        og_position = self.bonder_centroid()
        # Change the position of the centroid of the bonder atoms to
        # the origin so that the rotation occurs about this point.
        self.set_bonder_centroid([0, 0, 0])
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation on the original atomic coordinates to get
        # the new ones.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())
        # Set the atomic positions to the new coordinates.
        self.set_position_from_matrix(new_pos_mat)
        # Return the centroid to its original position.
        self.set_bonder_centroid(og_position)

    def set_bonder_centroid(self, position):
        """
        Move the molecule so that the bonder centroid is on `position`.

        Parameters
        ----------
        position : numpy.array
            A numpy array holding the desired the position. It holds
            the x, y and z coordinates, respectively.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The position of the molecule in this rdkit instance is
            changed, as described in this docstring.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule after it has been shifted. The same
            instance as held in `mol`.

        """

        center = self.bonder_centroid()
        shift = position - center
        new_conf = self.shift(shift).GetConformer()

        # Make sure the rkdit molecule has only one conformer.
        self.mol.RemoveAllConformers()
        self.mol.AddConformer(new_conf)

        return self.mol

    def _set_orientation2(self, start, end):
        """
        Note: The difference between this method and
        ``set_orientation()`` is about which point the rotation
        occurs: centroid of bonder atoms versus centroid of entire
        molecule, respectively.

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
        anything. This means that any geomteric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Parameters
        ----------
        start : numpy.array
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : numpy.array
            This array holds the vector, onto which `start` is rotated.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about its centroid.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule in `mol`.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)

        # Record the position of the molecule then translate the bonder
        # atom centroid to the origin. This is so that the rotation
        # occurs about this point.
        og_center = self.bonder_centroid()
        self.set_bonder_centroid(np.array([0, 0, 0]))

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())

        # Set the positions in the rdkit molecule.
        self.set_position_from_matrix(new_pos_mat)
        self.set_bonder_centroid(og_center)

        return self.mol

    def similar_molecules(self, database):
        """
        Returns molecules from `database` ordered by similarity.

        The most similar molecule is at index 0.

        This method uses the Morgan fingerprints of radius 4 to
        evaluate how similar the molecules in `database` are.

        Parameters
        ----------
        database : str
            The full path of the database from which the molecules are
            checked for similarity.

        Returns
        -------
        list of tuples of form (float, str)
            The float is the similarity of a given molecule in the
            database to `mol` while the str is the full path of the
            structure file of that molecule.

        """

        # First get the fingerprint of `self`.
        rdkit.GetSSSR(self.mol)
        self.mol.UpdatePropertyCache(strict=False)
        fp = rdkit.GetMorganFingerprint(self.mol, 4)

        # For every structure file in the database create a rdkit
        # molecule. Place these in a list.
        mols = []
        for file in os.listdir(database):
            path = os.path.join(database, file)
            # Ignore files which are not structure files and the
            # structure file of the molecule itself.
            _, ext = os.path.splitext(path)
            if ext not in self.init_funcs:
                continue

            mol = self.init_funcs[ext](path)
            rdkit.GetSSSR(mol)
            mol.UpdatePropertyCache(strict=False)
            mol_fp = rdkit.GetMorganFingerprint(mol, 4)
            similarity = DataStructs.DiceSimilarity(fp, mol_fp)
            mols.append((similarity, path))

        return sorted(mols, reverse=True)

    @classmethod
    def smarts_init(cls, smarts,
                    functional_group=None, note="", name=""):
        """
        Initialize from a SMARTS string.

        Parameters
        ----------
        smiles : str
            A SMARTS string of the molecule.

        functional_group : str (default = None)
            The name of the functional group which is to have atoms
            tagged. If no functional group is provided to this
            parameter, no tagging is done.

        note : str (default = "")
            A note or comment about the molecule.

        name : str (default = "")
            A name which can be optionally given to the molcule for
            easy identification.

        Returns
        -------
        StructUnit
            The StructUnit instance of the molecule represented by
            `smarts`.

        """

        mol = rdkit.MolFromSmarts(smarts)
        rdkit.SanitizeMol(mol)
        mol = rdkit.AddHs(mol)
        key = cls.gen_key(mol, functional_group)
        if key in cls.cache:
            return cls.cache[key]

        rdkit.EmbedMolecule(mol)
        obj = cls.__new__(cls)
        Molecule.__init__(obj, note, name)
        obj.file = smarts
        obj.key = key
        obj.mol = mol
        obj.func_grp = next((x for x in functional_groups if
                            x.name == functional_group), None)
        if obj.func_grp:
            obj.tag_atoms()

        cls.cache[key] = obj
        return obj

    def tag_atoms(self):
        """
        Adds bonding and deletion tags to atoms.

        All atoms which form the functional group of the molecule have
        the property 'fg' added. Its value is set to the name of the
        functional group.

        The atoms which form bonds during assembly have the property
        called 'bonder' added and set to '1'. Atoms which are deleted
        during reactions have the property 'del' set to '1'.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The atoms in this rdkit molecule have the properties 'fg',
            'bonder' and 'del' added, in accordance with the docstring.

        bonder_ids : set of ints
            Adds the ids of bonder atoms to this list.

        Returns
        -------
        None : NoneType

        """

        # Clear this list in case the method is being rerun.
        self.bonder_ids = []

        # Give all atoms in functional groups the tag 'fg' and set its
        # value to the name of the functional group.
        for atom_id in flatten(self.functional_group_atoms()):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetProp('fg', self.func_grp.name)

        # Give all atoms which form bonds during reactions the tag
        # 'bonder' and set its value to '1'. Add their ids to
        # `bonder_ids`.
        bond_mol = rdkit.MolFromSmarts(self.func_grp.bonder_smarts)
        bond_atoms = self.mol.GetSubstructMatches(bond_mol)
        for atom_id in flatten(bond_atoms):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetProp('bonder', '1')
            self.bonder_ids.append(atom_id)

        # Give all atoms which form bonds during reactions the tag
        # 'del' and set its value to '1'.
        del_mol = rdkit.MolFromSmarts(self.func_grp.del_smarts)
        del_atoms = self.mol.GetSubstructMatches(del_mol)
        for atom_id in flatten(del_atoms):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetProp('del', '1')

    def untag_atoms(self):
        """
        Removes the tags added by ``tag_atoms()``.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The atoms in this rdkit molecule have the properties
            'fg', 'bonder' and 'del' removed.

        bonder_ids : list of ints
            This list is set to [].

        Returns
        -------
        None : NoneType

        """

        self.bonder_ids = []

        for atom in self.mol.GetAtoms():
            atom.ClearProp('fg')
            atom.ClearProp('bonder')
            atom.ClearProp('del')

    def __str__(self):
        return "{} {}".format(self.__class__.__name__, list(self.key))

    def __repr__(self):
        return str(self)


class StructUnit2(StructUnit):
    """
    Represents building blocks with 2 functional groups.

    """

    def set_orientation2(self, end):
        """
        Rotate the molecule so that bonder atoms lie on `end`.

        The molecule is rotated about the centroid of the bonder atoms.
        It is rotated so that the direction vector running between the
        2 bonder atoms is aligned with the vector `end`.

        Parameters
        ----------
        end : numpy.array
            The vector with which the molecule's bonder atoms should be
            aligned.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the bonder
            atoms.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule in `mol`.

        """

        *_, start = next(self.bonder_direction_vectors())
        return self._set_orientation2(start, end)

    def minimize_theta2(self, vector, axis):
        """
        Rotates molecule about `axis` to minimze theta with `vector`.

        The molecule is rotated about `axis` passing through the bonder
        centroid. It is rotated so that the vector between the bonder
        and molecular centroids lies on the same plane as `vector`.

        Parameters
        ----------
        vector : numpy.array
            The vector to which the distance should be minimized.

        axis : numpy.array
            The direction vector along which the rotation happens.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the bonder
            atoms.

        Returns
        -------
        None : NoneType

        """

        self.minimize_theta(self.centroid_centroid_dir_vector(),
                            vector, axis, self.bonder_centroid())


class StructUnit3(StructUnit):
    """
    Represents building blocks with 3 functional groups.

    """

    def bonder_plane(self):
        """
        Returns the coefficients of the plane formed by bonder atoms.

        A plane is defined by the scalar plane equation,

            ax + by + cz = d.

        This method returns the a, b, c and d coefficients of this
        equation for the plane formed by the bonder atoms. The
        coefficents a, b and c decribe the normal vector to the plane.
        The coefficent d is found by substituting these coefficients
        along with the x, y and z variables in the scalar equation and
        solving for d. The variables x, y and z are substituted by the
        coordinate of some point on the plane. For example, the
        position of one of the bonder atoms.

        Returns
        -------
        numpy.array
            This array has the form [a, b, c, d] and represents the
            scalar equation of the plane formed by the bonder atoms.

        References
        ----------
        https://tinyurl.com/okpqv6

        """

        bonder_coord = self.atom_coords(self.bonder_ids[0])
        d = -np.sum(self.bonder_plane_normal() * bonder_coord)
        return np.append(self.bonder_plane_normal(), d)

    def bonder_plane_normal(self):
        """
        Returns the normal vector to the plane formed by bonder atoms.

        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the molecule.

        Returns
        -------
        numpy.array
            A unit vector which describes the normal to the plane of
            the bonder atoms.

        """

        if sum(1 for _ in self.bonder_direction_vectors()) < 2:
            raise ValueError(("StructUnit3 molecule "
                             "has fewer than 3 functional groups."))

        vgen = (v for *_, v in self.bonder_direction_vectors())
        v1, v2 = it.islice(vgen, 2)

        normal_v = normalize_vector(np.cross(v1, v2))

        theta = vector_theta(normal_v,
                             self.centroid_centroid_dir_vector())

        if theta > np.pi/2:
            normal_v *= -1

        return normal_v

    def minimize_theta2(self, atom, vector, axis):
        """
        Rotates molecule to minimize angle between `atom` and `vector`.

        The rotation is done about `axis` and is centered on the
        bonder centroid.

        Parameters
        ----------
        atom : int
            The id of atom which is to have angle minimized.

        vector : numpy.array
            A vector with which the angle is minimized.

        axis : numpy.array
            The vector about which the rotation happens.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The rdkit molecule is changed to reflect the rotation.

        Returns
        -------
        None : NoneType

        """

        v1 = self.atom_coords(atom) - self.bonder_centroid()
        self.minimize_theta(v1, vector, axis, self.bonder_centroid())

    def set_orientation2(self, end):
        """
        Rotates the molecule so the plane normal is aligned with `end`.

        Here ``plane normal`` referes to the normal of the plane formed
        by the bonder atoms in the substituted molecule. The molecule
        is rotated about the centroid of the bonder atoms. The rotation
        results in the normal of their plane being aligned with `end`.

        Parameters
        ----------
        end : numpy.array
            The vector with which the normal of plane of bonder atoms
            shoould be aligned.

        Modifies
        --------
        mol : rdkit.Chem.rdchem.Mol
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the bonder
            atoms.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit molecule in `mol`.

        """

        start = self.bonder_plane_normal()
        return self._set_orientation2(start, end)


@total_ordering
class MacroMolecule(Molecule, metaclass=Cached):
    """
    A class for assembled macromolecules.

    The goal of this class is to represent an individual used by the
    GA. As such, it holds attributes that are to be expected for this
    purpose. Mainly, it has a fitness value stored in its `fitness`
    attribute and a genetic code - as defined by its `building_blocks`
    and  `topology` attributes. If a change is made to either of these
    attributes, they should describe a different macromolecule. On the
    other hand, the same attributes should always describe the same
    macromolecule.

    Because of this, as well as the computational cost associated with
    macromolecule initialization, instances of this class are cached.
    This means that providing the same arguments to the initializer
    will not build a different instance with the same attribute values.
    It will yield the original instance, retrieved from memory.

    To prevent bloating this class, any information that can be
    categorized is. For example, storing information that concerns
    building blocks ``in a vacuum`` is stored in ``StructUnit``
    instances. Equally, manipulations of such data is also performed by
    those instances. Similarly, anything to do with topolgy should be
    held by a ``Topology`` instance in the topology attribute. There is
    a notable exception to this however. This happens when retrieving
    topological information directly from rdkit molecule instances
    of the macromolecule. Examples include the information about atomic
    coordinates, which can be access with the ``atom_coords()`` method.

    It should also be noted that only a single copy of each
    ``StructUnit`` instance representing a specific building block
    needs to be held. How many of such building blocks are need to
    assemble the molecule is the handled by the ``Topology`` class.

    This class is not intended to be used directly but should be
    inherited by subclasses representing specific macromolecules. The
    ``Cage`` and ``Polymer`` classes are examples of this. Any
    information or methods that apply generally to all macromolecules
    should be defined within this class while specific non-general data
    should be included in derived classes.

    This class also supports comparison operations, these act on the
    fitness value assiciated with a macromolecule. Comparison
    operations not explicitly defined are included via the
    ``total_ordering`` decorator. For other operations and methods
    supported by this class examine the rest of the class definition.

    Finally, a word of caution. The equality operator ``==`` compares
    fitness values. This means two macromolecules, made from different
    building blocks, can compare equal if they happen to have the same
    fitness. The operator is not to be used to check if one
    macromolecule is the same structurally as another. To do this check
    use the ``same()`` method. This method may be overwritten in
    derived classes, as necessary. In addition the ``is`` operator is
    implemented as is default in Python. It compares whether two
    objects are in the same location in memory. Because the
    ``MacroMolecule`` class is cached the ``is`` operator could in
    principle be used instead of the `same` method (including in
    derived classes). However, this is not intended use and is not
    guaranteed to work in future implementations. If caching stops
    being implemented such code would break.

    Attributes
    ----------
    building_blocks : list of ``StructUnit`` instances
        This attribute holds ``StructUnit`` instances which represent
        the monomers forming the macromolecule. Only one ``StructUnit``
        instance is needed per building block, even if multiples of a
        building block join up to form the macromolecule

    bb_counter : Counter
        A counter keeping track of how much of each building block is
        used to form the macromolecule.

    topology : Topology
        An object of a class derived from ``Topology``. It
        assembles the marcormolecule from the building_blocks.

    mol : rdkit.Chem.rdchem.Mol
        An rdkit instance representing the macromolecule.

    inchi : str
        The InChI of the molecule.

    optimized : bool (default = False)
        This is a flag to indicate if a molecule has been previously
        optimized. Optimization functions set this flag to ``True``
        after an optimization.

    energy : Energy
        This attribute handles information about the energy of the
        instance. See the documentation of ``Energy`` to see what is
        available.

    fitness : float (default = None)
        The fitness value of the macromolecule, as determined by the
        chosen fitness functions together with any normalization
        functions. Fitness functions do not place values into this
        attribute directly.

    unscaled_fitness : dict
        The dictionary hold the name of a fitness function as the
        key and the value it calculated for unscaled_fitness as the
        value.

    progress_params : dict (default = {})
        Holds the fitness parameters which the GA should track to make
        progress plots. The key is the name of a fitness function.

    bonds_made : int
        The number of bonds created during assembly.

    bonder_ids : list of ints
        The ids of atoms which have bonds added during assembly. This
        list is sorted from lowest to highest id.

    fg_ids : set of ints
        The ids of atoms which were parth of the functional group of
        the building blocks.

    key : MacroMolKey
        The key used for caching the molecule. Necessary for
        `update_cache` to work. This attribute is assigned by the
        `__call__()` method of the ``Cached`` metaclass.

    name : str (default = "")
        A name which can be optionally given to the molcule for easy
        identification.

    note : str (default = "")
        A note or comment about the molecule. Purely optional but can
        be useful for labelling and debugging.

    fragments : dict
        The key in this dictionary is a tuple of form (int, int). The
        first int is the index of a building block within
        `building_blocks`. The second int indentifies a molecule of
        that building block. For example, (1, 3) identifies the 4th
        molecule of type building_blocks[1] to be added to the
        macromolecule during assembly.

        The value in the dictionary is a set of ints. These hold the
        atom ids belonging to a particular molecule before assembly.

    """

    def __init__(self, building_blocks, topology, name="", note=""):
        """
        Initialize a ``MacroMolecule`` instance.

        When an exception occurs during initialization, all parameters
        which were provided to the initializer are saved to a file
        ``failures.txt`` which is in the ``output`` folder.

        Parameters
        ---------
        building_blocks : iterable of ``StructUnit`` instances
            A set of ``StructUnit`` instances which represent the
            monomers forming the macromolecule.

        topology : Topology
            An instance of a class derived from ``Topology``. It
            assembles the marcormolecule from the building_blocks.

        note : str (default = "")
            A note or comment about the molecule.

        name : str (default = "")
            A name which can be optionally given to the molcule for
            easy identification.


        """

        self.fitness = None
        self.unscaled_fitness = {}
        self.progress_params = {}
        self.building_blocks = building_blocks
        self.bb_counter = Counter()
        self.topology = topology
        self.fg_ids = set()
        self.bonds_made = 0
        self.fragments = defaultdict(set)

        try:
            # Ask the ``Topology`` instance to assemble/build the
            # macromolecule. This creates the `mol` attribute.
            topology.build(self)

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
            for bb in building_blocks:
                bb_blocks.append(
                    ('{0.__class__.__name__} {0.func_grp.name}\n'
                     '{1}').format(bb, bb.mdl_mol_block()))

            errormsg += '\n'.join(bb_blocks)

            logger.error(errormsg, exc_info=True)

        super().__init__(name, note)
        self.save_ids()

    def building_block_cores(self, bb):
        """
        Yields the "cores" of a building block molecule.

        The structure of the cores are representative of how they are
        found in the macromolecule.

        Parameters
        ----------
        bb : int
            The index of a building block molecule within
            `building_blocks`. The cores of this molecule are yielded.

        Yields
        ------
        rdkit.Chem.rdchem.Mol
            The core of a building block molecule, as found in the
            macromolecule.

        """

        for (bb_index, mol_index), atoms in self.fragments.items():
            # Ignore fragments which do not correspond to the molecule
            # `bb`.
            if bb_index != bb:
                continue

            # For each fragment make a new core. To preserver bond
            # information, start with the macromolecule.
            coremol = rdkit.EditableMol(self.mol)
            # Remove any atoms not part of the core - atoms not in the
            # fragment itself, hydrogens, atoms in the functional
            # group.
            for atom in reversed(self.mol.GetAtoms()):
                atomid = atom.GetIdx()
                if (atomid not in atoms or
                    atomid in self.fg_ids or
                   atom.GetAtomicNum() == 1):
                    coremol.RemoveAtom(atomid)
            yield coremol.GetMol()

    def json(self):
        """
        Returns a JSON representation of the object.

        The representation has the following form:

            {
                'class' : 'Polymer',
                'mol_block' : 'A string holding the V3000 mol
                               block of the molecule.'
                'building_blocks' : {bb1.json(), bb2.json()}
                'topology' : 'Copolymer(repeating_unit='AB')'
                'unscaled_fitness' : {'fitness_func1' : fitness1,
                                      'fitness_func2' : fitness2},
                'note' : 'A nice molecule.',
                'name' : 'Poly-Benzene'
            }

        Returns
        -------
        dict
            A dict which represents the molecule.

        """

        return {
            'bb_counter': [(key.json(), val) for key, val in
                           self.bb_counter.items()],
            'bonds_made': self.bonds_made,
            'bonder_ids': self.bonder_ids,
            'fg_ids': list(self.fg_ids),
            'class': self.__class__.__name__,
            'mol_block': self.mdl_mol_block(),
            'building_blocks': [x.json() for x in
                                self.building_blocks],
            'topology': repr(self.topology),
            'unscaled_fitness': repr(self.unscaled_fitness),
            'progress_params': self.progress_params,
            'note': self.note,
            'name': self.name,
            'fragments': dict(zip(
                          (str(x) for x in self.fragments.keys()),
                          (list(x) for x in self.fragments.values())))

        }

    @classmethod
    def _json_init(cls, json_dict):
        """
        Completes a JSON initialization.

        This function is not to be used. Use ``Molecule.load()``
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : dict
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : NoneType

        """

        bbs = [Molecule.fromdict(x) for x in
               json_dict['building_blocks']]

        topology = eval(json_dict['topology'],  topologies.__dict__)

        key = cls.gen_key(bbs, topology)
        if key in cls.cache:
            return cls.cache[key]

        obj = cls.__new__(cls)
        obj.mol = rdkit.MolFromMolBlock(json_dict['mol_block'],
                                        sanitize=False, removeHs=False)
        obj.topology = topology
        obj.unscaled_fitness = eval(json_dict['unscaled_fitness'],
                                    np.__dict__)
        obj.fitness = None
        obj.progress_params = json_dict['progress_params']
        obj.bb_counter = Counter({Molecule.fromdict(key): val for
                                  key, val in json_dict['bb_counter']})
        obj.bonds_made = json_dict['bonds_made']
        obj.energy = Energy(obj)
        obj.bonder_ids = json_dict['bonder_ids']
        obj.fg_ids = set(json_dict['fg_ids'])
        obj.optimized = json_dict['optimized']
        obj.note = json_dict['note']
        obj.name = json_dict['name'] if json_dict['load_names'] else ""
        obj.key = key
        obj.building_blocks = bbs
        obj.fragments = defaultdict(set)
        obj.fragments.update(zip(
            (tuple(int(y) for y in x.strip('()').split(','))
             for x in json_dict['fragments'].keys()),
            (set(x) for x in json_dict['fragments'].values())
        ))
        cls.cache[key] = obj
        return obj

    @staticmethod
    def gen_key(building_blocks, topology):
        """
        Generates the key for caching the molecule.

        Parameters
        ----------
        building_blocks : iterable of StructUnit instances
            The building blocks used to make the macromolecule.

        topology : Topology
            The topology used to make the macromolecule.

        Returns
        -------
        tuple
            The key used for caching the macromolecule.

        """

        return (frozenset(x.key for x in building_blocks),
                repr(topology))

    def bb_distortion(self):
        """
        Rmsd difference of building blocks before and after assembly.

        The function looks at each building block in the macromolecule
        and calculates the rmsd between the "free" verson and the one
        present in the macromolecule. The mean of these rmsds is
        returned.

        Atoms which form the functional group of the building blocks
        and hydrogens are excluded from the calculation.

        Returns
        -------
        float
            The mean rmsd of the macromole's building blocks to their
            "free" counterparts.

        """

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
                rmsd += rdkit.AlignMol(free, frag, atomMap=am)
                n += 1
        return rmsd / n

    def save_ids(self):
        """
        Updates `bonder_ids` and `fg_ids` attributes.

        Modifies
        --------
        bonder_ids : list of ints
            All molecules tagged 'bonder' have their ids added to this
            list.

        fg_ids : set of ints
            All molecules tagged 'fg' have their ids add to
            this set.

        Returns
        -------
        None : NoneType

        """

        self.bonder_ids = []
        self.fg_ids = set()
        for atom in self.mol.GetAtoms():
            if atom.HasProp('bonder'):
                self.bonder_ids.append(atom.GetIdx())
            if atom.HasProp('fg'):
                self.fg_ids.add(atom.GetIdx())
        self.bonder_ids = sorted(self.bonder_ids)

    def update_cache(self):
        """
        Set cached molecule with same 'key' to have equal attributes.

        When an instance of ``MacroMolecule`` is first created it
        is cached. Using multiprocessing to perform optimizations or
        calculate fitness returns modified copies of the cached
        molecules. In order to ensure that the cached molecules have
        their attributes updated to the values of the copies, this
        method must be run on the copies.

        Returns
        -------
        None : NoneType

        """

        self.__class__.cache[self.key].__dict__ = dict(vars(self))

    def update_fragments(self):
        """
        Saves rdkit atom properties in `fragments`.

        The properties saved are 'bb_index' and 'mol_index'. These
        should be added to the atoms by the ``place_mols()`` function
        during ``build()``.

        Modifies
        --------
        fragments : dict
            The dictionary is updated with the data from the 'bb_index'
            and 'mol_index' properties.

        Returns
        -------
        None : NoneType

        """

        self.fragments = defaultdict(set)
        for atom in self.mol.GetAtoms():
            if not atom.HasProp('bb_index'):
                continue

            molkey = (int(atom.GetProp('bb_index')),
                      int(atom.GetProp('mol_index')))
            self.fragments[molkey].add(atom.GetIdx())

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

    def window_difference(self):
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

        Returns
        -------
        float
            The total difference of window size when considering
            every combination of windows of the same type.

        None : NoneType
            If not all windows were found.


        Raises
        ------
        WindowError
            When the number of found windows is less than the
            number of expected windows. Likely due to a collapsed
            cage.

        """

        if (self.windows is None or
           len(self.windows) < self.topology.n_windows):
            return None

        # Cluster the windows into groups so that only size
        # differences between windows of the same type are taken
        # into account. To do this, first sort the windows by
        # size. If two windows types are present split the
        # windows at the two groups at the point where the window
        # sizes have the biggest difference. If there are three
        # types split it at the two biggest differences and so
        # on.
        windows = np.array(self.windows)

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

    @property
    def windows(self):
        """
        Returns window sizes found by pyWindow.

        Returns
        -------
        None : NoneType
            If the function for finding windows and their sizes
            found fewer than the required number of windows or
            if it failed for some other reason.

        list of floats
            Each float in the list represents the size of a
            window in the cage. If the window finding function
            found more than the expected number of windows, only
            the largest n windows are returned. Where n is the
            number of expected windows.

        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Load an RDKit molecule object to pyWINDOW.
            pw_molecule = pywindow.molecular.Molecule.load_rdkit_mol(self.mol)
            # Find windows and get a single array with windows' sizes.
            all_windows = pw_molecule.calculate_windows(output='windows')

        # If pyWindow failed, return ``None``.
        if all_windows is None:
            return None

        all_windows = sorted(all_windows,
                             reverse=True)[:self.topology.n_windows]
        for i, x in enumerate(all_windows):
            # Return ``None`` when pyWindow fucks up and outputs
            # a mistakenly large window size.
            if x > 500:
                return None

        return all_windows


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

    The method :meth:`~.MacroMolecule.save_ids` is extended to take
    into account ids of atoms in perdiodic bonds.

    Attributes
    ----------
    terminator_coords : :class:`dict`
        The key is an :class:`int` representing the id of a bonder
        atom. The value is a :class:`numpy.ndarray` which holds the x,
        y and z coordinates of a deleter atom attached to the bonder.
        The coordinates are relative to the bonder atom.

    periodic_bonds : :class:`list` of :class:`.PeriodicBond`
        When periodic topologies are assembled, periodic bonds
        do not get added to the rdkit molecule in the
        :attr:`~.MacroMolecule.mol` attribute. Instead,
        :meth:`~.PeriodicLattice.join_mols` adds
        :class:`.PeriodicBond` instances representing the bonds into
        this list.

    _ids_updated : :class:`bool`
        Indicates whether periodic bond ids have been updated already
        by :meth:`save_ids`.

    """

    def __init__(self, building_blocks, topology, name="", note=""):
        self.periodic_bonds = []
        self._ids_updated = False
        super().__init__(building_blocks, topology, name="", note="")

    def _is_subterminal(self, atom_id, bonder_map, bonded):
        """
        ``True`` if atom is bonded to a terminal one.

        Notes
        -----
        For internal use by :meth:`island`.

        Parameters
        ----------
        atom_id : :class:`int`
            The id of the atom being checked.

        bonder_map : :class:`dict`
            A dictionary which maps the atom id in the island to the
            equivalent bonder atom in the unit cell.

        bonded : :class:`set`
            A set of atom ids of all bonder atoms which have had a
            bond added by :meth:`_join_island`.

        Returns
        -------
        :class:`bool`
            ``True`` if atom is bonded to a terminal atom, ``False``
            otherwise.

        """

        # If `atom_id` is not in the `bonder_map` it means that the
        # atom is not a bonder and therefore no terminal atoms need to
        # be added to it.
        if atom_id not in bonder_map:
            return False

        # Get id of the equivalent bonder atom in the original unit
        # cell.
        bid = bonder_map[atom_id]
        # Make a set of all bonder atoms in the original unit cell
        # which were registered as having bonds crossing periodic
        # boundaries.
        periodic = {x.atom1 for x in self.periodic_bonds}
        periodic.update(x.atom2 for x in self.periodic_bonds)
        # If that atom does have a periodic bond and has not had a
        # bond added while building the island - it is subterminal and
        # needs to have a terminal atom attached.
        if bid not in periodic or atom_id in bonded:
            return False

        return True

    def island(self, dimensions, terminator=1, bond_type='1'):
        """
        Build a terminated supercell.

        Terminated means that the periodic bonds are replaced with
        bonds to terminating atoms.

        Parameters
        ----------
        dimensions : :class:`list` of :class:`int`
            A 3 member :class:`list`, holding the number of unit cells
            in the x, y and z directions used to make the supercell.

        terminator : :class:`int`, optional
            The atomic number of the terminating atoms added to the
            supercell.

        bond_type : :class:`str`, optional
            The bond type of bonds to terminating atoms. Valid options
            include ``'1'``, ``'2'``, ``'3'`` and ``'ar'``.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            An rdkit molecule of the island.

        """

        cells, island, bonder_map = self._place_island(dimensions)
        island, bonded = self._join_island(cells, island)
        return self._terminate_island(island, bonded, bonder_map,
                                      terminator, bond_dict[bond_type])

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
        :class:`tuple`
            The first member of the tuple an rdkit molecule of the
            island with all bonds between the unit cells added. The
            second member of the tuple is a :class:`set` of
            :class:`int`, where each :class:`int` is id of a bonder
            atom which had a bond added by this function.

        """

        emol = rdkit.EditableMol(island)
        bonded = set()
        # `self.periodic_bonds` holds objects of the
        # ``PeriodicBond`` class. Each ``PeriodicBond`` object has the
        # ids of two bonder atoms in the unit cell which are connected
        # by a bond running across the periodic boundary. The
        # `direction` attribute descibes the axes along which the
        # bond is periodic. For example, if `direction1` is [1, 0, 0]
        # it means that the bonder atom in `periodic_bond.atom1` has a
        # perdiodic bond connecting it to `periodic_bond.atom2` going
        # in the positive direction along the x-axis.

        # When iterating through all the unit cells composing the
        # island, you can use the `direction` vector to get index of
        # the unit cell which holds atoms connected the present cell.
        # Then just form bonds between the correct atoms by mapping
        # the atom ids in the unit cells to the ids of the equivalent
        # atoms in the original unit cell  and checking the
        # `periodic_bond` to see which atom ids are connected.
        for cell in flatten(cells):
            for periodic_bond in self.periodic_bonds:
                # Get the indices of the cell which holds the atom
                # bonded to the equivalent atom of
                # `periodic_bond.atom1` in the present `cell`.
                x, y, z = cell.id + periodic_bond.direction
                try:
                    # ccel as in "connected cell".
                    ccell = cells[x][y][z]
                except:
                    continue

                # `bonder1` is the id of a bonder atom, found in `cell`
                # and equivalent to `periodic_bond.atom1`, having a
                # bond added.
                bonder1 = cell.bonders[periodic_bond.atom1]
                # `bonder2` is the id of a bonder atom, found in
                # `ccell` and equivalent to `periodic_bond.atom2`,
                # having a bond added.
                bonder2 = ccell.bonders[periodic_bond.atom2]
                bond_type = self.topology.determine_bond_type(
                              self,
                              periodic_bond.atom1,
                              periodic_bond.atom2)
                emol.AddBond(bonder1, bonder2, bond_type)
                bonded.add(bonder1)
                bonded.add(bonder2)

        return emol.GetMol(), bonded

    def _terminate_island(self, island, bonded,
                          bonder_map, terminator, bond_type):
        """
        Adds atoms to all bonder atoms at the edges of `island`.

        Notes
        -----
        For internal use by :meth:`island`.

        Parameters
        ----------
        island : :class:`rdkit.Chem.rdchem.Mol`
            The island molecule which needs to have terminating atoms
            added.

        bonded : :class:`set` of :class:`int`
            Contains all the ids of all bonder atoms in `island` which
            have already had a bonded added by :methd:`_join_island`.

        bonder_map : :class:`dict`
            This is a mapping of the ids of bonder atoms in the island
            back to the id of the equivalent atom in the original unit
            cell.

        terminator : :class:`int`
            The atomic number of the terminating atoms added to the
            supercell.

        bond_type : :class:`str`
            A string holding the bond type of bonds to the terminating
            atoms. Valid options include ``'1'``, ``'2'``, ``'3'``
            and ``'ar'``.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The rdkit molecule of the final, assembled island.

        """

        iconf = island.GetConformer()
        emol = rdkit.EditableMol(island)
        coords = {}
        for atom in island.GetAtoms():
            atom_id = atom.GetIdx()
            # If the atom is not connected to a terminal atom, move on
            # to the next one.
            if not self._is_subterminal(atom_id, bonder_map, bonded):
                continue
            tid = emol.AddAtom(rdkit.Atom(terminator))
            emol.AddBond(tid, atom_id, bond_type)
            # Get the id of bonder atom in the original unit cell
            # which is equivalent to `atom`.
            bi = bonder_map[atom_id]
            # Using the equivalent bonder atom get the position
            # of the terminating atom relative to the bonder atom. Add
            # the relative position of the terminating atom and the
            # position of `atom` to get the final position of the
            # terminating atom.
            tcoords = (np.array(iconf.GetAtomPosition(atom_id)) +
                       self.terminator_coords[bi])
            rdkit_coords = rdkit_geo.Point3D(*tcoords)
            coords[tid] = rdkit_coords

        mol = emol.GetMol()
        conf = mol.GetConformer()
        # Update the positions of all the terminating atoms in the
        # conformer.
        for atom_id, atom_coords in coords.items():
            conf.SetAtomPosition(atom_id, atom_coords)

        return mol

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

            The second member is an rdkit molecule of the island being
            built. The third member is a :class:`dict` mapping the
            ids of bonder atoms in the island back to the id of the
            equivalent atom in the original unit cell.

        """

        a, b, c = self.topology.cell_dimensions
        cells = np.full(dimensions, None, object).tolist()
        island = rdkit.Mol()
        bonder_map = ChainMap()
        i = 0
        # For each dimension place a unit cell.
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    island = rdkit.CombineMols(
                                island, self.shift(x*a + y*b + z*c))
                    # `bonders` maps a bonder id in the original unit
                    # cell to the one currently being added to the
                    # island.
                    bonders = {bi: i*self.mol.GetNumAtoms() + bi for
                               bi in self.bonder_ids}
                    cells[x][y][z] = Cell((x, y, z), bonders)
                    # 'bonder_map' maps all bonder ids in the island
                    # to the bonder ids in the original unit cell.
                    inverse_bonders = dict(zip(bonders.values(),
                                               bonders.keys()))
                    bonder_map = bonder_map.new_child(inverse_bonders)
                    i += 1
        return cells, island, bonder_map

    def save_ids(self):
        super().save_ids()

        # If the ids of periodic bonds have already been updated, stop.
        if self._ids_updated:
            return

        # Update periodic bonds to hold atom ids directly instead of
        # indices of the atom ids within `bonder_ids`.
        for pb in self.periodic_bonds:
            pb.atom1 = self.bonder_ids[pb.atom1]
            pb.atom2 = self.bonder_ids[pb.atom2]

        # Update terminator_coords to hold atom ids directly instead
        # of indices of the atom ids within `bonder_ids`.
        terminator_coords = {}
        for index, coord in self.terminator_coords.items():
            bonder_id = self.bonder_ids[index]
            terminator_coords[bonder_id] = coord
        self.terminator_coords = terminator_coords

        self._ids_updated = True

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
            # Write the cell paramters.
            f.write('cell\n')
            # The sizes of cell vectors a, b and c are written first.
            for vector in self.topology.cell_dimensions:
                f.write(str(np.round(np.linalg.norm(vector), 6)) + ' ')
            # Then angles alpha, beta and gamma.
            a, b, c = self.topology.cell_dimensions
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
                a1 = bond.atom1 + 1
                a2 = bond.atom2 + 1
                dx, dy, dz = bond.direction
                f.write('connect {} {} {:+} {:+} {:+}\n'.format(a1, a2,
                                                                dx, dy,
                                                                dz))
