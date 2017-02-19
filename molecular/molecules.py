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
type of monomer form a macromolecule, only 2 StructUnit instances are in
`building_blocks`.

In its `topology` attribute a MacroMolecule holds a ``Topology` child
class instance. This instance is responsible for assembling the
macromolecule from the building blocks. The building should happen in
the ``__init__()`` method of the MacroMolecule via the Topology
instance's ``build()`` method. The ``build()`` method should place the
assembled macromoleclue in the `mol` attribute of the MacroMolecule
instance.

The StructUnit class labels atoms in the functional groups of the
building blocks as either ``bonder`` or ``del`` (see its documentation).
This tells the Topology intance which atoms form bonds and which are
removed during assembly.

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
The module ``/mmea/molecular_tools/fg_info.py`` defines a class called
FGInfo. Instances of this class are held in the list `functional_groups`
(also in that module). If you put an FGInfo instance in that list, the
functional group will be recognized.

    3) Place the FGInfo instance of the functional group in the
       `func_grp` attribute of the StructUnit.

    4) Using the FGInfo instance, tag atoms in the building block as
       either 'bonder' or 'del'. 'bonder' signifies that the atoms form
       a bond during macromolecular assembly, while 'del' means they are
       deleted.

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
       building blocks are orientated correctly when forming the
       macromolecule.

    10) Create bonds between all the disjoined building block molecules.
        This is where the tagging done by StructUnit is needed.

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

If you're adding a new class of macromolecules, it quite likely you want
to add a new ``Topology`` class. See the docstring of the
/mmea/molecular_tools/topologies/base.py module for guidance on adding
these. The topology class does the assembly of the macromolecule from
the building blocks. The assembled macromolecule instance is then held
by the macromolecular class.

"""

import numpy as np
from functools import total_ordering, partial
import itertools as it
from rdkit import Chem as chem
from rdkit.Chem import AllChem as ac
import rdkit.Geometry.rdGeometry as rdkit_geo
from rdkit import DataStructs
import os
import networkx as nx
from scipy.spatial.distance import euclidean
import pickle
import json
from uuid import uuid4

from ..convenience_tools import (flatten, periodic_table, MolError,
                                 normalize_vector, rotation_matrix,
                                 vector_theta, mol_from_mae_file,
                                 rotation_matrix_arbitrary_axis)
from .fg_info import functional_groups
from .energy import Energy

class CachedMacroMol(type):
    """
    A metaclass for creating classes which create cached instances.

    This class is tailored to the needs of createding cached
    ``MacroMolecule`` instances.

    Extending MMEA
    --------------
    If a MacroMolecule class is added such that one of its initializer
    args of kwargs should not be used for caching, the `__call__()`
    method in this class will need to be modified.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache = dict()

    def __call__(self, *args, **kwargs):
        # Sort the first argument which corresponds to an iterable of
        # StructUnit instances. You want (bb1, lk1) to register as the
        # same as (lk1, bb1). Because this creates the same cage.
        _, *other_args = args
        args = [sorted(args[0])]
        args.extend(other_args)
        # Do not take the last arg because this is the name of the
        # ``.mol`` file. This will likely be different even if the cage
        # is the same. However even in this case you want to return the
        # cached instance.
        key = str(args[:-1]) + str(kwargs)
        if key in self._cache.keys():
            return self._cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key

            dmp_file = os.path.splitext(obj.file)[0] + '.dmp'
            obj.dump(dmp_file)
            self._cache[key] = obj
            return obj

    def _update_cache(self, macro_mol):
        """
        Updates the cache of stored molecule.

        Parallel processes, such as optimization, return a copy of the
        optimized molecule. The original molecule, stored in `_cache`,
        does not have its attributes updated. In order to replace the
        molecule with the optimized copy within the cache this function
        should be used.

        Parameters
        ----------
        population : iterable of MacroMolecule instances
            A population holding molecules which should replace the
            ones held in `_cache` which share the same key.

        Returns
        -------
        None : NoneType

        """

        self._cache[macro_mol.key] = macro_mol


class Cached(type):
    """
    A metaclass for creating classes which create cached instances.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__cache = dict()

    def __call__(self, *args, **kwargs):
        key = str(args) + str(kwargs)
        if key in self.__cache.keys():
            return self.__cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            self.__cache[key] = obj
            return obj

class Molecule:
    """
    The most basic class representing molecules.

    This class defines defines the operations which any class describing
    molecules should inherit or may find useful. Examples of such
    classes are ``StructUnit`` and ``MacroMolecule``. This calls is
    unlikely to be useful as in and of itself. It lacks an `__init__()`
    method because it does not need one. Child classes should define it.

    It is assumed that child classes will have some basic attributes.

    Minimum required attributes of child classes
    --------------------------------------------
    mol : rdkit.Chem.rdchem.Mol
        A rdkit molecule instance representing the molecule.

    file : str
        The full path of the molecular structure file of the molecule.

    energy : Energy
        An instance of the ``Energy`` class. It handles all things
        energy.

    optimized : bool
        Indicates whether a Molecule has been passed through an
        optimization function or not.

    failed : bool
        True if the fitness function or optimization function failed
        when trying to evaluate the molecule.

    uid : int
        A unique number for each molecule. Generated by the initialzer
        defined in this class.

    """

    def __init__(self, file):
        self.file = file
        self.uid = uuid4()
        self.failed = False

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

        atom = self.mol.GetAtomWithIdx(atom_id)
        atomic_num = atom.GetAtomicNum()
        return periodic_table[atomic_num]

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

        center = np.array([0,0,0])
        total_mass = 0
        for atom_id, coord in self.all_atom_coords():
            mass = self.mol.GetAtomWithIdx(atom_id).GetMass()
            total_mass += mass
            center = np.add(center, np.multiply(mass, coord))

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

    def graph(self):
        """
        Returns a mathematical graph representing the molecule.

        Returns
        -------
        networkx.Graph
            A graph where the nodes are the ids of the atoms in the
            rdkit molecule and the edges are the bonds.

        """

        # Create a graph instance and add the atom ids as nodes. Use the
        # the atom ids from each end of a bond to define edges. Do this
        # for all bonds to account for all edges.

        graph = nx.Graph()

        for atom in self.mol.GetAtoms():
            graph.add_node(atom.GetIdx())

        for bond in self.mol.GetBonds():
            graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        return graph

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

        return np.matrix(pos_array.reshape(-1,3).T)

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
        self.set_position([0,0,0])
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())
        # Apply the rotation.
        self.set_position_from_matrix(new_pos_mat)
        # Return the centroid of the molecule to the origin position.
        self.set_position(og_position)

    def set_orientation(self, start, end):
        """
        Rotates the molecule by a rotation from `start` to `end`.

        Note: The difference between this method and
        `_set_orientation2()` is about which point the rotation
        occurs: centroid of entire molecule versus centroid of bonder
        atoms, respectively.

        Given two direction vectors, `start` and `end`, this method
        applies the rotation required transform `start` to `end` on
        the molecule. The rotation occurs about the centroid of the
        molecule.

        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you can
        associate a gemotric feature of the molecule with a vector, then
        the molecule can be roatated so that this vector is aligned with
        `end`. The defined vector can be virtually anything. This means
        that any geomteric feature of the molecule can be easily aligned
        with any arbitrary axis.

        Parameters
        ----------
        start : numpy.array
            A vector which is to be rotated so that it transforms to the
            `end` vector.

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
        self.set_position([0,0,0])

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield the
        # new atomic positions.
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

        # Replace the old rkdit conformer with one where the centroid is
        # at `position`.
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
            column sets the coordinate of the atom with id 1, and so on.

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
        rdkit.Chem.rdchem.Mol
            A copy of the molecule where the coordinates have been
            shifted by `shift`.

        """

        # The function does not modify the existing conformer, as a
        # result a new instance is created and used for modification.
        conformer = chem.Conformer(self.mol.GetConformer())

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
            # coordinates of an atom in its 'x', 'y' and 'z' attributes.
            atom_position = np.array(conformer.GetAtomPosition(atom_id))

            # Inducing the shift.
            new_atom_position = atom_position + shift

            # Creating a new geometry instance.
            new_coords = rdkit_geo.Point3D(*new_atom_position)

            # Changes the position of the atom in the conformer to the
            # values stored in the new geometry instance.
            conformer.SetAtomPosition(atom_id, new_coords)

        # Create a new copy of the rdkit molecule instance representing
        # the molecule - the original instance is not to be modified.
        new_mol = chem.Mol(self.mol)

        # The new rdkit molecule was copied from the one held in the
        # `mol` attribute, as result it has a copy of its conformer. To
        # prevent the rdkit molecule from holding multiple conformers
        # the `RemoveAllConformers` method is run first. The shifted
        # conformer is then given to the rdkit molecule, which is
        # returned.
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(conformer)
        return new_mol

    def update_from_mae(self, mae_path=None):
        """
        Updates data in attributes to match what is held a .mae file.

        The molecule can be updated from a random .mae file held
        anywhere. Alterntatively, if no path is specified it is assumed
        the .mae file has the same path and name as `file` only with
        the extension replaced with .mae.

        Parameters
        ----------
        mae_path : str (default = None)
            The full path of the .mae file from which the attributes
            should be updated. If ``None`` the .mae file is assumed to
            have the same path and name as `file` only with the
            extension replaced with .mae.

        Modifies
        --------
        self.mol : rdkit.Chem.rdchem.Mol
            The rdkit molecule held in this attribute is changed so that
            it matches the moleclue held in the .mae file.

        self.file's content
            The content in this file is changed to match the content in
            the .mae file.

        Returns
        -------
        None : NoneType

        """

        if mae_path is None:
            mae_path = os.path.splitext(self.file)[0] + '.mae'

        self.mol = mol_from_mae_file(mae_path)
        self.write()

    def write(self, path=None):
        """
        Writes a molecular structure file of the macromolecule.

        The molecule is written to the location in the `file`
        attribute. The function uses the structure of the rdkit molecule
        held in `mol` and as the basis for what is written to the file.

        This bypasses the need to use rdkit's writing functions, which
        have issues with macromolecules due to poor ring finding and
        sanitization issues.

        Parameters
        ----------
        path : str (default = None)
            If the file is to be written to a directory other than
            the one in `file`, it should be written here.

        Modifies
        --------
        file's content
            The conetent in this file will be replaced with the current
            structure of the molecule.

        Returns
        -------
        None : NoneType

        """

        write_funcs = {'.mol' : self._write_mdl_mol_file,
                       '.pdb' : self._write_pdb_file}

        if path is None:
            path = self.file

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
        chem.MolToPDBFile(self.mol, self.file)

        # Edit the file because rkdit does poor atom labelling.
        new_content = ''
        with open(self.file, 'r') as pdb:
            for line in pdb:
                if 'HETATM' in line:
                    words = line.split()
                    lbl_word = words[2]
                    rpl_word = words[-1]
                    rpl_word += " "*(len(lbl_word)-len(rpl_word))
                    line = line.replace(lbl_word, rpl_word)

                new_content += line

        with open(self.file, 'w') as pdb:
            pdb.write(new_content)

@total_ordering
class StructUnit(Molecule, metaclass=Cached):
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
        `init_funcs` dictionary. As long as a file of one of these types
        is provided, MMEA will automatically use the correct
        initialization function.

    mol : rdkit.Chem.rdchem.Mol
        The rdkit instance of the molecule held in `file`.

    func_grp : FGInfo
        The ``FGInfo`` instance holding information about the functional
        group which will react when the building block assembles to form
        macromolecules.

    bonder_ids : list of ints
        A list holding the atom ids of the atoms which form bonds during
        macromolecular assembly.

    energy : Energy
        This attribute handles information about the energy of the
        instance. See the documentation of ``Energy`` to see what is
        available.

    optimized : bool (default = False)
        A flag to monitor whether an optimization has been performed on
        the molecule.

    failed : bool
        True if an optimization function failed to run on the molecule.

    uid : int
        A unique id number.

    """

    init_funcs = {'.mol' : partial(chem.MolFromMolFile,
                                   sanitize=False, removeHs=False),

                  '.sdf' : partial(chem.MolFromMolFile,
                                   sanitize=False, removeHs=False),

                  '.mol2' : partial(chem.MolFromMol2File,
                                 sanitize=False, removeHs=False),
                  '.mae' : mol_from_mae_file,

                  '.pdb' : partial(chem.MolFromPDBFile,
                                 sanitize=False, removeHs=False)}

    def __init__(self, file, functional_group=None):
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

        """

        super().__init__(file)
        _, ext = os.path.splitext(file)

        if ext not in self.init_funcs:
            raise TypeError(
            'Unable to initialize from "{}" files.'.format(ext))

        self.mol = self.init_funcs[ext](file)
        self.bonder_ids = []
        self.energy = Energy(self)
        self.optimized = False

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

        return np.matrix(pos_array.reshape(-1,3).T)

    def centroid_centroid_dir_vector(self):
        """
        Returns the direction vector between the 2 molecular centroids.

        The first molecule centroid is the centroid of the entire
        molecule. The second molecular centroid is the centroid of the
        bonder atoms.

        Returns
        -------
        numpy.array
            The normalized direction vector running from the centroid of
            the bonder atoms to the molecular centroid.

        """

        return normalize_vector(self.centroid() -
                                self.bonder_centroid())

    def json(self):
        """
        Returns a JSON representation of the object.

        The representation has the following form:

            {
                'id' : 5125
                'class' : 'StructUnit',
                'mol_block' : 'A string holding the V3000 mol
                               block of the molecule.'
            }

        Returns
        -------
        str
            A JSON string which represents the molecule.

        """

        return json.dumps({

        'id' : self.uid,
        'class' : self.__class__.__name__,
        'mol_block' : self.mdl_mol_block()

        })

    def functional_group_atoms(self):
        """
        Returns a container of atom ids of atoms in functional groups.

        Returns
        -------
        tuple of tuples of ints
            The form of the returned tuple is:
            ((1,2,3), (4,5,6), (7,8,9)). This means that all atoms with
            ids 1 to 9 are in a functional group and that the atoms 1, 2
            and 3 all form one functional group together. So do 4, 5 and
            5 and so on.

        """

        # Generate a ``rdkit.Chem.rdchem.Mol`` instance which represents
        # the functional group of the molecule.
        func_grp_mol = chem.MolFromSmarts(self.func_grp.fg_smarts)

        # Do a substructure search on the the molecule in `mol` to find
        # which atoms match the functional group. Return the atom ids of
        # those atoms.
        return self.mol.GetSubstructMatches(func_grp_mol)

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
        # Change the position of the centroid of the bonder atoms to the
        # origin so that the rotation occurs about this point.
        self.set_bonder_centroid([0,0,0])
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
            A numpy array holding the desired the position. It holds the
            x, y and z coordinates, respectively.

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

        The great thing about this method is that you as long as you can
        associate a gemotric feature of the molecule with a vector, then
        the molecule can be roatated so that this vector is aligned with
        `end`. The defined vector can be virtually anything. This means
        that any geomteric feature of the molecule can be easily aligned
        with any arbitrary axis.

        Parameters
        ----------
        start : numpy.array
            A vector which is to be rotated so that it transforms to the
            `end` vector.

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
        self.set_bonder_centroid(np.array([0,0,0]))

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield the
        # new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix())

        # Set the positions in the rdkit molecule.
        self.set_position_from_matrix(new_pos_mat)
        self.set_bonder_centroid(og_center)

        return self.mol

    def similar_molecules(self, database):
        """
        Returns molecules from `database` ordered by similarity.

        This method uses the Morgan fingerprints of radius 4 to evaluate
        how similar the molecules in `database` are to `self`.

        Parameters
        ----------
        database : str
            The full path of the database from which the molecules are
            checked for similarity.

        Returns
        -------
        list of tuples of form (float, str)
            The float is the similarity of a given molecule in the
            database to `self` while the str is the full path of the
            .mol file of that molecule.

        """

        # First get the fingerprint of `self`.
        chem.GetSSSR(self.mol)
        self.mol.UpdatePropertyCache(strict=False)
        fp = ac.GetMorganFingerprint(self.mol, 4)

        # For every structure file in the database create a rdkit
        # molecule. Place these in a list.
        mols = []
        for file in os.listdir(database):
            path = os.path.join(database, file)
            # Ignore files which are not structure files and the
            # structure file of the molecule itself.
            _, ext = os.path.splitext(path)
            if ext not in self.init_funcs or path == self.file:
                continue

            mol = self.init_funcs[ext](path)
            chem.GetSSSR(mol)
            mol.UpdatePropertyCache(strict=False)
            mol_fp = ac.GetMorganFingerprint(mol, 4)
            similarity = DataStructs.DiceSimilarity(fp, mol_fp)
            mols.append((similarity, path))

        return sorted(mols, reverse=True)

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

        bonder_ids : list of ints
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
        bond_mol = chem.MolFromSmarts(self.func_grp.target_smarts)
        bond_atoms = self.mol.GetSubstructMatches(bond_mol)
        for atom_id in flatten(bond_atoms):
            atom = self.mol.GetAtomWithIdx(atom_id)
            atom.SetProp('bonder', '1')
            self.bonder_ids.append(atom_id)

        # Give all atoms which form bonds during reactions the tag
        # 'del' and set its value to '1'.
        del_mol = chem.MolFromSmarts(self.func_grp.del_smarts)
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
            The atoms in this rdkit molecule have the properties 'fg',
            'bonder' and 'del' removed.

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


    def __eq__(self, other):
        return self.file == other.file

    def __lt__(self, other):
        return self.file < other.file

    def __hash__(self):
        return id(self)

    def __str__(self):
        return self.file

    def __repr__(self):
        repr_ =  "{0!r}".format(type(self))
        repr_ = repr_.replace(">",
        ", file={0.file!r}>".format(self))
        repr_ = repr_.replace("class ", "class=")
        return repr_

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

    def minimize_theta(self, vector, axis, step=0.17):
        """
        Rotates molecule about `axis` to minimze theta with `vector`.

        The molecule is iteratively rotated about `axis` so that the
        vector between its bonder atoms is as close as possible to
        `vector`.

        Parameters
        ----------
        vector : numpy.array
            The vector to which the distance should be minimized.

        axis : numpy.array
            The direction vector along which the rotation happens.

        step : float
            The size of the iterative step in radians.

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

        vector = normalize_vector(vector)
        axis = normalize_vector(axis)

        theta = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)

        # First determine the direction in which iteration should occur.
        self.rotate2(step, axis)
        theta2 = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)
        if theta2 > theta:
            axis = np.multiply(axis, -1)

        prev_theta = theta2
        while True:
            self.rotate2(step, axis)
            theta = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)

            if theta > prev_theta:
                axis = np.multiply(axis, -1)
                self.rotate2(step, axis)
                break

            prev_theta = theta


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
        coordinate of some point on the plane. For example, the position
        of one of the bonder atoms.

        Returns
        -------
        numpy.array
            This array has the form [a, b, c, d] and represents the
            scalar equation of the plane formed by the bonder atoms.

        References
        ----------
        http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx

        """

        bonder_coord = self.atom_coords(self.bonder_ids[0])
        d = np.multiply(np.sum(np.multiply(self.bonder_plane_normal(),
                                           bonder_coord)), -1)
        return np.append(self.bonder_plane_normal(), d)

    def bonder_plane_normal(self):
        """
        Returns the normal vector to the plane formed by bonder atoms.

        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the molecule.

        Returns
        -------
        numpy.array
            A unit vector which describes the normal to the plane of the
            bonder atoms.

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
            normal_v = np.multiply(normal_v, -1)

        return normal_v

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
class MacroMolecule(Molecule, metaclass=CachedMacroMol):
    """
    A class for assembled macromolecules.

    The goal of this class is to represent an individual used by the GA.
    As such, it holds attributes that are to be expected for this
    purpose. Mainly, it has a fitness value stored in its `fitness`
    attribute and a genetic code - as defined by its `building_blocks`
    and  `topology` attributes. If a change is made to either of these
    attributes, they should describe a different macromolecule. On the
    other hand, the same attributes should always describe the same
    macromolecule.

    Because of this, as well as the computational cost associated with
    macromolecule initialization, instances of this class are cached.
    This means that providing the same arguments to the initializer will
    not build a different instance with the same attribute values. It
    will yield the original instance, retrieved from memory.

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
    ``StructUnit`` instance representing a specific building block needs
    to be held. How many of such building blocks are need to assemble
    the molecule is the handled by the ``Topology`` class.

    This class is not intended to be used directly but should be
    inherited by subclasses representing specific macromolecules. The
    ``Cage`` and ``Polymer`` classes are examples of this. Any
    information or methods that apply generally to all macromolecules
    should be defined within this class while specific non-general data
    should be included in derived classes.

    This class also supports comparison operations, these act on the
    fitness value assiciated with a macromolecule. Comparison operations
    not explicitly defined are included via the ``total_ordering``
    decorator. For other operations and methods supported by this class
    examine the rest of the class definition.

    Finally, a word of caution. The equality operator ``==`` compares
    fitness values. This means two macromolecules, made from different
    building blocks, can compare equal if they happen to have the same
    fitness. The operator is not to be used to check if one
    macromolecule is the same structurally as another. To do this check
    use the ``same()`` method. This method may be overwritten in derived
    classes, as necessary. In addition the ``is`` operator is
    implemented as is default in Python. It compares whether two objects
    are in the same location in memory. Because the ``MacroMolecule``
    class is cached the ``is`` operator could in principle be used
    instead of the `same` method (including in derived classes).
    However, this is not intended use and is not guuranteed to work in
    future implementations. If caching stops being implemented such code
    would break.

    Attributes
    ----------
    building_blocks : list of ``StructUnit`` instances
        This attribute holds ``StructUnit`` instances which represent
        the monomers forming the macromolecule. Only one ``StructUnit``
        instance is needed per monomer, even if multiples of a monomer
        join up to form the macromolecule

    topology : A child class of ``Topology``
        This instance represents the topology of the macromolecule. Any
        information to do with how individual building blocks of the
        macromolecule are organized and joined up in space is held by
        this attribute. For more details about what information and
        functions this entails see the docstring of the ``Topology``
        class and its derived classes.

    topology_args : list (default = [])
        This attribue holds the initializer arguments for the topology
        instance. This is stored so that exceptions can print all
        values required to make an identical copy of a ``MacroMolecule``
        instance.

    file : str
        The full path of the molecule structure file holding the
        macromolecule.

    mol : rdkit.Chem.rdchem.Mol
        An rdkit instance representing the macromolecule.

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

    failed : bool (default = False)
        ``True`` if the fitness or optimization function failed to
        evaluate the molecule.

    progress_params : list (default = None)
        Holds the fitness parameters which the GA should track to make
        progress plots. If the default ``None`` is used, the fitness
        value will be used.

    key : str
        The key used for caching the molecule. Necessary for
        `update_cache` to work. This attribute is assigned by the
        `__call__()` method of the ``CachedMacroMol`` class.

    uid : int
        A unique id number.

    """

    def __init__(self, building_blocks, topology, file,
                 topology_args=None):
        """
        Initialize a ``MacroMolecule`` instance.

        When an exception occurs during initialization, all parameters
        which were provided to the initializer are saved to a file
        ``failures.txt`` which is in the ``output`` folder.

        Parameters
        ---------
        building_blocks : list of ``StructUnit`` instances
            A list of ``StructUnit`` instances which represent the
            monomers forming the macromolecule.

        topology : A child class of ``Topology``
            The class which defines the topology of the macromolecule.
            Such classes are defined in the topology module. The class
            will be a child class which inherits ``Topology``.

        file : str
            The full path of the structure file where the macromolecule
            will be stored.

        topology_args : dict (default = None)
            Any additional arguments needed to initialize the topology
            class supplied in the `topology` argument.

        """

        try:
            self._std_init(building_blocks,
                           topology, file, topology_args)

        except Exception as ex:
            self.building_blocks = building_blocks
            self.topology = topology
            self.mol = chem.Mol()
            self.file = file
            self.topology_args = topology_args
            self.fail()
            MolError(ex, self, 'During initialization.')

    def _std_init(self, building_blocks, topology, file,
                 topology_args):

        super().__init__(file)

        if topology_args is None:
            topology_args = {}

        self.optimized = False
        self.fitness = None
        self.unscaled_fitness = {}
        self.progress_params = None
        self.building_blocks = tuple(building_blocks)
        self.topology_args = topology_args
        self.energy = Energy(self)

        self.topology = topology(self, **topology_args)
        # Ask the ``Topology`` instance to assemble/build the cage. This
        # creates the cage's structure file with all the building blocks
        # and linkers joined up.
        self.topology.build()

        # Write the structure file of the assembled molecule.
        self.write()

    def json(self):
        """
        Returns a JSON representation of the object.

        The representation has the following form:

            {
                'id' : 5125
                'class' : 'Polymer',
                'mol_block' : 'A string holding the V3000 mol
                               block of the molecule.'
                'building_blocks' : [bb1.json(), bb2.json()]
                'topology' : 'Copolymer'
                'topology_args' : {'arg1' : arg1_val,
                                   'arg2' : arg2_val},
                'unscaled_fitness' : {'fitness_func1' : fitness1,
                                      'fitness_func2' : fitness2}

            }

        Returns
        -------
        str
            A JSON string which represents the molecule.

        """

        return json.dumps({

        'id' : self.uid,
        'class' : self.__class__.__name__,
        'mol_block' : self.mdl_mol_block()

        })

    def same(self, other):
        """
        Check if the `other` instance has the same molecular structure.

        Parameterss.getcwd()
        ----------
        other : MacroMolecule
            The ``MacroMolecule`` instance you are checking has the same
            structure.

        Returns
        -------
        bool
            Returns ``True`` if the building blocks and topology of the
            macromolecules are all the same.

        """

        # Compare the building blocks and topology making up the
        # macromolecule. If these are the same then the molecules have
        # the same structure.
        return (self.building_blocks == other.building_blocks and
                type(self.topology) == type(other.topology) and
                self.topology_args == other.topology_args)

    def update_cache(self):
        """
        Updates the caching dictionary so that it contains `self`.

        When an instance of ``MacroMolecule`` is first created it is
        cached. Using multiprocessing to perform optimizations returns
        modified copies of the cached molecules. In order to ensure that
        the cache points to the modified copies not to the originally
        initialized molecule, this method must be run.

        Returns
        -------
        None : NoneType

        """

        # The caching is done by the class. Access it and use its
        # ``_update_cache()`` method.
        cls = type(self)
        cls._update_cache(self)

    def __eq__(self, other):
        return self.fitness == other.fitness

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __str__(self):
        return str({key: value for key, value in
                                    self.__dict__.items() if
                                    key in {'file',
                                            'topology',
                                            'fitness',
                                            'optimized',
                                            'building_blocks'}}) + "\n"

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
    def testing_init(cls, bb_str, lk_str, topology_str):
        key = (bb_str, lk_str, topology_str)
        if key in MacroMolecule._cache.keys():
            return MacroMolecule._cache[key]
        else:
            cage = cls.__new__(cls)
            cage.building_blocks = (bb_str, lk_str)
            cage.topology = topology_str
            cage.topology_args = {'a' : 1}
            cage.fitness = 3.14
            MacroMolecule._cache[key] = cage
            return cage


class Cage(MacroMolecule):
    """
    Used to represent molecular cages.

    """

    @classmethod
    def init_fixed_bb(cls, bb_file, lk_db, topologies, file):
        bb = StructUnit3(bb_file)

        while True:
            try:
                lk_file = np.random.choice(os.listdir(lk_db))
                lk_file = os.path.join(lk_db, lk_file)
                lk = StructUnit2(lk_file)
                break
            except TypeError:
                continue

        topology = np.random.choice(topologies)

        return cls((bb, lk), topology, file)

    @classmethod
    def init_random(cls, bb_db, lk_db, topologies, file):
        """
        Makes ``Cage`` from random building blocks and topology.

        Parameters
        ----------
        bb_db : str

        lk_db : str

        topologies : list of ``Topology`` child classes.

        file : str

        """

        while True:
            try:
                bb_file = np.random.choice(os.listdir(bb_db))
                bb_file = os.path.join(bb_db, bb_file)
                bb = StructUnit3(bb_file)
                break

            except TypeError:
                continue

        while True:
            try:
                lk_file = np.random.choice(os.listdir(lk_db))
                lk_file = os.path.join(lk_db, lk_file)
                lk = StructUnit(lk_file)

                if len(lk.bonder_ids) >= 3:
                    lk = StructUnit3(lk_file)
                else:
                    lk = StructUnit2(lk_file)

                break

            except TypeError:
                continue

        topology = np.random.choice(topologies)
        return cls((bb, lk), topology, file)

class Polymer(MacroMolecule):
    """
    Used to represent polymers.

    """
    pass
