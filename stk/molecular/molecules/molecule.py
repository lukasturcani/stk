import json
import os
import numpy as np
import networkx as nx
import rdkit.Geometry.rdGeometry as rdkit_geo
import rdkit.Chem.AllChem as rdkit
from rdkit.Chem import rdMolTransforms
from scipy.spatial.distance import euclidean
from sklearn.metrics.pairwise import euclidean_distances
from collections import defaultdict

from ...utilities import (normalize_vector,
                          rotation_matrix,
                          mol_from_mae_file,
                          rotation_matrix_arbitrary_axis,
                          atom_vdw_radii,
                          remake)


class MoleculeSubclassError(Exception):
    ...


class Molecule:
    """
    The most basic class representing molecules.

    This class defines the operations which any class
    describing molecules should inherit or may find useful. Examples of
    such are :class:`StructUnit` and :class:`MacroMolecule`. This class
    should not be used directly.

    Attributes
    ----------
    mol : :class:`rdkit.Mol`
        An ``rdkit`` molecule instance representing the molecule.

    inchi : :class:`str`
        The InChI of the molecule.

    atom_props : :class:`dict`
        Maps atom id to a :class:`dict` holding the properties of
        that atom. For example

        .. code-block:: python

            atom_props = {0: {'prop1': 0,
                              'prop2': 'value1',
                              'prop3': 10.},

                          5: {'prop1': 'value2',
                              'prop5': 2.0}}

    name : :class:`str`
        A name which can be optionally given to the molecule for easy
        identification.

    note : :class:`str`
        A note or comment about the molecule. Purely optional but can
        be useful for labelling and debugging.

    """

    subclasses = {}

    def __init__(self, name="", note=""):
        self.name = name
        self.note = note
        self.atom_props = defaultdict(dict)

    def __init_subclass__(cls, **kwargs):
        if cls.__name__ in cls.subclasses:
            msg = 'Subclass with this name already exists.'
            raise MoleculeSubclassError(msg)
        cls.subclasses[cls.__name__] = cls
        super().__init_subclass__(**kwargs)

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

    def atom_centroid(self, atom_ids, conformer=-1):
        """
        Return the centroid of a group of atoms.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`
            The ids of atoms which which are used to calculate the
            centroid.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of atoms specified by `atom_ids`.

        """

        coords = (self.atom_coords(a, conformer) for a in atom_ids)
        return sum(coords) / len(atom_ids)

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
        :class:`numpy.ndarray`
            An array holding the x, y and z coordinates of the atom.

        """

        conf = self.mol.GetConformer(conformer)
        return np.array(conf.GetAtomPosition(atom_id))

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

    def center_of_mass(self, conformer=-1):
        """
        Returns the centre of mass of the molecule.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
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
        :class:`numpy.ndarray`
            A numpy array holding the position of the centroid.

        """

        coords = (
            coord for _, coord in self.all_atom_coords(conformer)
        )
        return sum(coords) / self.mol.GetNumAtoms()

    def core(self):
        """
        Return the molecule with no H or functional group atoms.

        Returns
        -------
        :class:`rdkit.Mol`
            The "core" of the molecule.

        """

        emol = rdkit.EditableMol(self.mol)
        for atom in reversed(self.mol.GetAtoms()):
            atomid = atom.GetIdx()
            if not self.is_core_atom(atomid):
                emol.RemoveAtom(atomid)
        return emol.GetMol()

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
                # Calculate the dihedral angle.
                dihedral_value = rdMolTransforms.GetDihedralDeg(
                                    self.mol.GetConformer(conformer),
                                    atoms_group[0],
                                    atoms_group[1],
                                    atoms_group[2],
                                    atoms_group[3])
                # Check that the dihedral is calculated in the right
                # direction.
                if abs(dihedral_value) > 90:
                    dihedral_value = abs(dihedral_value)
                else:
                    dihedral_value = 180 - abs(dihedral_value)

                dihedral_info.append(dihedral_value)

            # Calculate the average dihedral value.
            avg_dihedral = np.mean([abs(x) for x in dihedral_info])
            # Calculate the relative diff with the target dihedral
            # value.
            diff = (abs(target - avg_dihedral) / target) * 100
        else:
            # If the molecule does not contain the bond, give 1%
            # strain.
            diff = 1

        return diff

    def dump(self, path, include_attrs=None):
        """
        Writes a JSON :class:`dict` of the molecule to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file to which the JSON dict should be
            written.

        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecule to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            json.dump(self.json(include_attrs), f, indent=4)

    @classmethod
    def from_dict(self, json_dict, load_names=True):
        """
        Creates a :class:`Molecule` from a JSON :class:`dict`.

        The :class:`Molecule` returned has the class specified in
        `json_dict`, not :class:`Molecule`.

        Parameters
        ----------
        json_dict : :class:`dict`
            A :class:`dict` holding the JSON representation of a
            molecule.

        load_names : :class:`bool`, optional
            If ``True`` then the ``name`` key stored in `json_dict`
            is loaded.

        Returns
        -------
        :class:`Molecule`
            The molecule represented by `json_dict`.

        """

        # Get the class of the object.
        c = self.subclasses[json_dict['class']]
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

    def is_core_atom(self, atom_id):
        """
        Returns ``True`` if atom is not H or part of a fg.

        Parameters
        ----------
        atom_id : :class:`int`
            The id of the atom being queried.

        Returns
        -------
        :class:`bool`
            Indicates whether the atom with `atom_id` is part of the
            core.

        """

        atom = self.mol.GetAtomWithIdx(atom_id)
        if atom.GetAtomicNum() == 1:
            return False
        return all(
            atom_id not in fg.atom_ids for fg in self.func_groups
        )

    @classmethod
    def load(cls, path, load_names=True):
        """
        Creates a :class:`Molecule` from a JSON file.

        The returned :class:`Molecule` has the class specified in the
        JSON file, not :class:`Molecule`.

        Parameters
        ----------
        path : :class:`str`
            The full path holding a JSON representation to a molecule.

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

        return cls.from_dict(json_dict, load_names)

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

        coords = self.mol.GetConformer(conformer).GetPositions()
        dist = euclidean_distances(coords, coords)
        vdw = np.array([[atom_vdw_radii[self.atom_symbol(i)] for
                        i in range(self.mol.GetNumAtoms())]])
        dist = dist + vdw + vdw.T
        maxid1, maxid2 = np.unravel_index(dist.argmax(), dist.shape)
        return dist[maxid1, maxid2], int(maxid1), int(maxid2)

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

        # Kekulize the mol, which means that each aromatic bond is
        # converted to a single or double. This is necessary because
        # .mol V3000 only supports integer bonds. However, this fails
        # sometimes on big molecules.
        try:
            rdkit.Kekulize(self.mol)
        except ValueError:
            pass

        n_atoms = self.mol.GetNumAtoms()
        n_bonds = self.mol.GetNumBonds()

        dtype = np.dtype(object)

        atom_ids = np.array([np.arange(1, n_atoms+1)],
                            dtype=dtype).T

        atom_symbols = np.array([[self.atom_symbol(i)]
                                 for i in range(n_atoms)],
                                dtype=dtype)

        pos_mat = self.mol.GetConformer(conformer).GetPositions()

        charges = np.array([[f' CHG={a.GetFormalCharge()}' if
                             a.GetFormalCharge() else '']
                            for a in self.mol.GetAtoms()],
                           dtype=dtype)

        atom_data = np.concatenate(
                        [atom_ids, atom_symbols, pos_mat, charges],
                        axis=1).reshape((-1, ))
        atom_block = "M  V30 {} {} {:.4f} {:.4f} {:.4f} 0{}\n"*n_atoms
        atom_block = atom_block.format(*atom_data)

        bond_data = [prop for bond in self.mol.GetBonds() for prop in
                     (bond.GetIdx(),
                      int(bond.GetBondTypeAsDouble()),
                      bond.GetBeginAtomIdx()+1,
                      bond.GetEndAtomIdx()+1)]
        bond_block = "M  V30 {} {} {} {}\n"*n_bonds
        bond_block = bond_block.format(*bond_data)

        return ("\n"
                "     RDKit          3D\n"
                "\n"
                "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                "M  V30 BEGIN CTAB\n"
                f"M  V30 COUNTS {n_atoms} {n_bonds} 0 0 0\n"
                "M  V30 BEGIN ATOM\n"
                f"{atom_block}"
                "M  V30 END ATOM\n"
                "M  V30 BEGIN BOND\n"
                f"{bond_block}"
                "M  V30 END BOND\n"
                "M  V30 END CTAB\n"
                "M  END\n"
                "\n"
                "$$$$\n")

    def position_matrix(self, conformer=-1):
        """
        Returns a matrix holding the atomic positions of a conformer.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            The array has the shape ``[3, n]``. Each column holds the
            x, y and z coordinates of a bonder centroid. The index of
            the column corresponds to the atom id.

        """

        return self.mol.GetConformer(conformer).GetPositions().T

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

    def save_rdkit_atom_props(self, prop_names):
        """
        Updates :attr:`~.Molecule.atom_props` with rdkit atom tags.

        Parameters
        ----------
        prop_names : :class:`set` of :class:`str`
            The names of atom properties which should be saved.

        Returns
        -------
        None : :class:`NoneType`

        """

        for atom in self.mol.GetAtoms():
            atom_id = atom.GetIdx()
            props = atom.GetPropsAsDict()
            valid_props = props.keys() & prop_names
            for prop_name in valid_props:
                self.atom_props[atom_id][prop_name] = props[prop_name]

    def rotate(self, theta, axis, conformer=-1):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the molecular centroid.

        Parameters
        ----------
        theta : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.ndarray`
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
        pos_mat = self.mol.GetConformer(conformer).GetPositions()
        new_pos_mat = rot_mat @ pos_mat.T
        # Apply the rotation.
        self.set_position_from_matrix(new_pos_mat, conformer)
        # Return the centroid of the molecule to the origin position.
        self.set_position(og_position, conformer)

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
        start : :class:`numpy.ndarray`
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : :class:`numpy.ndarray`
            This array holds the vector, onto which `start` is rotated.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
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
        pos_mat = self.mol.GetConformer(conformer).GetPositions().T
        new_pos_mat = np.dot(rot_mat, pos_mat)

        # Set the positions of the molecule.
        self.set_position_from_matrix(new_pos_mat, conformer)
        self.set_position(og_center, conformer)

        return self.mol

    def set_position(self, position, conformer=-1):
        """
        Sets the centroid of the molecule to `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            This array holds the position on which the centroid of the
            molecule should be placed.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`rdkit.Mol`
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
        pos_mat : :class:`numpy.ndarray`
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
        shift : :class:`numpy.ndarray`
            A numpy array holding the value of the shift along each
            axis.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
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

    def update_cache(self):
        """
        Update attributes of cached molecule.

        Using ``multiprocessing`` returns modified copies of molecules.
        In order to ensure that the cached molecules have
        their attributes updated to the values of the copies, this
        method must be run on the copies.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.key in self.__class__.cache:
            self.__class__.cache[self.key].__dict__ = dict(vars(self))

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

        mol = remake(rdkit.MolFromMolFile(molFileName=path,
                                          sanitize=False,
                                          removeHs=False))
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
