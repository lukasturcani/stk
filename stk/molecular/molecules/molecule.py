import json
import os
import numpy as np
import rdkit.Geometry.rdGeometry as rdkit_geo
import rdkit.Chem.AllChem as rdkit
from scipy.spatial.distance import euclidean

from ...utilities import (
    normalize_vector,
    rotation_matrix,
    mol_from_mae_file,
    rotation_matrix_arbitrary_axis,
    remake,
    periodic_table
)


class MoleculeSubclassError(Exception):
    ...


class Atom:
    """
    Base class for atoms.

    A subclass is made for each element. The name of each elemental
    class is the periodic table symbol.

    Attributes
    ----------
    atomic_number : :class:`int`
        A class attribute. Specifies the atomic number.

    """

    ...


class Bond:
    """
    Represents an atomic bond.

    Attributes
    ----------
    order : :class:`float`
        The bond order.

    atom1 : :class:`Atom`
        The first atom in the bond.

    atom2 : :class:`Atom`
        The second atom in the bond.

    """

    def __init__(self, order, atom1, atom2):
        self.order = order
        self.atom1 = atom1
        self.atom2 = atom2


class Molecule:
    """
    The most basic class representing molecules.

    This class defines the operations which any class
    describing molecules should inherit or may find useful. Examples of
    such are :class:`.BuildingBlock` and :class:`.ConstructedMolecule`.
    This class should not be used directly.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms which compose the molecule.

    bonds : :class:`tuple` of :class:`.Bond`
        The bonds of the molecule.

    inchi : :class:`str`
        The InChI of the molecule.

    name : :class:`str`
        A name which can be optionally given to the molecule for easy
        identification.

    key : :class:`object`
        A hashable :class:`object`. This attribute will be the same
        for molecules of the same class, which have the same structure.
        A private method :meth:`_generate_key` must be defined for
        each subclass of :class:`.Molecule` and it will be used to
        generate the :attr:`key`.

    cache : :class:`dict`
        This is a class attribute. Which maps :attr:`key` to the
        :class:`.Molecule` instance with that :attr:`key`.

    _mol : :class:`rdkit.Mol`
        A :mod:`rdkit` molecule instance representing the molecule.

    Methods
    -------

    """

    subclasses = {}

    def __init__(self, name):
        self.name = name

    @classmethod
    def init_from_dict(self, json_dict, load_names=True):
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

    def __init_subclass__(cls, **kwargs):
        if cls.__name__ in cls.subclasses:
            msg = 'Subclass with this name already exists.'
            raise MoleculeSubclassError(msg)
        cls.subclasses[cls.__name__] = cls
        cls.cache = {}
        super().__init_subclass__(**kwargs)

    def apply_displacement(self, displacement, conformer=-1):
        """
        Shift the molecule by `displacement`.

        Parameters
        ----------
        displacement : :class:`numpy.ndarray`
            A displacement vector applied to the molecule.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        position_matrix = self.get_position_matrix(False, conformer)
        self.set_position_matrix(
            position_matrix=(position_matrix+displacement).T,
            conformer=conformer
        )
        return self

    def apply_rotation(self, theta, axis, origin, conformer=-1):
        """
        Rotate the molecule by `theta` about `axis` on the `origin`.

        Parameters
        ----------
        theta : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.ndarray`
            The axis about which the rotation happens.

        origin : :class:`numpy.ndarray`
            The origin about which the rotation happens.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        # Set the origin of the rotation to the origin.
        self.apply_displacement(-origin, conformer)
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        pos_mat = self.get_position_matrix(conformer=conformer)
        # Apply the rotation.
        self.set_position_matrix(rot_mat @ pos_mat, conformer)
        # Return the centroid of the molecule to the original position.
        self.apply_displacement(origin, conformer)
        return self

    def get_atom_coords(self, atom_ids=None, conformer=-1):
        """
        Yield the coordinates of atoms.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`, optional
            The ids of the atoms whose coordinates are desired.
            If ``None``, then the coordinates of all atoms will be
            yielded.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Yields
        ------
        :class:`numpy.ndarray`
            An array holding the x, y and z coordinates of the
            next atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        conf = self._mol.GetConformer(conformer)
        for atom_id in atom_ids:
            yield np.array(conf.GetAtomPosition(atom_id))

    def get_atom_distance(self, atom1_id, atom2_id, conformer=-1):
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
        :class:`float`
            The distance between the first and second atoms.

        """

        coords = self.get_atom_coords(
            atom_ids=(atom1_id, atom2_id),
            conformer=conformer
        )
        return float(euclidean(*coords))

    def get_center_of_mass(self, atom_ids=None, conformer=-1):
        """
        Return the centre of mass of a group of atoms.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            center of mass. If ``None``, then all atoms will be used.

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

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        center = np.array([0., 0., 0.])
        total_mass = 0.
        coords = self.get_atom_coords(atom_ids, conformer)
        for atom_id, coord in zip(atom_ids, coords):
            mass = self.atoms[atom_id].mass
            total_mass += mass
            center += mass*coord
        return np.divide(center, total_mass)

    def get_centroid(self, atom_ids=None, conformer=-1):
        """
        Return the centroid of a group of atoms.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`, optional
            The ids of atoms which are used to calculate the
            centroid. If ``None``, then all atoms will be used.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of atoms specified by `atom_ids`.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        coords = self.get_atom_coords(atom_ids, conformer)
        return sum(coords) / len(atom_ids)

    def get_direction(self, atom_ids=None, conformer=-1):
        """
        Return a vector of best fit through the molecule's atoms.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            vector. If ``None``, then all atoms will be used.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The vector of best fit.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        pos = self.get_position_matrix(False, conformer)[:, atom_ids]
        return np.linalg.svd(pos - pos.mean(axis=1))[-1][0]

    def get_plane_normal(self, atom_ids=None, conformer=-1):
        """
        Return the normal to the plane of best fit.

        Parameters
        ----------
        atom_ids : :class:`list` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            plane. If ``None``, then all atoms will be used.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`numpy.ndarray`
            Vector orthonormal to the plane of the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        pos = self.get_position_matrix(False, conformer)[:, atom_ids]
        centroid = self.get_centroid(atom_ids, conformer)
        return np.linalg.svd(pos - centroid)[-1][2, :]

    def get_position_matrix(self, atom_columns=True, conformer=-1):
        """
        Returns a matrix holding the atomic positions of a conformer.

        Parameters
        ----------
        atom_columns : :class:`bool`, optional
            The matrix has a shape of ``[3, n]`` if this is ``True``
            and ``[n, 3]`` if it is ``False``.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            The array has the shape ``[3, n]``. Each column holds the
            x, y and z coordinates of a bonder centroid. The index of
            the column corresponds to the atom id.

        """

        if atom_columns:
            self._mol.GetConformer(conformer).GetPositions().T
        else:
            self._mol.GetConformer(conformer).GetPositions()

    def dump(self, path, include_attrs=None):
        """
        Write a JSON :class:`dict` of the molecule to a file.

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
            json.dump(self.to_json(include_attrs), f, indent=4)

    @property
    def inchi(self):
        """
        Returns the InChI of the molecule.

        Returns
        -------
        :class:`str`
            The InChI of the molecule.

        """

        self._update_stereochemistry()
        return rdkit.MolToInchi(self.mol)

    @staticmethod
    def _generate_key(*args, **kwargs):
        """

        """

        raise NotImplementedError()

    @classmethod
    def load(cls, path, load_names=True):
        """
        Create a :class:`Molecule` from a JSON file.

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

        return cls.init_from_dict(json_dict, load_names)

    def to_mdl_mol_block(self, atoms=None, conformer=-1):
        """
        Return a V3000 mol block of the molecule.

        Parameters
        ----------
        atoms : :class:`set` of :class:`int`, optional
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

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

        if atoms is None:
            atoms = range(self.mol.GetNumAtoms())

        n_atoms = len(atoms)
        atom_lines = []
        for atom_id in atoms:
            x, y, z = self.atom_coords(atom_id, conformer)
            symbol = self.atom_symbol(atom_id)
            charge = self.mol.GetAtomWithIdx(atom_id).GetFormalCharge()
            charge = f' CHG={charge}' if charge else ''
            atom_lines.append(
                'M  V30 {} {} {:.4f} {:.4f} {:.4f} 0{}\n'.format(
                    atom_id+1, symbol, x, y, z, charge
                )
            )
        atom_block = ''.join(atom_lines)

        # Convert to set because membership is going to be checked by
        # bonds.
        atoms = set(atoms)
        bond_lines = []
        for bond in self.mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in atoms and a2 in atoms:
                # Keep bond ids if all bonds are getting written.
                if n_atoms == self.mol.GetNumAtoms():
                    bond_id = bond.GetIdx()
                else:
                    bond_id = len(bond_lines)
                bond_type = int(bond.GetBondTypeAsDouble())
                bond_lines.append(
                    f'M  V30 {bond_id+1} {bond_type} {a1+1} {a2+1}\n'
                )

        n_bonds = len(bond_lines)
        bond_block = ''.join(bond_lines)

        return (
            '\n'
            '     RDKit          3D\n'
            '\n'
            '  0  0  0  0  0  0  0  0  0  0999 V3000\n'
            'M  V30 BEGIN CTAB\n'
            f'M  V30 COUNTS {n_atoms} {n_bonds} 0 0 0\n'
            'M  V30 BEGIN ATOM\n'
            f'{atom_block}'
            'M  V30 END ATOM\n'
            'M  V30 BEGIN BOND\n'
            f'{bond_block}'
            'M  V30 END BOND\n'
            'M  V30 END CTAB\n'
            'M  END\n'
            '\n'
            '$$$$\n'
        )

    def is_same_molecule(self, other):
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
        :meth:`.BuildingBlock._set_orientation2` is about which point
        the rotation occurs: centroid of the entire molecule versus
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

    def set_centroid(self, position, conformer=-1):
        """
        Set the centroid of the molecule to `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            This array holds the position on which the centroid of the
            molecule is going to be placed.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`Molecule`
            The molecule is returned.

        """

        centroid = self.get_centroid(conformer=conformer)
        self.apply_displacement(-centroid, conformer)
        return self

    def set_position_matrix(self, position_matrix, conformer=-1):
        """
        Set the molecule's coordinates to those in `position_matrix`.

        Parameters
        ----------
        position_matrix : :class:`numpy.ndarray`
            A position matrix of the molecule. The shape of the matrix
            is ``[3, n]``.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`Molecule`
            The molecule is returned.

        """

        conf = self.mol.GetConformer(conformer)
        for i, coord_mat in enumerate(pos_mat.T):
            coord = rdkit_geo.Point3D(coord_mat.item(0),
                                      coord_mat.item(1),
                                      coord_mat.item(2))
            conf.SetAtomPosition(i, coord)

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

        mol = Molecule()
        mol.mol = mol_from_mae_file(path)
        self.set_position_from_matrix(mol.position_matrix(), conformer)

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

        mol = Molecule()
        mol.mol = remake(
            rdkit.MolFromMolFile(
                molFileName=path,
                sanitize=False,
                removeHs=False
            )
        )
        self.set_position_from_matrix(mol.position_matrix(), conformer)

    def update_from_xyz(self, path, conformer=-1):
        """
        Updates molecular structure to match a ``.xyz`` file.

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

        Raises
        ------
        :class:`RuntimeError`
            If the number of atoms in the file does not match the
            number of atoms in the molecule or if atom elements in the
            file do not agree with the atom elements in the molecule.

        """

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()

        with open(path, 'r') as f:
            atom_count, _, *content = f.readlines()

        # Check the atom count is correct.
        num_atoms = self.mol.GetNumAtoms()
        if int(atom_count) != num_atoms:
            raise RuntimeError(
                f'The number of atoms in the xyz file, {atom_count}, '
                'does not match the number of atoms in the '
                f'molecule, {num_atoms}.'
            )

        # Save all the coords in the file.
        new_coords = []
        for i, line in enumerate(content):
            element, *coords = line.split()
            if element.isnumeric():
                element = periodic_table[int(element)]

            if element != self.atom_symbol(i):
                raise RuntimeError(
                    f'Atom {i} element does not match file.'
                )

            new_coords.append([float(i) for i in coords])

        # Check that the correct number of atom
        # lines was present in the file.
        if i+1 != num_atoms:
            raise RuntimeError(
                f'The number of atom lines in the xyz file, {i+1}, '
                'does not match the number of atoms in the '
                f'molecule, {num_atoms}.'
            )

        # Update the structure.
        new_coords = np.array(new_coords).T
        self.set_position_from_matrix(new_coords, conformer=conformer)

    def update_from_turbomole(self, path, conformer=-1):
        """
        Updates molecular structure from a Turbomole ``.coord`` file.

        Note that coordinates in ``.coord`` files are given in Bohr.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.coord`` file from which the
            structure should be updated.

        conformer : :class:`int`, optional
            The conformer to be updated.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`RuntimeError`
            If the number of atoms in the file does not match the
            number of atoms in the molecule or if atom elements in the
            file do not agree with the atom elements in the molecule.

        """
        bohr_to_ang = 0.5291772105638411

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()

        with open(path, 'r') as f:
            _, *content, __ = f.readlines()

        # Check the atom count is correct.
        num_atoms = self.mol.GetNumAtoms()
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

            if element != self.atom_symbol(i):
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
        new_coords = np.array(new_coords).T
        self.set_position_from_matrix(new_coords, conformer=conformer)

    def _update_stereochemistry(self, conformer=-1):
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

        for atom in self._mol.GetAtoms():
            atom.UpdatePropertyCache()
        rdkit.AssignAtomChiralTagsFromStructure(self._mol, conformer)
        rdkit.AssignStereochemistry(self._mol, True, True, True)

    def write(self, path, atoms=None, conformer=-1):
        """
        Writes a molecular structure file of the molecule.

        This bypasses the need for writing functions in ``rdkit``.
        These have issues with large molecules due to poor ring finding
        and sanitization issues.

        Parameters
        ----------
        path : :class:`str`
            The `path` to which the molecule should be written.

        atoms : :class:`list` of :class:`int`, optional
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        write_funcs = {
            '.mol': self._write_mdl_mol_file,
            '.sdf': self._write_mdl_mol_file,
            '.xyz': self._write_xyz_file,
            '.pdb': self._write_pdb_file
        }

        _, ext = os.path.splitext(path)
        write_func = write_funcs[ext]
        write_func(path, atoms, conformer)

    def _write_mdl_mol_file(self, path, atoms, conformer):
        """
        Writes a V3000 ``.mol`` file of the molecule

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atoms : :class:`list` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            f.write(self.mdl_mol_block(atoms, conformer))

    def _write_xyz_file(self, path, atoms, conformer):
        """
        Writes a ``.xyz`` file of the molecule

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atoms : :class:`list` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()
        if atoms is None:
            atoms = range(self.mol.GetNumAtoms())

        num_atoms = str(len(atoms))

        content = [f'{num_atoms}\n\n']
        for atom_id in atoms:
            x, y, z = self.atom_coords(atom_id, conformer)
            symbol = self.atom_symbol(atom_id)
            content.append(f'{symbol} {x:f} {y:f} {z:f}\n')

        with open(path, "w") as xyz:
            xyz.write(''.join(content))

    def _write_pdb_file(self, path, atoms, conformer):
        """
        Writes a ``.pdb`` file of the molecule

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atoms : :class:`list` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atoms is None:
            atoms = range(self.mol.GetNumAtoms())

        if conformer == -1:
            conformer = self.mol.GetConformer(conformer).GetId()

        lines = []
        atom_counts = {}
        hetatm = 'HETATM'
        alt_loc = ''
        res_name = 'UNL'
        chain_id = ''
        res_seq = '1'
        i_code = ''
        occupancy = '1.00'
        temp_factor = '0.00'
        for atom in atoms:
            serial = atom+1
            element = self.atom_symbol(atom)
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f'{element}{atom_counts[element]}'
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = (
                f'{i}'[:8]
                for i in self.atom_coords(atom, conformer)
            )
            charge = self.mol.GetAtomWithIdx(atom).GetFormalCharge()
            lines.append(
                f'{hetatm:<6}{serial:>5} {name:<4}'
                f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
                f'{res_seq:>4}{i_code:<1}   '
                f'{x:>8}{y:>8}{z:>8}'
                f'{occupancy:>6}{temp_factor:>6}          '
                f'{element:>2}{charge:>2}\n'
            )

        # Convert to set because membership is going to be checked by
        # bonds.
        atoms = set(atoms)
        conect = 'CONECT'
        for bond in self.mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in atoms and a2 in atoms:
                lines.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        lines.append('END\n')
        with open(path, 'w') as f:
            f.write(''.join(lines))
