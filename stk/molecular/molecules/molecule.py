import json
import os
import numpy as np
import rdkit.Chem.AllChem as rdkit
from scipy.spatial.distance import euclidean
from inspect import signature

from ...utilities import (
    normalize_vector,
    vector_theta,
    rotation_matrix,
    rotation_matrix_arbitrary_axis,
    mol_from_mae_file,
    remake,
    periodic_table
)


class MoleculeSubclassError(Exception):
    ...


class _Cached(type):
    def __call__(cls, *args, **kwargs):
        sig = signature(cls.__init__)
        sig = sig.bind_partial(cls, *args, **kwargs)
        sig.apply_defaults()
        key = cls._get_key(**sig.arguments)
        if sig.arguments['use_cache'] and key in cls._cache:
            return cls._cache[key]
        obj = super().__call__(*args, **kwargs)
        obj._key = key
        if sig.arguments['use_cache']:
            cls._cache[key] = obj
        return obj


class Molecule(metaclass=_Cached):
    """
    A base class for representing molecules.

    This class defines the operations which any class
    describing molecules should inherit or may find useful. Examples of
    such are :class:`.BuildingBlock` and :class:`.ConstructedMolecule`.
    This class should not be used directly, use an instance of a
    derived class instead.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms which compose the molecule.

    bonds : :class:`tuple` of :class:`.Bond`
        The bonds of the molecule.

    _cache : :class:`dict`
        This is a class attribute. It maps a :attr:`_key` to the
        :class:`.Molecule` instance with that :attr:`_key`.

    _position_matrix : :class:`numpy.ndarray`
        A ``(3, n)`` matrix holding the position of every atom in the
        molecule.

    _key : :class:`object`
        A hashable :class:`object`. This attribute will be the same
        for molecules of the same class, which have the same structure.
        A private method :meth:`_get_key` must be defined for
        each subclass of :class:`.Molecule` and it will be used to
        generate the :attr:`_key`.

    Methods
    -------
    :meth:`__init__`
    :meth:`apply_displacement`
    :meth:`apply_rotation_about_axis`
    :meth:`apply_rotation_between_vectors`
    :meth:`apply_rotation_to_minimize_theta`
    :meth:`get_atom_coords`
    :meth:`get_atom_distance`
    :meth:`get_center_of_mass`
    :meth:`get_centroid`
    :meth:`get_direction`
    :meth:`get_maximum_diameter`
    :meth:`get_plane_normal`
    :meth:`get_position_matrix`
    :meth:`is_identical`
    :meth:`set_centroid`
    :meth:`set_position_matrix`
    :meth:`dump`
    :meth:`load`
    :meth:`to_rdkit_mol`
    :meth:`update_cache`
    :meth:`update_from_rdkit_mol`
    :meth:`update_from_file`
    :meth:`write`

    """

    subclasses = {}

    def __init__(self, atoms, bonds, position_matrix):
        """
        Initialize a :class:`Molecule`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms which compose the molecule.

        bonds : :class:`tuple` of :class:`.Bond`
            The bonds of the molecule.

        position_matrix : :class:`numpy.ndarray`
            A ``(n, 3)`` matrix holding the position of every atom in
            the :class:`.Molecule`.

        """

        self.atoms = atoms
        self.bonds = bonds
        self._position_matrix = position_matrix.T

    @classmethod
    def _init_from_dict(self, mol_dict, use_cache=False):
        """
        Create a :class:`Molecule` from a :class:`dict`.

        The :class:`Molecule` returned has the class specified in
        `mol_dict`, not :class:`Molecule`. You can use this if you
        don't know what class the instance in `mol_dict` is or should
        be.

        Parameters
        ----------
        mol_dict : :class:`dict`
            A :class:`dict` holding the :class:`dict` representation
            of a molecule.

        use_cache : :class:`bool`, optional
            If ``True``, a new instance will not be made if a cached
            and identical one already exists, the one which already
            exists will be returned. If ``True`` and a cached,
            identical instance does not yet exist the created one will
            be added to the cache.

        Returns
        -------
        :class:`Molecule`
            The molecule represented by `mol_dict`.

        """

        c = self.subclasses[mol_dict['class']]
        return c._init_from_dict(mol_dict, use_cache=use_cache)

    def __init_subclass__(cls, **kwargs):
        if cls.__name__ in cls.subclasses:
            msg = 'Subclass with this name already exists.'
            raise MoleculeSubclassError(msg)
        cls.subclasses[cls.__name__] = cls
        cls._cache = {}
        super().__init_subclass__(**kwargs)

    def apply_displacement(self, displacement):
        """
        Shift the molecule by `displacement`.

        Parameters
        ----------
        displacement : :class:`numpy.ndarray`
            A displacement vector applied to the molecule.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        self._position_matrix += displacement[:, np.newaxis]
        return self

    def apply_rotation_about_axis(self, theta, axis, origin):
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

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        # Set the origin of the rotation to "origin".
        self.apply_displacement(-origin)
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)

        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        self._position_matrix = rot_mat @ self._position_matrix

        # Return the centroid of the molecule to the original position.
        self.apply_displacement(origin)
        return self

    def apply_rotation_between_vectors(self, start, target, origin):
        """
        Rotate the molecule by a rotation from `start` to `target`.

        Given two direction vectors, `start` and `target`, this method
        applies the rotation required transform `start` to `target`
        onto the molecule. The rotation occurs about the `origin`.

        For example, if the `start` and `target` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you
        can associate a geometric feature of the molecule with a
        vector, then the molecule can be rotated so that this vector is
        aligned with `target`. The defined vector can be virtually
        anything. This means that any geometric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            A vector which is to be rotated so that it transforms to
            the `target` vector.

        targert : :class:`numpy.ndarray`
            This array holds the vector, onto which `start` is rotated.

        origin : :class:`numpy.ndarray`
            The point about which the rotation occurs.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        target = normalize_vector(target)

        # Set the origin of the rotation to "origin".
        self.apply_displacement(-origin)
        rot_mat = rotation_matrix(start, target)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        self._position_matrix = rot_mat @ self._position_matrix

        # Restore original position.
        self.apply_displacement(origin)
        return self

    def apply_rotation_to_minimize_theta(
        self,
        start,
        target,
        axis,
        origin
    ):
        """
        Rotates the molecule to minimize angle between vectors.

        Note that this function will not necessarily overlay the
        `start` and `target` vectors. This is because the possible
        rotation is restricted to the `axis`.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            The vector which is rotated.

        target : :class:`numpy.ndarray`
            The vector which is stationary.

        axis : :class:`numpy.ndarray`
            The vector about which the rotation happens.

        origin : :class:`numpy.ndarray`
            The origin about which the rotation happens.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        # If the vector being rotated is not finite, exit. This is
        # probably due to a planar molecule.
        if not all(np.isfinite(x) for x in start):
            return self

        self.apply_displacement(-origin)

        # 1. First transform the problem.
        # 2. The rotation axis is set equal to the z-axis.
        # 3. Apply this transformation to all vectors in the problem.
        # 4. Take only the x and y components of `start` and `target`.
        # 5. Work out the angle between them.
        # 6. Apply that rotation along the original rotation axis.

        rotmat = rotation_matrix(axis, [0, 0, 1])
        tstart = np.dot(rotmat, start)
        tstart = np.array([tstart[0], tstart[1], 0])

        # If the `tstart` vector is 0 after these transformations it
        # means that it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], 1e-8):
            self.apply_displacement(origin)
            return self

        tend = np.dot(rotmat, target)
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

        rot_mat = rotation_matrix_arbitrary_axis(angle, axis)
        self._position_matrix = rot_mat @ self._position_matrix
        self.apply_displacement(origin)
        return self

    def get_atom_coords(self, atom_ids=None):
        """
        Yield the coordinates of atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of the atoms whose coordinates are desired.
            If ``None``, then the coordinates of all atoms will be
            yielded.

        Yields
        ------
        :class:`numpy.ndarray`
            An array holding the x, y and z coordinates of the
            next atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        coords = self._position_matrix[:, atom_ids].T
        for atom_coords in coords:
            yield atom_coords

    def get_atom_distance(self, atom1_id, atom2_id):
        """
        Return the distance between 2 atoms.

        This method does not account for the van der Waals radius of
        atoms.

        Parameters
        ----------
        atom1_id : :class:`int`
            The id of the first atom.

        atom2_id : :class:`int`
            The id of the second atom.

        Returns
        -------
        :class:`float`
            The distance between the first and second atoms.

        """

        distance = euclidean(
            u=self._position_matrix[:, atom1_id],
            v=self._position_matrix[:, atom2_id]
        )
        return float(distance)

    def get_center_of_mass(self, atom_ids=None):
        """
        Return the centre of mass of a group of atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            center of mass. If ``None``, then all atoms will be used.

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
        elif not isinstance(atom_ids, (list, tuple)):
            # Iterable gets used twice, once in get_atom_coords
            # and once in zip.
            atom_ids = list(atom_ids)

        center = 0
        total_mass = 0.
        coords = self.get_atom_coords(atom_ids)
        for atom_id, coord in zip(atom_ids, coords):
            mass = self.atoms[atom_id].mass
            total_mass += mass
            center += mass*coord
        return np.divide(center, total_mass)

    def get_centroid(self, atom_ids=None):
        """
        Return the centroid of a group of atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which are used to calculate the
            centroid. If ``None``, then all atoms will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of atoms specified by `atom_ids`.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        s = self._position_matrix[:, atom_ids].sum(axis=1)
        return s / len(atom_ids)

    def get_direction(self, atom_ids=None):
        """
        Return a vector of best fit through the molecule's atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            vector. If ``None``, then all atoms will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The vector of best fit.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        pos = self._position_matrix[:, atom_ids].T
        return normalize_vector(
            np.linalg.svd(pos - pos.mean(axis=0))[-1][0]
        )

    def get_maximum_diameter(self, atom_ids=None):
        """
        Get the maximum diameter in the molecule.

        This method does not account for the van der Waals radius of
        atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`
            The ids of atoms which are considered when looking for the
            maximum diamater. If ``None`` then all atoms in the
            molecule are considered.

        Returns
        -------
        :class:`float`
            The maximum diameter in the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        coords = self._position_matrix[:, atom_ids]
        return float(euclidean(coords.min(axis=1), coords.max(axis=1)))

    def get_plane_normal(self, atom_ids=None):
        """
        Return the normal to the plane of best fit.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            plane. If ``None``, then all atoms will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            Vector orthonormal to the plane of the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        pos = self._position_matrix[:, atom_ids].T
        centroid = self.get_centroid(atom_ids)
        return np.linalg.svd(pos - centroid)[-1][2, :]

    def get_position_matrix(self):
        """
        Return a matrix holding the atomic positions.

        Returns
        -------
        :class:`numpy.ndarray`
            The array has the shape ``(n, 3)``. Each row holds the
            x, y and z coordinates of an atom.

        """

        return self._position_matrix.T

    def is_identical(self, other):
        """
        Check if `other` is identical to `self`.

        Identical means that the two molecules are treated in the same
        way by every ``stk`` operation.

        Identical molecules do not have to have the same atomic
        positions, but must have the same atoms, bonds and charges.

        Parameters
        ----------
        other : :class:`Molecule`
            The molecule to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `other` is an identical molecule, from
            an ``stk`` perspective.

        """

        return (
            self.__class__ is other.__class__
            and self._key == other._key
        )

    def set_centroid(self, position, atom_ids=None):
        """
        Set the centroid of the molecule to `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            This array holds the position on which the centroid of the
            molecule is going to be placed.

        atom_ids : :class:`iterable` of :class:`int`
            The ids of atoms which should have their centroid set to
            `position`. If ``None`` then all atoms are used.

        Returns
        -------
        :class:`Molecule`
            The molecule is returned.

        """

        centroid = self.get_centroid(atom_ids=atom_ids)
        self.apply_displacement(position-centroid)
        return self

    def set_position_matrix(self, position_matrix):
        """
        Set the molecule's coordinates to those in `position_matrix`.

        Parameters
        ----------
        position_matrix : :class:`numpy.ndarray`
            A position matrix of the molecule. The shape of the matrix
            is ``(n, 3)``.

        Returns
        -------
        :class:`Molecule`
            The molecule is returned.

        """

        self._position_matrix = position_matrix.T
        return self

    def dump(self, path, include_attrs=None):
        """
        Write a :class:`dict` of the molecule to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file to which the  :class:`dict`
            should be written.

        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecule to be added to
            the :class:. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            json.dump(self._to_dict(include_attrs), f, indent=4)

    def _get_key(*args, **kwargs):
        raise NotImplementedError()

    @classmethod
    def load(cls, path, use_cache=False):
        """
        Create a :class:`Molecule` from a dump file.

        The :class:`Molecule` returned has the class specified in
        in the file, not :class:`Molecule`. You can use this if you
        don't know what class the instance in the loaded molecule is or
        should be.

        Parameters
        ----------
        path : :class:`str`
            The full path holding a dumped molecule.

        use_cache : :class:`bool`, optional
            If ``True``, a new instance will not be made if a cached
            and identical one already exists, the one which already
            exists will be returned. If ``True`` and a cached,
            identical instance does not yet exist the created one will
            be added to the cache.

        Returns
        -------
        :class:`Molecule`
            The molecule held in `path`.

        """

        with open(path, 'r') as f:
            mol_dict = json.load(f)

        return cls._init_from_dict(mol_dict, use_cache=use_cache)

    def _to_mdl_mol_block(self, atom_ids=None):
        """
        Return a V3000 mol block of the molecule.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        Returns
        -------
        :class:`str`
            The V3000 mol block representing the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        atom_lines = []
        # This set gets used by bonds.
        atoms = set()
        for i, atom_id in enumerate(atom_ids, 1):
            atoms.add(atom_id)

            x, y, z = self._position_matrix[:, atom_id]
            atom = self.atoms[atom_id]
            symbol = atom.__class__.__name__
            charge = atom.charge
            charge = f' CHG={charge}' if charge else ''
            atom_lines.append(
                'M  V30 {} {} {:.4f} {:.4f} {:.4f} 0{}\n'.format(
                    atom_id+1, symbol, x, y, z, charge
                )
            )
        atom_block = ''.join(atom_lines)
        num_atoms = i

        bond_lines = []
        for bond_idx, bond in enumerate(self.bonds):
            a1 = bond.atom1.id
            a2 = bond.atom2.id
            if a1 in atoms and a2 in atoms:
                # Keep bond ids if all bonds are getting written.
                if num_atoms == len(self.atoms):
                    bond_id = bond_idx
                else:
                    bond_id = len(bond_lines)

                bond_lines.append(
                    f'M  V30 {bond_id+1} {bond.order} {a1+1} {a2+1}\n'
                )

        num_bonds = len(bond_lines)
        bond_block = ''.join(bond_lines)
        return (
            '\n'
            '     RDKit          3D\n'
            '\n'
            '  0  0  0  0  0  0  0  0  0  0999 V3000\n'
            'M  V30 BEGIN CTAB\n'
            f'M  V30 COUNTS {num_atoms} {num_bonds} 0 0 0\n'
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

    def to_rdkit_mol(self):
        """
        Return a :mod:`rdkit` version of the molecule.

        Returns
        -------
        :class:`rdkit.Mol`
            The molecule in :mod:`rdkit` format.

        """

        mol = rdkit.EditableMol(rdkit.Mol())
        for atom in self.atoms:
            rdkit_atom = rdkit.Atom(atom.atomic_number)
            rdkit_atom.SetFormalCharge(atom.charge)
            mol.AddAtom(rdkit_atom)

        for bond in self.bonds:
            mol.AddBond(
                beginAtomIdx=bond.atom1.id,
                endAtomIdx=bond.atom2.id,
                order=rdkit.BondType(bond.order)
            )

        mol = mol.GetMol()
        rdkit_conf = rdkit.Conformer(len(self.atoms))
        for atom_id, atom_coord in enumerate(self._position_matrix.T):
            rdkit_conf.SetAtomPosition(atom_id, atom_coord)
        mol.AddConformer(rdkit_conf)
        return mol

    def update_cache(self):
        """
        Update attributes of the cached molecule.

        If there is no identical molecule in the cache, then this
        molecule is added.

        When using multiprocessing, modified copies of the original
        molecules are created. In order to ensure that the cached
        molecules have their attributes updated to the values of the
        copies, this method must be run on the copies.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self._key in self.__class__._cache:
            d = dict(vars(self))
            self.__class__._cache[self._key].__dict__ = d
        else:
            self.__class__._cache[self._key] = self

    def update_from_rdkit_mol(self, mol):
        """
        Updates the structure to match an :mod:`rdkit` molecule.

        Parameters
        ----------
        mol : :class:`rdkit.Mol`
            The :mod:`rdkit` molecule to use for the structure update.

        Returns
        -------
        :class:`.Molecule`
            The molecule is returned.

        """

        pos_mat = mol.GetConformer().GetPositions()
        self.set_position_matrix(pos_mat)
        return self

    def update_from_file(self, path):
        """
        Update the molecular structure from a file.

        Multiple file types are supported, namely:

        #. ``.mol``, ``.sdf`` - MDL V2000 and V3000 files
        #. ``.xyz`` - XYZ files
        #. ``.mae`` - Schrodinger Maestro files
        #. ``.coord`` - Turbomole files

        Parameters
        ----------
        path : :class:`str`
            The path to a molecular structure file holding updated
            coordinates for the :class:`.Molecule`.

        Returns
        -------
        :class:`.Molecule`
            The molecule.

        """

        update_fns = {
            '.mol': self._update_from_mol,
            '.sdf': self._update_from_mol,
            '.mae': self._update_from_mae,
            '.xyz': self._update_from_xyz,
            '.coord': self._update_from_turbomole
        }
        _, ext = os.path.splitext(path)
        update_fns[ext](path=path)
        return self

    def _update_from_mae(self, path):
        """
        Update the molecular structure to match an ``.mae`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mae`` file from which the structure
            should be updated.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.update_from_rdkit_mol(mol_from_mae_file(path))

    def _update_from_mol(self, path):
        """
        Update the molecular structure to match an ``.mol`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mol`` file from which the structure
            should be updated.

        Returns
        -------
        None : :class:`NoneType`

        """

        mol = remake(
            rdkit.MolFromMolFile(
                molFileName=path,
                sanitize=False,
                removeHs=False
            )
        )
        self.update_from_rdkit_mol(mol=mol)

    def _update_from_xyz(self, path):
        """
        Update the molecular structure to match a ``.xyz`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mol`` file from which the structure
            should be updated.

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

        with open(path, 'r') as f:
            atom_count, _, *content = f.readlines()

        # Check the atom count is correct.
        num_atoms = len(self.atoms)
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

            if element != self.atoms[i].__class__.__name__:
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
        new_coords = np.array(new_coords)
        self.set_position_matrix(new_coords)

    def _update_from_turbomole(self, path):
        """
        Update the structure from a Turbomole ``.coord`` file.

        Note that coordinates in ``.coord`` files are given in Bohr.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.coord`` file from which the
            structure should be updated.

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

        with open(path, 'r') as f:
            _, *content, __ = f.readlines()

        # Check the atom count is correct.
        num_atoms = len(self.atoms)
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

            if element != self.atoms[i].__class__.__name__:
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
        new_coords = np.array(new_coords)
        self.set_position_matrix(new_coords)

    def write(self, path, atom_ids=None):
        """
        Write a molecular structure file of the molecule.

        This function will write the format based on the extension
        of `path`.

        #. ``.mol``, ``.sdf`` - MDL V3000 file
        #. ``.xyz`` - XYZ file
        #. ``.pdb`` - PDB file

        Parameters
        ----------
        path : :class:`str`
            The `path` to which the molecule should be written.

        atom_ids : :class:`iterable` of :class:`int`, optional
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

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
        write_func(path, atom_ids)

    def _write_mdl_mol_file(self, path, atom_ids):
        """
        Write a V3000 ``.mol`` file of the molecule.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            f.write(self._to_mdl_mol_block(atom_ids))

    def _write_xyz_file(self, path, atom_ids):
        """
        Write a ``.xyz`` file of the molecule.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

        content = [0]
        for i, atom_id in enumerate(atom_ids, 1):
            x, y, z = self._position_matrix[:, atom_id]
            symbol = self.atoms[atom_id].__class__.__name__
            content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
        # Set first line to the atom_count.
        content[0] = f'{i}\n\n'

        with open(path, 'w') as xyz:
            xyz.write(''.join(content))

    def _write_pdb_file(self, path, atom_ids):
        """
        Write a ``.pdb`` file of the molecule.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. If ``None`` then all atoms
            are written.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atom_ids is None:
            atom_ids = range(len(self.atoms))

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

        coords = self._position_matrix
        # This set will be used by bonds.
        atoms = set()
        for atom in atom_ids:
            atoms.add(atom)

            serial = atom+1
            element = self.atoms[atom].__class__.__name__
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f'{element}{atom_counts[element]}'
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = (f'{i}'[:8] for i in coords[:, atom])
            lines.append(
                f'{hetatm:<6}{serial:>5} {name:<4}'
                f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
                f'{res_seq:>4}{i_code:<1}   '
                f'{x:>8}{y:>8}{z:>8}'
                f'{occupancy:>6}{temp_factor:>6}          '
                f'{element:>2}{self.atoms[atom].charge:>2}\n'
            )

        conect = 'CONECT'
        for bond in self.bonds:
            a1 = bond.atom1.id
            a2 = bond.atom2.id
            if a1 in atoms and a2 in atoms:
                lines.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        lines.append('END\n')
        with open(path, 'w') as f:
            f.write(''.join(lines))
