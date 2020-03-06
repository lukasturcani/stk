import os
import numpy as np
from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit

from stk.utilities import (
    vector_angle,
    rotation_matrix,
    rotation_matrix_arbitrary_axis,
    periodic_table,
    mol_from_mae_file,
    remake,
)


class Molecule:
    """
    An abstract base class for all molecules.

    """

    def __init__(self, atoms, bonds, position_matrix):
        """
        Initialize a :class:`.Molecule`.

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

        self._atoms = atoms
        self._bonds = bonds
        # Take the transpose because it will make some matrix
        # multiplications faster.
        self._position_matrix = np.array(position_matrix.T)

    def _with_displacement(self, displacement):
        """
        Modify molecule.

        """

        self._position_matrix = (
            self._position_matrix.T + displacement
        ).T
        return self

    def with_displacement(self, displacement):
        """
        Return a displaced clone.

        Parameters
        ----------
        displacement : :class:`numpy.ndarray`
            The displacement vector applied.

        Returns
        -------
        :class:`.Molecule`
            A displaced clone.

        """

        return self.clone()._with_displacement(displacement)

    def _with_rotation_about_axis(self, angle, axis, origin):
        """
        Modify molecule.

        """

        # Set the origin of the rotation to "origin".
        self._with_displacement(-origin)
        rot_mat = rotation_matrix_arbitrary_axis(angle, axis)

        # Apply the rotation matrix on the position matrix, to get the
        # new position matrix.
        self._position_matrix = rot_mat @ self._position_matrix

        # Return the centroid of the molecule to the original position.
        self._with_displacement(origin)
        return self

    def with_rotation_about_axis(self, angle, axis, origin):
        """
        Return a rotated clone.

        The clone is rotated by `angle` about `axis` on the
        `origin`.

        Parameters
        ----------
        angle : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.ndarray`
            The axis about which the rotation happens.

        origin : :class:`numpy.ndarray`
            The origin about which the rotation happens.

        Returns
        -------
        :class:`.Molecule`
            A rotated clone.

        """

        return self.clone()._with_rotation_about_axis(
            angle=angle,
            axis=axis,
            origin=origin
        )

    def _with_rotation_between_vectors(self, start, target, origin):
        """
        Modify molecule.

        """

        # Set the origin of the rotation to "origin".
        self._with_displacement(-origin)
        rot_mat = rotation_matrix(start, target)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        self._position_matrix = rot_mat @ self._position_matrix

        # Restore original position.
        self._with_displacement(origin)
        return self

    def with_rotation_between_vectors(self, start, target, origin):
        """
        Return a rotated clone.

        The rotation is equal to a rotation from `start` to `target`.

        Given two direction vectors, `start` and `target`, this method
        applies the rotation required transform `start` to `target`
        onto the clone. The rotation occurs about the `origin`.

        For example, if the `start` and `target` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the clone. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you
        can associate a geometric feature of the molecule with a
        vector, then the clone can be rotated so that this vector is
        aligned with `target`. The defined vector can be virtually
        anything. This means that any geometric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            A vector which is to be rotated so that it transforms into
            the `target` vector.

        target : :class:`numpy.ndarray`
            The vector onto which `start` is rotated.

        origin : :class:`numpy.ndarray`
            The point about which the rotation occurs.

        Returns
        -------
        :class:`.Molecule`
            A rotated clone.

        """

        return self.clone()._with_rotation_between_vectors(
            start=start,
            target=target,
            origin=origin,
        )

    def _with_rotation_to_minimize_angle(
        self,
        start,
        target,
        axis,
        origin,
    ):
        # If the vector being rotated is not finite, exit. This is
        # probably due to a planar molecule.
        if not all(np.isfinite(x) for x in start):
            return self
        if np.allclose(target, [0, 0, 0], atol=1e-15):
            raise ValueError(
                'target has a magnitude of 0. It is therefore not '
                'possible to calculate an angle.'
            )

        self._with_displacement(-origin)

        # 1. Remove any component of the start and target vectors long
        # the axis. This puts them both on the same plane.
        # 2. Calculate the angle between them.
        # 3. Apply the rotation.

        axis_dot = np.dot(axis, axis)
        tstart = start - np.dot(start, axis)*axis/axis_dot

        # If `tstart` is 0, it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], 1e-8):
            self._with_displacement(origin)
            return self

        tend = target - np.dot(target, axis)*axis/axis_dot
        # If `tend` is 0, it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], 1e-8):
            self._with_displacement(origin)
            return self

        angle = vector_angle(tstart, tend)

        # Check in which direction the rotation should go.
        # This is done by applying the rotation in each direction and
        # seeing which one leads to a smaller angle.
        r1 = rotation_matrix_arbitrary_axis(angle, axis)
        t1 = vector_angle(np.dot(r1, tstart), tend)
        r2 = rotation_matrix_arbitrary_axis(-angle, axis)
        t2 = vector_angle(np.dot(r2, tstart), tend)

        if t1 < t2:
            self._position_matrix = r1 @ self._position_matrix
        else:
            self._position_matrix = r2 @ self._position_matrix

        self._with_displacement(origin)
        return self

    def with_rotation_to_minimize_angle(
        self,
        start,
        target,
        axis,
        origin,
    ):
        """
        Return a rotated clone.

        The clone is rotated by the rotation required to minimize
        the angle between `start` and `target`.

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
            A rotated clone.

        Raises
        ------
        :class:`ValueError`
            If `target` has a magnitude of 0. In this case it is
            not possible to calculate an angle between `start` and
            `target`.

        """

        return self.clone()._with_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=axis,
            origin=origin,
        )

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Molecule`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        Molecule.__init__(
            self=clone,
            atoms=self._atoms,
            bonds=self._bonds,
            position_matrix=self._position_matrix.T,
        )
        return clone

    def get_atomic_positions(self, atom_ids=None):
        """
        Yield the positions of atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of the atoms whose positions are desired.
            If ``None``, then the positions of all atoms will be
            yielded. Can be a single :class:`int`, if the position of
            a single atom is desired.

        Yields
        ------
        :class:`numpy.ndarray`
            The x, y and z coordinates of an atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        coords = self._position_matrix[:, atom_ids].T
        for atom_coords in coords:
            yield atom_coords

    def get_atoms(self, atom_ids=None):
        """
        Yield the atoms in the molecule, ordered by id.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms to yield. Can be a single
            :class:`int` if a single atom is wanted, or ``None`` if
            all atoms are wanted.

        Yields
        ------
        :class:`.Atom`
            An atom in the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        for atom_id in atom_ids:
            yield self._atoms[atom_id]

    def get_num_atoms(self):
        """
        Return the number of atoms in the molecule.

        Returns
        -------
        :class:`int`
            The number of atoms in the molecule.

        """

        return len(self._atoms)

    def get_bonds(self):
        """
        Yield the bond in the molecule.

        Yields
        ------
        :class:`.Bond`
            A bond in the molecule.

        """

        for bond in self._bonds:
            yield bond

    def get_num_bonds(self):
        """
        Return the number of bonds in the molecule.

        Returns
        -------
        :class:`int`
            The number of bonds in the molecule.

        """

        return len(self._bonds)

    def get_centroid(self, atom_ids=None):
        """
        Return the centroid.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which are used to calculate the
            centroid. Can be a single :class:`int`, if a single
            atom is to be used, or ``None`` if all atoms are to be
            used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of atoms specified by `atom_ids`.

        Raises
        ------
        :class:`ValueError`
            If `atom_ids` has a length of ``0``.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        if len(atom_ids) == 0:
            raise ValueError('atom_ids was of length 0.')

        return np.divide(
            self._position_matrix[:, atom_ids].sum(axis=1),
            len(atom_ids)
        )

    def get_direction(self, atom_ids=None):
        """
        Return a vector of best fit through the atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            vector. Can be a single :class:`int`, if a single atom
            is to be used, or ``None``, if all atoms are to be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The vector of best fit.

        Raises
        ------
        :class:`ValueError`
            If `atom_ids` has a length of ``0``.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        if len(atom_ids) == 0:
            raise ValueError('atom_ids was of length 0.')

        pos = self._position_matrix[:, atom_ids].T
        return np.linalg.svd(pos - pos.mean(axis=0))[-1][0]

    def get_maximum_diameter(self, atom_ids=None):
        """
        Return the maximum diameter.

        This method does not account for the van der Waals radius of
        atoms.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which are considered when looking for the
            maximum diameter. Can be a single :class:`int`, if a
            single atom is to be used, or ``None``, if all atoms are to
            be used.

        Returns
        -------
        :class:`float`
            The maximum diameter in the molecule.

        Raises
        ------
        :class:`ValueError`
            If `atom_ids` has a length of ``0``.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        if len(atom_ids) == 0:
            raise ValueError('atom_ids was of length 0.')

        coords = self._position_matrix[:, atom_ids]
        return float(euclidean(coords.min(axis=1), coords.max(axis=1)))

    def get_plane_normal(self, atom_ids=None):
        """
        Return the normal to the plane of best fit.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should be used to calculate the
            plane. Can be a single :class:`int`, if a
            single atom is to be used, or ``None``, if all atoms are to
            be used.

        Returns
        -------
        :class:`numpy.ndarray`
            Vector orthonormal to the plane of the molecule.

        Raises
        ------
        :class:`ValueError`
            If `atom_ids` has a length of ``0``.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        if len(atom_ids) == 0:
            raise ValueError('atom_ids was of length 0.')

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

        return np.array(self._position_matrix.T)

    def _with_centroid(self, position, atom_ids):
        centroid = self.get_centroid(atom_ids=atom_ids)
        self._with_displacement(position-centroid)
        return self

    def with_centroid(self, position, atom_ids=None):
        """
        Return a clone with its centroid at `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            This array holds the position on which the centroid of the
            clone is going to be placed.

        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which should have their centroid set to
            `position`. Can be a single :class:`int`, if a
            single atom is to be used, or ``None``, if all atoms are to
            be used.

        Returns
        -------
        :class:`Molecule`
            A clone with its centroid at `position`.

        """

        return self.clone()._with_centroid(position, atom_ids)

    def _with_position_matrix(self, position_matrix):
        """
        Modify molecule.

        """

        self._position_matrix = np.array(position_matrix.T)
        return self

    def with_position_matrix(self, position_matrix):
        """
        Return a clone with atomic positions set by `position_matrix`.

        Parameters
        ----------
        position_matrix : :class:`numpy.ndarray`
            A position matrix of the clone. The shape of the matrix
            is ``(n, 3)``.

        Returns
        -------
        :class:`Molecule`
            The clone.

        """
        return self.clone()._with_position_matrix(position_matrix)

    def _to_mdl_mol_block(self, atom_ids=None):
        """
        Return a V3000 mol block of the molecule.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or
            ``None``, if all atoms are to be used.

        Returns
        -------
        :class:`str`
            The V3000 mol block representing the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        atom_lines = []
        # This set gets used by bonds.
        atoms = set()
        for i, atom_id in enumerate(atom_ids, 1):
            atoms.add(atom_id)

            x, y, z = self._position_matrix[:, atom_id]
            atom = self._atoms[atom_id]
            symbol = atom.__class__.__name__
            charge = atom.get_charge()
            charge = f' CHG={charge}' if charge else ''
            atom_lines.append(
                'M  V30 {} {} {:.4f} {:.4f} {:.4f} 0{}\n'.format(
                    atom_id+1, symbol, x, y, z, charge
                )
            )
        atom_block = ''.join(atom_lines)
        num_atoms = i

        bond_lines = []
        for bond_idx, bond in enumerate(self._bonds):
            a1 = bond.get_atom1().get_id()
            a2 = bond.get_atom2().get_id()
            if a1 in atoms and a2 in atoms:
                # Keep bond ids if all bonds are getting written.
                if num_atoms == len(self._atoms):
                    bond_id = bond_idx
                else:
                    bond_id = len(bond_lines)

                bond_lines.append(
                    f'M  V30 {bond_id+1} '
                    f'{int(bond.get_order())} {a1+1} {a2+1}\n'
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
        Return an :mod:`rdkit` representation.

        Returns
        -------
        :class:`rdkit.Mol`
            The molecule in :mod:`rdkit` format.

        """

        mol = rdkit.EditableMol(rdkit.Mol())
        for atom in self._atoms:
            rdkit_atom = rdkit.Atom(atom.get_atomic_number())
            rdkit_atom.SetFormalCharge(atom.get_charge())
            mol.AddAtom(rdkit_atom)

        for bond in self._bonds:
            mol.AddBond(
                beginAtomIdx=bond.get_atom1().get_id(),
                endAtomIdx=bond.get_atom2().get_id(),
                order=rdkit.BondType(bond.get_order())
            )

        mol = mol.GetMol()
        rdkit_conf = rdkit.Conformer(len(self._atoms))
        for atom_id, atom_coord in enumerate(self._position_matrix.T):
            rdkit_conf.SetAtomPosition(atom_id, atom_coord)
            mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
        mol.AddConformer(rdkit_conf)
        return mol

    def with_structure_from_file(self, path, extension=None):
        """
        Return a clone, with its structure taken from a file.

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

        extension : :class:`str, optional
            If you want to treat the file as though it has a
            particular extension, put it here. Include the dot.

        Returns
        -------
        :class:`.Molecule`
            A clone with atomic positions found in `path`.

        """

        if extension is None:
            _, extension = os.path.splitext(path)

        clone = self.clone()
        update_fns = {
            '.mol': clone._with_structure_from_mol,
            '.sdf': clone._with_structure_from_mol,
            '.mae': clone._with_structure_from_mae,
            '.xyz': clone._with_structure_from_xyz,
            '.coord': clone._with_structure_from_turbomole,
        }
        return update_fns[extension](path=path)

    def _with_structure_from_mae(self, path):
        """
        Change structure to match an ``.mae`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mae`` file from which the structure
            should be updated.

        Returns
        -------
        :class:`.Molecule`
            The molecule.

        """

        molecule = mol_from_mae_file(path)
        return self._with_position_matrix(
            position_matrix=molecule.GetConformer().GetPositions()
        )

    def _with_structure_from_mol(self, path):
        """
        Change structure to match a ``.mol`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mol`` file from which the structure
            should be updated.

        Returns
        -------
        :class:`.Molecule`
            The molecule.

        """

        molecule = remake(
            rdkit.MolFromMolFile(
                molFileName=path,
                sanitize=False,
                removeHs=False,
            )
        )
        return self._with_position_matrix(
            position_matrix=molecule.GetConformer().GetPositions()
        )

    def _with_structure_from_xyz(self, path):
        """
        Return a clone, with its structure taken from an ``.xyz`` file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the ``.mol`` file from which the structure
            should be updated.

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

        with open(path, 'r') as f:
            atom_count, _, *content = f.readlines()

        # Check the atom count is correct.
        num_atoms = len(self._atoms)
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

            if element != self._atoms[i].__class__.__name__:
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
        return self._with_position_matrix(new_coords)

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

    def write(self, path, atom_ids=None):
        """
        Write the structure to a file.

        This function will write the format based on the extension
        of `path`.

        #. ``.mol``, ``.sdf`` - MDL V3000 MOL file
        #. ``.xyz`` - XYZ file
        #. ``.pdb`` - PDB file

        Parameters
        ----------
        path : :class:`str`
            The `path` to which the molecule should be written.

        atom_ids : :class:`iterable` of :class:`int`, optional
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or
            ``None``, if all atoms are to be used.

        Returns
        -------
        :class:`.Molecule`
            The molecule.

        """

        writers = {
            '.mol': self._write_mdl_mol_file,
            '.sdf': self._write_mdl_mol_file,
            '.xyz': self._write_xyz_file,
            '.pdb': self._write_pdb_file,
        }

        _, ext = os.path.splitext(path)
        write_func = writers[ext]
        write_func(path, atom_ids)
        return self

    def _write_mdl_mol_file(self, path, atom_ids):
        """
        Write to a V3000 ``.mol`` file.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            f.write(self._to_mdl_mol_block(atom_ids))

    def _write_xyz_file(self, path, atom_ids):
        """
        Write to a ``.xyz`` file.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        content = [0]
        for i, atom_id in enumerate(atom_ids, 1):
            x, y, z = self._position_matrix[:, atom_id]
            symbol = self._atoms[atom_id].__class__.__name__
            content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
        # Set first line to the atom_count.
        content[0] = f'{i}\n\n'

        with open(path, 'w') as xyz:
            xyz.write(''.join(content))

    def _write_pdb_file(self, path, atom_ids):
        """
        Write to a ``.pdb`` file.

        This function should not be used directly, only via
        :meth:`write`.

        Parameters
        ----------
        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

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
            element = self._atoms[atom].__class__.__name__
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
                f'{element:>2}{self._atoms[atom].get_charge():>2}\n'
            )

        conect = 'CONECT'
        for bond in self._bonds:
            a1 = bond.get_atom1().get_id()
            a2 = bond.get_atom2().get_id()
            if a1 in atoms and a2 in atoms:
                lines.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        lines.append('END\n')
        with open(path, 'w') as f:
            f.write(''.join(lines))
