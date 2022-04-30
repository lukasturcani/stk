"""
Molecule
========

.. toctree::
    :maxdepth: 2

    Building Block <stk.molecular.molecules.building_block>
    Constructed Molecule <stk.molecular.molecules.constructed_molecule>

"""


from __future__ import annotations

import os
import typing

import numpy as np
import rdkit.Chem.AllChem as rdkit
from scipy.spatial.distance import euclidean

from stk.utilities import (
    rotation_matrix,
    rotation_matrix_arbitrary_axis,
    vector_angle,
)

from ....utilities import OneOrMany
from ...atoms import Atom
from ...bonds import Bond
from ..utilities import get_bond_atom_ids, sort_bond_atoms_by_id
from .utilities import updaters, writers

_T = typing.TypeVar('_T', bound='Molecule')


class Molecule:
    """
    An abstract base class for molecules.

    Notes:

        You might notice that some of the methods of this abstract base
        class are implemented. This is purely for convenience when
        implementing subclasses. The implemented public methods are
        simply default implementations, which can be safely ignored
        or overridden, when implementing subclasses. Any private
        methods are implementation details of these default
        implementations.

    Examples:

        *Aligning a Molecule with a Vector*

        You want to rotate a molecule, such that it is aligned along
        with a specific direction.

        .. testcode:: aligning-a-molecule-with-a-vector

            import stk
            import numpy as np

            molecule1 = stk.BuildingBlock('CCCCC')
            # Align molecule1 along x-axis.
            molecule1 = molecule1.with_rotation_between_vectors(
                start=molecule1.get_direction(),
                target=np.array([1., 0., 0.]),
                origin=molecule1.get_centroid(),
            )

            molecule2 = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCCBr',
                            functional_groups=(stk.BromoFactory(), ),
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=15,
                ),
            )
            # Align molecule2 along the [1, 4, -3] vector.
            molecule2 = molecule2.with_rotation_between_vectors(
                start=molecule2.get_direction(),
                target=np.array([1., 4., -3.]),
                origin=molecule2.get_centroid(),
            )

        *Aligning a Molecule along a Plane*

        You want to place the benzene flat along the xy plane.

        .. testcode:: aligning-a-molecule-along-a-plane

            import stk
            import numpy as np

            benzene = stk.BuildingBlock('c1ccccc1')
            benzene = benzene.with_rotation_between_vectors(
                start=benzene.get_plane_normal(),
                target=np.array([0., 0., 1.]),
                origin=benzene.get_centroid(),
            )

    """

    def __init__(
        self,
        atoms: tuple[Atom, ...],
        bonds: tuple[Bond, ...],
        position_matrix: np.ndarray,
    ) -> None:
        """
        Initialize a :class:`.Molecule`.

        Parameters:

            atoms:
                The atoms which compose the molecule.

            bonds:
                The bonds of the molecule.

            position_matrix:
                A ``(n, 3)`` matrix holding the position of every atom
                in the :class:`.Molecule`.

        """

        self._atoms = atoms
        self._bonds = bonds
        # Take the transpose because it will make some matrix
        # multiplications faster.
        self._position_matrix = np.array(
            position_matrix.T,
            dtype=np.float64,
        )

    def _with_displacement(
        self: _T,
        displacement: np.ndarray,
    ) -> _T:
        """
        Modify molecule.

        """

        self._position_matrix = (
            self._position_matrix.T + displacement
        ).T
        return self

    def with_displacement(
        self,
        displacement: np.ndarray,
    ) -> Molecule:
        """
        Return a displaced clone.

        Parameters:

            displacement:
                The displacement vector to be applied.

        Returns:

            A displaced clone.

        """

        return self.clone()._with_displacement(displacement)

    def _with_rotation_about_axis(
        self: _T,
        angle: float,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> _T:
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

    def with_rotation_about_axis(
        self,
        angle: float,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> Molecule:
        """
        Return a rotated clone.

        The clone is rotated by `angle` about `axis` on the
        `origin`.

        Parameters:

            angle:
                The size of the rotation in radians.

            axis:
                The axis about which the rotation happens. Must have
                unit magnitude.

            origin:
                The origin about which the rotation happens.

        Returns:

            A rotated clone.

        """

        return self.clone()._with_rotation_about_axis(
            angle=angle,
            axis=axis,
            origin=origin,
        )

    def _with_rotation_between_vectors(
        self: _T,
        start: np.ndarray,
        target: np.ndarray,
        origin: np.ndarray,
    ) -> _T:
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

    def with_rotation_between_vectors(
        self,
        start: np.ndarray,
        target: np.ndarray,
        origin: np.ndarray,
    ) -> Molecule:
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
        can be easily aligned with any arbitrary direction.

        Parameters:

            start:
                A vector which is to be rotated so that it transforms
                into the `target` vector.

            target:
                The vector onto which `start` is rotated.

            origin:
                The point about which the rotation occurs.

        Returns:

            A rotated clone.

        """

        return self.clone()._with_rotation_between_vectors(
            start=start,
            target=target,
            origin=origin,
        )

    def _with_rotation_to_minimize_angle(
        self: _T,
        start: np.ndarray,
        target: np.ndarray,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> _T:

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
        tstart = start - np.dot(start, axis)*axis

        # If `tstart` is 0, it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], 1e-8):
            self._with_displacement(origin)
            return self

        tend = target - np.dot(target, axis)*axis
        # If `tend` is 0, it is parallel to the rotation axis, stop.
        if np.allclose(tend, [0, 0, 0], 1e-8):
            self._with_displacement(origin)
            return self

        angle = vector_angle(tstart, tend)

        projection = tstart @ np.cross(axis, tend)
        if projection > 0:
            angle = 2*np.pi - angle

        rotation_matrix = rotation_matrix_arbitrary_axis(angle, axis)
        self._position_matrix = rotation_matrix @ self._position_matrix
        self._with_displacement(origin)
        return self

    def with_rotation_to_minimize_angle(
        self,
        start: np.ndarray,
        target: np.ndarray,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> Molecule:
        """
        Return a rotated clone.

        The clone is rotated by the rotation required to minimize
        the angle between `start` and `target`.

        Note that this function will not necessarily overlay the
        `start` and `target` vectors. This is because the possible
        rotation is restricted to the `axis`.

        Parameters:

            start:
                The vector which is rotated.

            target:
                The vector which is stationary.

            axis:
                The vector about which the rotation happens. Must have
                unit magnitude.

            origin:
                The origin about which the rotation happens.

        Returns:

            A rotated clone.

        Raises:

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

    def _clone(self: _T) -> _T:
        clone = self.__class__.__new__(self.__class__)
        Molecule.__init__(
            self=clone,
            atoms=self._atoms,
            bonds=self._bonds,
            position_matrix=self._position_matrix.T,
        )
        return clone

    def clone(self) -> Molecule:
        """
        Return a clone.

        Returns:
            The clone.

        """

        return self._clone()

    def get_atomic_positions(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> typing.Iterator[np.ndarray]:
        """
        Yield the positions of atoms.

        Parameters:

            atom_ids:
                The ids of the atoms whose positions are desired.
                If ``None``, then the positions of all atoms will be
                yielded. Can be a single :class:`int`, if the position
                of a single atom is desired.

        Yields:

            The x, y and z coordinates of an atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        positions = self._position_matrix[:, atom_ids].T
        for atomic_position in positions:
            yield atomic_position

    def get_atoms(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> typing.Iterator[Atom]:
        """
        Yield the atoms in the molecule, ordered by id.

        Parameters:

            atom_ids:
                The ids of atoms to yield. Can be a single
                :class:`int` if a single atom is wanted, or ``None`` if
                all atoms are wanted.

        Yields:

            An atom in the molecule.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        for atom_id in atom_ids:
            yield self._atoms[atom_id]

    def get_num_atoms(self) -> int:
        """
        Return the number of atoms in the molecule.

        Returns:

            The number of atoms in the molecule.

        """

        return len(self._atoms)

    def get_bonds(self) -> typing.Iterator[Bond]:
        """
        Yield the bond in the molecule.

        Yields:

            A bond in the molecule.

        """

        for bond in self._bonds:
            yield bond

    def get_num_bonds(self) -> int:
        """
        Return the number of bonds in the molecule.

        Returns:

            The number of bonds in the molecule.

        """

        return len(self._bonds)

    def get_centroid(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> np.ndarray:
        """
        Return the centroid.

        Parameters:

            atom_ids:
                The ids of atoms which are used to calculate the
                centroid. Can be a single :class:`int`, if a single
                atom is to be used, or ``None`` if all atoms are to be
                used.

        Returns:

            The centroid of atoms specified by `atom_ids`.

        Raises:

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

    def get_direction(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> np.ndarray:
        """
        Return a vector of best fit through the atoms.

        Parameters:

            atom_ids:
                The ids of atoms which should be used to calculate the
                vector. Can be a single :class:`int`, if a single atom
                is to be used, or ``None``, if all atoms are to be
                used.

        Returns:

            The vector of best fit.

        Raises:

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
        return np.around(
            a=np.linalg.svd(pos - pos.mean(axis=0))[-1][0],
            decimals=14,
        )

    def get_maximum_diameter(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> float:
        """
        Return the maximum diameter.

        This method does not account for the van der Waals radius of
        atoms.

        Parameters:

            atom_ids:
                The ids of atoms which are considered when looking for
                the maximum diameter. Can be a single :class:`int`, if
                a single atom is to be used, or ``None``, if all atoms
                are to be used.

        Returns:

            The maximum diameter in the molecule.

        Raises:

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

    def get_plane_normal(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> np.ndarray:
        """
        Return the normal to the plane of best fit.

        Parameters:

            atom_ids:
                The ids of atoms which should be used to calculate the
                plane. Can be a single :class:`int`, if a
                single atom is to be used, or ``None``, if all atoms
                are to be used.

        Returns:

            Vector orthonormal to the plane of the molecule.

        Raises:

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
        return np.around(np.linalg.svd(pos - centroid)[-1][2, :], 14)

    def get_position_matrix(self) -> np.ndarray:
        """
        Return a matrix holding the atomic positions.

        Returns:

            The array has the shape ``(n, 3)``. Each row holds the
            x, y and z coordinates of an atom.

        """

        return np.array(self._position_matrix.T)

    def _with_centroid(
        self: _T,
        position: np.ndarray,
        atom_ids: typing.Optional[OneOrMany[int]],
    ) -> _T:

        centroid = self.get_centroid(atom_ids=atom_ids)
        self._with_displacement(position-centroid)
        return self

    def with_centroid(
        self,
        position: np.ndarray,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> Molecule:
        """
        Return a clone with its centroid at `position`.

        Parameters:

            position:
                This array holds the position on which the centroid of
                the clone is going to be placed.

            atom_ids:
                The ids of atoms which should have their centroid set
                to `position`. Can be a single :class:`int`, if a
                single atom is to be used, or ``None``, if all atoms
                are to be used.

        Returns:

            A clone with its centroid at `position`.

        """

        return self.clone()._with_centroid(position, atom_ids)

    def _with_position_matrix(
        self: _T,
        position_matrix: np.ndarray,
    ) -> _T:
        """
        Modify molecule.

        """

        self._position_matrix = np.array(position_matrix.T)
        return self

    def with_position_matrix(
        self,
        position_matrix: np.ndarray,
    ) -> Molecule:
        """
        Return a clone with atomic positions set by `position_matrix`.

        Parameters:

            position_matrix:
                The position matrix of the clone. The shape of the
                matrix is ``(n, 3)``.

        Returns:

            The clone.

        """
        return self.clone()._with_position_matrix(position_matrix)

    def to_rdkit_mol(self) -> rdkit.Mol:
        """
        Return an :mod:`rdkit` representation.

        Returns:

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
                order=(
                    rdkit.BondType.DATIVE if bond.get_order() == 9
                    else rdkit.BondType(bond.get_order())
                ),
            )

        mol = mol.GetMol()
        rdkit_conf = rdkit.Conformer(len(self._atoms))
        for atom_id, atom_coord in enumerate(self._position_matrix.T):
            rdkit_conf.SetAtomPosition(atom_id, atom_coord)
            mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
        mol.AddConformer(rdkit_conf)
        return mol

    def with_structure_from_file(
        self,
        path: str,
        extension: typing.Optional[str] = None,
    ) -> Molecule:
        """
        Return a clone, with its structure taken from a file.

        Multiple file types are supported, namely:

        #. ``.mol``, ``.sdf`` - MDL V2000 and V3000 files
        #. ``.xyz`` - XYZ files
        #. ``.mae`` - Schrodinger Maestro files
        #. ``.coord`` - Turbomole files
        #. ``.pdb`` - PDB files

        Parameters:

            path:
                The path to a molecular structure file holding updated
                coordinates for the :class:`.Molecule`.

            extension:
                If you want to treat the file as though it has a
                particular extension, put it here. Include the dot.

        Returns:

            A clone with atomic positions found in `path`.

        """

        if extension is None:
            _, extension = os.path.splitext(path)

        return {
            '.mol': updaters._with_structure_from_mol,
            '.sdf': updaters._with_structure_from_mol,
            '.mae': updaters._with_structure_from_mae,
            '.xyz': updaters._with_structure_from_xyz,
            '.coord': updaters._with_structure_from_turbomole,
            '.pdb': updaters._with_structure_from_pdb,
        }[extension](self.clone(), path)

    def with_canonical_atom_ordering(self) -> Molecule:
        """
        Return a clone, with canonically ordered atoms.

        Returns:

            The clone.

        """

        return self.clone()._with_canonical_atom_ordering()

    def _with_canonical_atom_ordering(self: _T) -> _T:
        """
        Modify the molecule.

        """

        atom_map = {
            atom.get_id(): atom.with_id(new_id)
            for new_id, atom in zip(
                rdkit.CanonicalRankAtoms(self.to_rdkit_mol()),
                self._atoms,
            )
        }
        self._atoms = tuple(sorted(
            atom_map.values(),
            key=lambda atom: atom.get_id()
        ))
        self._bonds = tuple(sorted(
            (
                sort_bond_atoms_by_id(bond.with_atoms(atom_map))
                for bond in self._bonds
            ),
            key=get_bond_atom_ids,
        ))
        old_ids = {
            atom.get_id(): old_id for old_id, atom in atom_map.items()
        }
        self._position_matrix = np.array(np.array([
            self._position_matrix.T[old_ids[new_id]]
            for new_id in range(len(self._atoms))
        ]).T)
        return self

    def get_canonical_atom_ids(self) -> dict[int, int]:
        """
        Map the id of each atom to its id under canonical ordering.

        Returns:

            Maps the id of each atom in the molecule to the id it would
            have under canonical ordering.

        """

        return {
            old_id: new_id
            for old_id, new_id in enumerate(
                rdkit.CanonicalRankAtoms(self.to_rdkit_mol()),
            )
        }

    def write(
        self,
        path: str,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> Molecule:
        """
        Write the structure to a file.

        This function will write the format based on the extension
        of `path`.

        #. ``.mol``, ``.sdf`` - MDL V3000 MOL file
        #. ``.xyz`` - XYZ file
        #. ``.pdb`` - PDB file

        Parameters:

            path:
                The `path` to which the molecule should be written.

            atom_ids:
                The atom ids of atoms to write. Can be a single
                :class:`int`, if a single atom is to be used, or
                ``None``, if all atoms are to be used. If you use this
                parameter, the atom ids in the file may not correspond
                to the atom ids in the molecule.

        Returns:

            The molecule.

        """

        _, extension = os.path.splitext(path)
        {
            '.mol': writers._write_mdl_mol_file,
            '.sdf': writers._write_mdl_mol_file,
            '.xyz': writers._write_xyz_file,
            '.pdb': writers._write_pdb_file,
        }[extension](self, path, atom_ids)
        return self

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__} at {id(self)}>'
