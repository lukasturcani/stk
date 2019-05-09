import numpy as np
import rdkit.Chem.AllChem as rdkit
import itertools as it
import math

from .struct_unit import StructUnit
from .macro_molecule import MacroMolecule
from ..functional_groups import Reactor
from ...utilities import flatten, vector_theta, Cell


class Periodic(MacroMolecule):
    """
    Represents periodic structures.

    This class is essentially the same as :class:`.MacroMolecule`,
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
        :class:`numpy.ndarray`, the element is an :class:`int` and the
        bond type is an :class:`rdkit.Chem.rdchem.BondType`.

    periodic_bonds : :class:`list` of :class:`.PeriodicBond`
        When periodic topologies are assembled, periodic bonds
        do not get added to the ``rdkit`` molecule in the
        :attr:`~.MacroMolecule.mol` attribute. Instead,
        :meth:`~.PeriodicLattice.join_mols` adds
        :class:`.PeriodicBond` instances representing the bonds into
        this list.

    cell_dimensions : :class:`list` of :class:`numpy.ndarray`
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
        :class:`rdkit.Mol`
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

        island : :class:`rdkit.Mol`
            The island molecule holding unit cells placed side by
            side like in a supercell but with no bonds running between
            them.

        Returns
        -------
        :class:`rdkit.Mol`
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

        reactor = Reactor(island)

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

                # `fg1` is found in `cell` and equivalent to
                # `periodic_bond.fg1`, having a bond added.
                fg1 = cell.fgs[periodic_bond.fg1]
                # `fg2`  is found in `ccell` and equivalent to
                # `periodic_bond.fg2`, having a bond added.
                fg2 = ccell.fgs[periodic_bond.fg2]

                reactor.react(fg1, fg2)

        return reactor.result(True)

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

                >>> cells, island = periodic._place_island([4, 4, 4])
                >>> cells[2][1][3]
                <Cell at 0x7fa0155d54e0>

            where the returned :class:`Cell` object represents the
            3rd unit cell along the x axis, the second along the y axis
            and the fourth along the z axis.

            The second member is an ``rdkit`` molecule of the island
            being built.

        """

        a, b, c = self.cell_dimensions
        cells = np.full(dimensions, None, object).tolist()
        island = rdkit.Mol()

        xdim, ydim, zdim = (range(d) for d in dimensions)
        nfgs = len(self.func_groups)

        # For each dimension place a unit cell.
        for i, (x, y, z) in enumerate(it.product(xdim, ydim, zdim)):
            unit_cell = self.shift(x*a + y*b + z*c)

            fgs = {}
            for fg in self.func_groups:
                id_ = fg.id + i*nfgs
                new_fg = fg.shifted_fg(id_, island.GetNumAtoms())
                fgs[fg] = new_fg

            cells[x][y][z] = Cell((x, y, z), fgs)
            island = rdkit.CombineMols(island, unit_cell)

        return cells, island

    def periodic_mol(self):
        """
        Creates a molecule by reacting periodic functional groups.

        Returns
        -------
        :class:`tuple`
            The first member is a :class:`rdkit.Mol`
            where the functional groups involved in a periodic bond
            have been reacted. This means the deleter atoms have been
            removed.

            The second member is :class:`list` of
            :class:`.AtomicPeriodicBond` holding the periodic bonds
            created as a result of the reactions.

        """

        reactor = Reactor()
        mol = rdkit.Mol(self.mol)
        periodic_bonds = []
        for pb in self.periodic_bonds:
            mol, _, new_bonds = reactor.periodic_react(
                                               mol,
                                               True,
                                               pb.direction,
                                               pb.fg1,
                                               pb.fg2)
            periodic_bonds.extend(new_bonds)

        return mol, periodic_bonds

    def write_gulp_input(self,
                         path,
                         keywords,
                         cell_fix=None,
                         atom_fix=None):
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

        atom_fix : :class:`numpy.ndarray` of :class:`int`, optional
            An n by 3 array where n is the number of atoms in the
            unit cell. Each row has the fix parameters for a given
            atom.

        Returns
        -------
        None : :class:`NoneType`

        """

        if cell_fix is None:
            cell_fix = [0, 0, 0, 0, 0, 0]

        pmol, pbs = self.periodic_mol()
        mol = StructUnit.__new__(StructUnit)
        mol.mol = pmol

        if atom_fix is None:
            atom_fix = np.ones([mol.mol.GetNumAtoms(), 3])

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
            for (id_, coords), fix in zip(mol.all_atom_coords(),
                                          atom_fix):

                x, y, z = [round(x, 4) for x in coords]
                fx, fy, fz = [int(x) for x in fix]
                f.write('{} core {} {} {} {} {} {}\n'.format(
                         mol.atom_symbol(id_), x, y, z, fx, fy, fz))
            f.write('\n')
            # Add bonds.
            for bond in mol.mol.GetBonds():
                a1 = bond.GetBeginAtomIdx() + 1
                a2 = bond.GetEndAtomIdx() + 1
                f.write('connect {} {} 0 0 0\n'.format(a1, a2))

            # Add periodic bonds.
            for bond in pbs:
                a1 = bond.atom1(mol.mol) + 1
                a2 = bond.atom2(mol.mol) + 1
                dx, dy, dz = bond.direction
                f.write(f'connect {a1} {a2} {dx:+} {dy:+} {dz:+}\n')
