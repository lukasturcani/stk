"""
PdbWriter
==========

.. toctree::
    :maxdepth: 2

    PdbWriter <stk.molecular.writers.pdb>

"""

from .utilities import cell_matrix_to_lengths_angles


class PdbWriter:
    """
    A writer class for `.pdb` files.

    Examples
    --------
    *Writing to File with Unit Cell*

    This writer can write to file for visualisation with the unit cell
    included for periodic molecules. Note that this always assumes P1
    space group.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
        topology_graph = stk.cof.Honeycomb(
            building_blocks=(bb1, bb2),
            lattice_size=(3, 3, 1),
            periodic=True,
        )
        cof = stk.ConstructedMolecule(topology_graph)
        writer = stk.PdbWriter()
        writer.write(
            mol=cof,
            file='cof.pdb',
            periodic_cell=topology_graph.get_periodic_cell()
        )

    """

    def _write_string(self, mol, atom_ids, periodic_cell=None):

        if atom_ids is None:
            atom_ids = range(mol.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        lines = []
        if periodic_cell is not None:
            # Input unit cell information.
            lengths_and_angles = cell_matrix_to_lengths_angles(
                periodic_cell
            )
            a, b, c, alpha, beta, gamma = lengths_and_angles
            lines.append(
                f'CRYST1 {a:>8.3f} {b:>8.3f} {c:>8.3f}'
                f' {alpha:>6.2f} {beta:>6.2f} {gamma:>6.2f} '
                f'P 1\n'
            )

        atom_counts = {}
        hetatm = 'HETATM'
        alt_loc = ''
        res_name = 'UNL'
        chain_id = ''
        res_seq = '1'
        i_code = ''
        occupancy = '1.00'
        temp_factor = '0.00'

        coords = mol.get_position_matrix()
        # This set will be used by bonds.
        atoms = set()
        for atom_id in atom_ids:
            atoms.add(atom_id)
            atom, = mol.get_atoms(atom_ids=atom_id)
            serial = atom_id+1
            element = atom.__class__.__name__
            charge = atom.get_charge()
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f'{element}{atom_counts[element]}'
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = (i for i in coords[atom_id])

            lines.append(
                f'{hetatm:<6}{serial:>5} {name:<4}'
                f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
                f'{res_seq:>4}{i_code:<1}   '
                f' {x:>7.3f} {y:>7.3f} {z:>7.3f}'
                f'{occupancy:>6}{temp_factor:>6}          '
                f'{element:>2}{charge:>2}\n'
            )

        conect = 'CONECT'
        for bond in mol.get_bonds():
            a1 = bond.get_atom1().get_id()
            a2 = bond.get_atom2().get_id()
            # # Do not form periodic bonds in a PDB.
            # if periodic_cell is not None and bond.is_periodic():
            #     continue
            if a1 in atoms and a2 in atoms:
                lines.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        lines.append('END\n')

        return ''.join(lines)

    def write(self, mol, file=None, atom_ids=None, periodic_cell=None):
        """
        Write `mol` to `.pdb` file format.

        Parameters
        ----------
        mol : :class:`.Molecule`
            Molecule to write to `.pdb` format.

        file : :class:`str`, optional
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_cell : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        Returns
        -------
        string :class:`string`
            The string containing all information for a `.pdb` file.

        None : :class:`NoneType`
            A file is written if :attr:`file` is not `None`.

        """

        string = self._write_string(mol, atom_ids, periodic_cell)

        if file is None:
            return string
        else:
            with open(file, 'w') as f:
                f.write(string)
            return None
