"""
Turbomole Writer
================

"""

from __future__ import annotations

from collections import abc
import pathlib
import typing

from .. import molecule as _molecule
from .. import periodic_info as _periodic_info
from stk.utilities import typing as _typing


__all__ = (
    'TurbomoleWriter',
)


class TurbomoleWriter:
    """
    A writer class for ``Turbomole`` files.

    Examples:

        *Writing to a File with a Unit Cell*

        This writer can write to a file with the unit
        cell included for periodic molecules. Note that this always
        assumes P1 space group.

        .. testcode:: writing-to-a-file-with-a-unit-cell

            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
            bb2 = stk.BuildingBlock(
                smiles='BrCC(CBr)CBr',
                functional_groups=[stk.BromoFactory()],
            )
            topology_graph = stk.cof.PeriodicHoneycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 3, 1),
            )
            construction_result = topology_graph.construct()
            cof = (
                stk.ConstructedMolecule.init_from_construction_result(
                    construction_result=construction_result,
                )

            )
            writer = stk.TurbomoleWriter()
            writer.write(
                molecule=cof,
                path='cof.coord',
                periodic_info=construction_result.get_periodic_info(),
            )

        .. testcode:: writing-to-a-file-with-a-unit-cell
            :hide:

            import os

            assert os.path.exists('cof.coord')

        .. testcleanup:: writing-to-a-file-with-a-unit-cell

            os.remove('cof.coord')

    """

    def _write_content(
        self,
        molecule: _molecule.Molecule,
        atom_ids: typing.Optional[_typing.OneOrMany[int]],
        periodic_info: typing.Optional[_periodic_info.PeriodicInfo],
    ) -> abc.Iterable[str]:

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        if periodic_info is not None:
            # Input unit cell information.
            a = periodic_info.get_a()
            b = periodic_info.get_b()
            c = periodic_info.get_c()
            alpha = periodic_info.get_alpha()
            beta = periodic_info.get_beta()
            gamma = periodic_info.get_gamma()
            yield (
                '$periodic 3\n'
                '$cell angs\n'
                f' {a:>8.3f} {b:>8.3f} {c:>8.3f} '
                f'{alpha:>6.2f} {beta:>6.2f} {gamma:>6.2f}\n'
            )

        coords = molecule.get_position_matrix()
        yield '$coord angs\n'
        for atom_id in atom_ids:
            atom, = molecule.get_atoms(atom_ids=atom_id)
            element = atom.__class__.__name__
            x, y, z = (i for i in coords[atom_id])
            yield (
                f' {round(x, 4)} {round(y, 4)} {round(z, 4)} '
                f'{element}\n'
            )

        yield '$end\n'

    def to_string(
        self,
        molecule: _molecule.Molecule,
        atom_ids: typing.Optional[_typing.OneOrMany[int]] = None,
        periodic_info:
            typing.Optional[_periodic_info.PeriodicInfo] = None
    ) -> str:
        """
        Get a ``Turbomole`` file format string of `molecule`.

        Parameters:

            molecule:
                Molecule to write to ``Turbomole`` format.

            atom_ids:
                The atom ids of atoms to write. Can be a single
                :class:`int`, if a single atom is to be used, or
                ``None``, if all atoms are to be used.

            periodic_info:
                Information about the periodic cell.

        Returns:

            The content of a `.coord` file.

        """

        return ''.join(self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_info=periodic_info,
        ))

    def write(
        self,
        molecule: _molecule.Molecule,
        path: typing.Union[pathlib.Path, str],
        atom_ids: typing.Optional[_typing.OneOrMany[int]] = None,
        periodic_info:
            typing.Optional[_periodic_info.PeriodicInfo] = None,
    ) -> None:
        """
        Write `molecule` to ``Turbomole`` file format.

        Parameters:

            molecule:
                Molecule to write to ``Turbomole`` format.

            path:
                The full path to the file being written.

            atom_ids:
                The atom ids of atoms to write. Can be a single
                :class:`int`, if a single atom is to be used, or
                ``None``, if all atoms are to be used.

            periodic_info:
                Information about the periodic cell.

        Returns:

            A file is written.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_info=periodic_info,
        )

        with open(path, 'w') as f:
            f.write(''.join(content))
