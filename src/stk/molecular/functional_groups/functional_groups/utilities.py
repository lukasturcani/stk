from __future__ import annotations

from typing import Iterable

from ...atoms import Atom


def get_atom_map(
    id_map: dict[int, int],
    atoms: Iterable[Atom],
) -> dict[int, Atom]:
    """
    Get an atom map from an `id_map`.

    The atom map maps the id of an existing atom to the atom it
    should be replaced by.

    Parameters:

        id_map:
            Maps the id of an atom to its new id.

        atoms:
            The atoms which should have their ids updated as
            specified by the `id_map`.

    """

    atom_map = {}
    for atom in atoms:
        atom_id = atom.get_id()
        new_id = id_map.get(atom_id, atom_id)
        atom_map[atom_id] = atom.with_id(new_id)
    return atom_map
