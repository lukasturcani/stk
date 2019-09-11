"""
Functional Groups
=================

The purpose of functional groups is to identify which atoms of a
building block are modified during construction. See
:class:`.FunctionalGroup` for more details.

See the documentation of :class:`.Reactor` to see how reactions between
functional groups are performed.

.. _`adding functional groups`:

Adding Functional Groups
------------------------

During initialization of a :class:`.BuildingBlock` instance, the
names of :class:`.FGType` instances are supplied. These are used to
create the :class:`.FunctionalGroup` instances
found in :attr:`.BuildingBlock.func_groups` and therefore the ones
used in the construction of :class:`.ConstructedMolecule` instances.
If you want to be able to use a new :class:`.FGType` with building
blocks in this way, simply add a new :class:`.FGType` instance to
:data:`fg_types`. This can be done in this file or externally.
For example, if you want make ``stk`` use a new :class:`.FGType`
without touching the source code you can just do

.. code-block:: python

    import stk

    new_fg_type = FGType(
        name='my_fg_type',
        func_group_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        deleter_smarts=['[$([H][N][H])]']*2
    )
    stk.fg_types[new_fg_type.name] = new_fg_type

    # You can now use the new_fg_type with BuildingBlocks.
    bb = stk.BuildingBlock('NCCCN', ['my_fg_type'])

Note that when adding SMARTS, if you want to make a SMARTS that targets
an atom in an environment, for example, a bromine connected to a
carbon::

    [$([Br][C])]

The atom you are targeting needs to be written first. The above SMARTS
works but::

    [$([C][Br])]

does not.

"""

import rdkit.Chem.AllChem as rdkit
from collections import Counter

from ..utilities import flatten


class FGType:
    """
    Creates :class:`.FunctionalGroup` instances.

    Attributes
    ----------
    name : :class:`str`
        A name for the :class:`.FGType`. Helps with identification.

    """

    def __init__(
        self,
        name,
        func_group_smarts,
        bonder_smarts,
        deleter_smarts
    ):
        """
        Initialize a :class:`FGType` instance.

        Parameters
        ----------
        name : :class:`str`
            A name for the :class:`.FGType`. Helps with identification.

        func_group_smarts : :class:`str`
            A SMARTS string which matches all atoms in a functional
            group.

        bonder_smarts : :class:`list` of :class:`str`
            A :class:`list` of SMARTS strings, each of which matches a
            single atom in a functional group. The matched atom is
            added to :attr:`.FunctionalGroup.bonders`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which appears in
            the same :class:`.FunctionalGroup`.

        deleter_smarts : :class:`list` of :class:`str`
            A :class:`list` of SMARTS strings, each of which matches a
            single atom in a functional group. The matched atom is
            added to :attr:`.FunctionalGroup.deleters`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which appears in
            the same :class:`.FunctionalGroup`.

        """

        self.name = name
        self._func_group = rdkit.MolFromSmarts(func_group_smarts)
        self._bonders = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(bonder_smarts).items()
        ]
        self._deleters = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(deleter_smarts).items()
        ]

    def get_functional_groups(self, mol):
        """
        Yield the functional groups in `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule which is to have functional groups
            identified.

        Yields
        ------
        :class:`.FunctionalGroup`
            A :class:`.FunctionalGroup` for every matched
            functional group in the `mol`.

        Examples
        --------
        .. code-block:: python

            import stk

            mol = stk.BuildingBlock('NCCN')
            amine = stk.FGType(
                func_group_smarts='[N]([H])[H]',
                bonder_smarts=['[$([N]([H])[H])]'],
                deleter_smarts=['[$([H][N][H])]']*2
            )

            # Make mol use amine functional groups during construction
            # of ConstructedMolecule.
            mol.func_groups = tuple(amine.get_functional_groups(mol))

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        func_groups = rdkit_mol.GetSubstructMatches(self._func_group)

        # All the bonder atoms, grouped by fg.
        bonders = [[] for i in range(len(func_groups))]

        for bonder, count in self._bonders:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(bonder)
            ))

            matched_bonders = [
                [aid for aid in fg if aid in matches]
                for fg in func_groups
            ]

            for fg_id, fg in enumerate(func_groups):
                bonders[fg_id].extend(matched_bonders[fg_id][:count])

        # All the deleter atoms, grouped by fg.
        deleters = [[] for i in range(len(func_groups))]
        for deleter, count in self._deleters:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(deleter)
            ))

            matched_deleters = [
                [aid for aid in fg if aid in matches]
                for fg in func_groups
            ]

            for fg_id, fg in enumerate(func_groups):
                deleters[fg_id].extend(matched_deleters[fg_id][:count])

        for atom_ids in zip(func_groups, bonders, deleters):
            fg, fg_bonders, fg_deleters = atom_ids
            yield FunctionalGroup(
                atoms=tuple(mol.atoms[id_] for id_ in fg),
                bonders=tuple(mol.atoms[id_] for id_ in fg_bonders),
                deleters=tuple(mol.atoms[id_] for id_ in fg_deleters),
                fg_type=self
            )

    def __repr__(self):
        func_group_smarts = rdkit.MolToSmarts(self._func_group)
        bonder_smarts = [
            rdkit.MolToSmarts(mol) for mol, _ in self._bonders
        ]
        deleter_smarts = [
            rdkit.MolToSmarts(mol) for mol, _ in self._deleters
        ]
        return (
            f'FGType(\n'
            f'    name={self.name!r}\n'
            f'    func_group_smarts={func_group_smarts!r},\n'
            f'    bonder_smarts={bonder_smarts!r},\n'
            f'    deleter_smarts={deleter_smarts!r}\n'
            ')'
        )

    def __str__(self):
        return f'FGType({self.name!r})'


class FunctionalGroup:
    """
    Represents a functional group in a molecule.

    Instances of this class should only be made via
    :meth:`.FGType.get_functional_groups`.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms in the functional group.

    bonders : :class:`tuple` of :class:`.Atom`
        The bonder atoms in the functional group. These are atoms which
        have bonds added during the construction of a
        :class:`.ConstructedMolecule`.

    deleters : :class:`tuple` of :class:`.Atom`
        The deleter atoms in the functional group. These are atoms
        which are deleted during construction of a
        :class:`.ConstructedMolecule`.

    fg_type : :class:`.FGType`
        The :class:`.FGType` instance which created the
        :class:`.FunctionalGroup`.

    """

    def __init__(self, atoms, bonders, deleters, fg_type):
        """
        Initialize a :class:`.FunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms in the functional group.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms in the functional group.

        fg_type : :class:`.FGType`
            The :class:`.FGType` instance which created the
            :class:`.FunctionalGroup`.

        """

        self.atoms = atoms
        self.bonders = bonders
        self.deleters = deleters
        self.fg_type = fg_type

    def clone(self, atom_map=None):
        """
        Return a clone.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold different :class:`.Atom`
            instances, then a :class:`dict` should be provided, which
            maps atoms in the current :class:`.FunctionalGroup` to the
            atoms which should be used in the clone. Only atoms which
            need to be remapped need to be present in the `atom_map`.

        Returns
        -------
        :class:`FunctionalGroup`
            A clone.

        Examples
        --------
        .. code-block:: python

            import stk

            c0 = stk.C(0)
            c1 = stk.C(1)
            bb = stk.BuildingBlock('NCCN', ['amine'])
            fg = bb.func_groups[0]
            a0, a1 = fg.atoms[:2]

            # fg_clone is a clone of fg, except that in all places
            # fg holds a0 fg_clone holds c0 and in all places where
            # fg holds a1, fg_clone holds c1.
            fg_clone = fg.clone({
                a0: c0
                a1: c1
            })

        """

        if atom_map is None:
            atom_map = {}

        return self.__class__(
            atoms=tuple(atom_map.get(a, a) for a in self.atoms),
            bonders=tuple(atom_map.get(a, a) for a in self.bonders),
            deleters=tuple(atom_map.get(a, a)for a in self.deleters),
            fg_type=self.fg_type
        )

    def get_atom_ids(self):
        """
        Yield the ids of :attr:`atoms` in order.

        Yields
        ------
        :class:`int`
            The id an :class:`.Atom` in :attr:`atoms`.

        """

        yield from (a.id for a in self.atoms)

    def get_bonder_ids(self):
        """
        Yield the ids of :attr:`bonders` in order.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom` in :attr:`bonders`.

        """
        yield from (a.id for a in self.bonders)

    def get_deleter_ids(self):
        """
        Yield the ids of :attr:`deleters` in order.

        Yields
        -------
        :class:`int`
            The id of an :class:`.Atom` in :attr:`deleters`.

        """

        yield from (a.id for a in self.deleters)

    def __repr__(self):
        atoms = list(self.atoms)
        if len(atoms) == 1:
            atoms.append(None)
        atoms = ', '.join(
            repr(atom) if atom is not None else '' for atom in atoms
        )

        bonders = list(self.bonders)
        if len(bonders) == 1:
            bonders.append(None)
        bonders = ', '.join(
            repr(atom) if atom is not None else '' for atom in bonders
        )

        deleters = list(self.deleters)
        if len(deleters) == 1:
            deleters.append(None)
        deleters = ', '.join(
            repr(atom) if atom is not None else '' for atom in deleters
        )

        return (
            f'FunctionalGroup(\n'
            f'    atoms=( {atoms} ), \n'
            f'    bonders=( {bonders} ), \n'
            f'    deleters=( {deleters} ), \n'
            f'    fg_type={self.fg_type}\n'
            ')'
        )

    def __str__(self):
        atoms = list(self.atoms)
        if len(atoms) == 1:
            atoms.append('')
        atoms = ', '.join(str(atom) for atom in atoms)

        bonders = list(self.bonders)
        if len(bonders) == 1:
            bonders.append('')
        bonders = ', '.join(str(atom) for atom in bonders)

        deleters = list(self.deleters)
        if len(deleters) == 1:
            deleters.append('')
        deleters = ', '.join(str(atom) for atom in deleters)

        return (
            f'FunctionalGroup(\n'
            f'    atoms=( {atoms} ), \n'
            f'    bonders=( {bonders} ), \n'
            f'    deleters=( {deleters} ), \n'
            f'    fg_type={self.fg_type}\n'
            ')'
        )


_fg_types = (

    FGType(
        name='amine',
        func_group_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        deleter_smarts=['[$([H][N][H])]']*2
    ),

    FGType(
        name='primary_amine',
        func_group_smarts=(
            '[N;$(NC);!$(NC=*);!$(NC#*);!$(NC(-[O,N,S]))]([H])[H]'
        ),
        bonder_smarts=['[$([N]([H])[H])]'],
        deleter_smarts=['[$([H][N][H])]']*2
    ),

    FGType(
        name='aldehyde',
        func_group_smarts='[C](=[O])[H]',
        bonder_smarts=['[$([C](=[O])[H])]'],
        deleter_smarts=['[$([O]=[C][H])]']
    ),

    FGType(
        name='carboxylic_acid',
        func_group_smarts='[C](=[O])[O][H]',
        bonder_smarts=['[$([C](=[O])[O][H])]'],
        deleter_smarts=[
            '[$([H][O][C](=[O]))]',
            '[$([O]([H])[C](=[O]))]'
        ]
    ),

    FGType(
        name='amide',
        func_group_smarts='[C](=[O])[N]([H])[H]',
        bonder_smarts=['[$([C](=[O])[N]([H])[H])]'],
        deleter_smarts=(
            ['[$([N]([H])([H])[C](=[O]))]'] +
            ['[$([H][N]([H])[C](=[O]))]']*2
        )
    ),

    FGType(
        name='thioacid',
        func_group_smarts='[C](=[O])[S][H]',
        bonder_smarts=['[$([C](=[O])[S][H])]'],
        deleter_smarts=[
            '[$([H][S][C](=[O]))]',
            '[$([S]([H])[C](=[O]))]'
        ]
    ),

    FGType(
        name='alcohol',
        func_group_smarts='[O][H]',
        bonder_smarts=['[$([O][H])]'],
        deleter_smarts=['[$([H][O])]']
    ),

    FGType(
        name='thiol',
        func_group_smarts="[S][H]",
        bonder_smarts=['[$([S][H])]'],
        deleter_smarts=['[$([H][S])]']
    ),

    FGType(
        name='bromine',
        func_group_smarts='*[Br]',
        bonder_smarts=['[$(*[Br])]'],
        deleter_smarts=['[$([Br]*)]']
    ),

    FGType(
        name='iodine',
        func_group_smarts='*[I]',
        bonder_smarts=['[$(*[I])]'],
        deleter_smarts=['[$([I]*)]']
    ),

    FGType(
        name='alkyne',
        func_group_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]([H])#[C])]'],
        deleter_smarts=['[$([H][C]#[C])]']
    ),

    FGType(
        name='terminal_alkene',
        func_group_smarts='[C]=[C]([H])[H]',
        bonder_smarts=['[$([C]=[C]([H])[H])]'],
        deleter_smarts=(
            ['[$([H][C]([H])=[C])]']*2 +
            ['[$([C](=[C])([H])[H])]']
        )
    ),

    FGType(
        name='boronic_acid',
        func_group_smarts='[B]([O][H])[O][H]',
        bonder_smarts=['[$([B]([O][H])[O][H])]'],
        deleter_smarts=(
            ['[$([O]([H])[B][O][H])]']*2 +
            ['[$([H][O][B][O][H])]']*2
        )
    ),

    # This amine functional group only deletes one of the
    # hydrogen atoms when a bond is formed.
    FGType(
        name='amine2',
        func_group_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        deleter_smarts=['[$([H][N][H])]']
    ),

    FGType(
        name='secondary_amine',
        func_group_smarts='[H][N]([#6])[#6]',
        bonder_smarts=[
            '[$([N]([H])([#6])[#6])]'
        ],
        deleter_smarts=['[$([H][N]([#6])[#6])]']
    ),

    FGType(
        name='diol',
        func_group_smarts='[H][O][#6]~[#6][O][H]',
        bonder_smarts=['[$([O]([H])[#6]~[#6][O][H])]']*2,
        deleter_smarts=['[$([H][O][#6]~[#6][O][H])]']*2
    ),

    FGType(
        name='difluorene',
        func_group_smarts='[F][#6]~[#6][F]',
        bonder_smarts=['[$([#6]([F])~[#6][F])]']*2,
        deleter_smarts=['[$([F][#6]~[#6][F])]']*2
    ),

    FGType(
        name='dibromine',
        func_group_smarts='[Br][#6]~[#6][Br]',
        bonder_smarts=['[$([#6]([Br])~[#6][Br])]']*2,
        deleter_smarts=['[$([Br][#6]~[#6][Br])]']*2
    ),

    FGType(
        name='alkyne2',
        func_group_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]#[C][H])]'],
        deleter_smarts=['[$([H][C]#[C])]', '[$([C](#[C])[H])]']
    ),

    FGType(
        name='ring_amine',
        func_group_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
        bonder_smarts=[
            '[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
            '[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
        ],
        deleter_smarts=(
                ['[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]']*2 +
                ['[$([H][#6](~[#6R1])~[#6][N]([H])[H])]']
        )
    ),

)

fg_types = {fg_type.name: fg_type for fg_type in _fg_types}
