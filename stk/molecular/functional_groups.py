"""
Defines tools for dealing with functional groups.

See the documentation of :class:`.Reactor` to see how reactions between
functional groups are performed.

.. _`adding functional groups`:

Extending stk: Adding more functional groups.
---------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:`FGType` instance should be added to
:data:`functional_groups`, which is defined in this module.

Adding a new :class:`FGType` instance to :data:`functional_groups` will
allow :meth:`.TopologyGraph.construct` to connect the functional group
to all others during construction. In most cases, nothing except adding
this instance should be necessary in order to incorporate new
functional groups.

Note that when adding SMARTS, if you want to make a SMARTS that targets
an atom in an environment, for example, a bromine connected to a
carbon::

    [$([Br][C])]

The atom you are targeting needs to be written first. The above SMARTS
works but::

    [$([C][Br])]

does not.

If a new functional group is to connect to another functional group
with a bond other than a single, the names of the functional groups
should be added to :attr:`Reactor.bond_orders`, along with the desired
bond order.

"""

import rdkit.Chem.AllChem as rdkit

from ..utilities import flatten


class FGType:
    """
    Creates :class:`.FunctionalGroup` instances.

    Attributes
    ----------
    name : :class:`str`
        A name for the :class:`.FGType`. Helps with identification.

    """

    __slots__ = ['name', '_func_group', '_bonders', '_deleters']

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
            added to :class:`.FunctionalGroup.bonders`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which has bonds
            created during construction. The number of times the string
            needs to be repeated, is equal to number of atoms which
            need to be tagged.

        deleter_smarts : :class:`list` of :class:`str`
            A :class:`list` of SMARTS strings, each of which matches a
            single atom in a functional group. The matched atom is
            added to :class:`.FunctionalGroup.deleters`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which is deleted
            during construction. The number of times the string needs
            to be repeated, is equal to number of atoms which need to
            be tagged.

        """

        self.name = name
        self._func_group = rdkit.MolFromSmarts(func_group_smarts)
        self._bonders = [
            rdkit.MolFromSmarts(smarts) for smarts in bonder_smarts
        ]
        self._deleters = [
            rdkit.MolFromSmarts(smarts) for smarts in deleter_smarts
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
            A :class:`.FunctionalGroup` instance for every matched
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

        func_groups = mol.GetSubstructMatches(self._func_group)

        # All the bonder atoms, grouped by fg.
        bonders = [[] for i in range(len(func_groups))]

        for bonder in self._bonders:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(bonder)
            ))

            for fg_id, fg in enumerate(func_groups):
                for atom_id in fg:
                    if atom_id in matches:
                        bonders[fg_id].append(atom_id)
                        break

        # All the deleter atoms, grouped by fg.
        deleters = [[] for i in range(len(func_groups))]
        for deleter in self._deleters:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(deleter)
            ))

            for fg_id, fg in enumerate(func_groups):
                for atom_id in fg:
                    if atom_id in matches:
                        deleters[fg_id].append(atom_id)
                        break

        for atom_ids in zip(func_groups, bonders, deleters):
            fg, fg_bonders, fg_deleters = atom_ids
            yield FunctionalGroup(
                atoms=tuple(self.atoms[id_] for id_ in fg),
                bonders=tuple(self.atoms[id_] for id_ in fg_bonders),
                deleters=tuple(self.atoms[id_] for id_ in fg_deleters),
                fg_type=self
            )

    def __repr__(self):
        func_group_smarts = rdkit.MolToSmarts(self._func_group)
        bonder_smarts = [
            rdkit.MolToSmarts(mol) for mol in self._bonders
        ]
        deleter_smarts = [
            rdkit.MolToSmarts(mol) for mol in self._deleters
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
        return f'FGType({self.name})'


class FunctionalGroup:
    """
    Represents a functional group in a molecule.

    Instances of this class should only by made by using
    :class:`.FGType.get_functional_groups`.

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

    def __init__(self, atoms, bonders, deleters, fg_type=None):
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
        Get the ids of :attr:`atoms`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`atoms`.

        """

        return (a.id for a in self.atoms)

    def get_bonder_ids(self):
        """
        Get the ids of :attr:`bonders`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`bonders`.

        """
        return (a.id for a in self.bonders)

    def get_deleter_ids(self):
        """
        Get the ids of :attr:`deleters`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`deleters`.

        """

        return (a.id for a in self.deleters)

    def __repr__(self):
        atoms = list(self.atoms)
        if len(atoms) == 1:
            atoms.append('')
        atoms = ', '.join(repr(atom) if atom else '' for atom in atoms)

        bonders = list(self.bonders)
        if len(bonders) == 1:
            bonders.append('')
        bonders = ', '.join(
            repr(atom) if atom else '' for atom in bonders
        )

        deleters = list(self.deleters)
        if len(deleters) == 1:
            deleters.append('')
        deleters = ', '.join(
            repr(atom) if atom else '' for atom in deleters
        )

        return (
            f'FunctionalGroup(\n'
            f'    atoms=( {atoms} ), \n'
            f'    bonders=( {bonders} ), \n'
            f'    deleters=( {deleters} ), \n'
            f'    fg_type={self.fg_type.name}\n'
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
            f'    fg_type={self.fg_type.name}\n'
            ')'
        )


fg_types = {

    'amine': FGType(
        fg_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        del_smarts=['[$([H][N][H])]']*2
    ),

    'aldehyde': FGType(
        fg_smarts='[C](=[O])[H]',
        bonder_smarts=['[$([C](=[O])[H])]'],
        del_smarts=['[$([O]=[C][H])]']
    ),

    'carboxylic_acid': FGType(
        fg_smarts='[C](=[O])[O][H]',
        bonder_smarts=['[$([C](=[O])[O][H])]'],
        del_smarts=[
            '[$([H][O][C](=[O]))]',
            '[$([O]([H])[C](=[O]))]'
        ]
    ),

    'amide': FGType(
        fg_smarts='[C](=[O])[N]([H])[H]',
        bonder_smarts=['[$([C](=[O])[N]([H])[H])]'],
        del_smarts=(
            ['[$([N]([H])([H])[C](=[O]))]'] +
            ['[$([H][N]([H])[C](=[O]))]']*2
        )
    ),

    'thioacid': FGType(
        fg_smarts='[C](=[O])[S][H]',
        bonder_smarts=['[$([C](=[O])[S][H])]'],
        del_smarts=['[$([H][S][C](=[O]))]', '[$([S]([H])[C](=[O]))]']
    ),

    'alcohol': FGType(
        fg_smarts='[O][H]',
        bonder_smarts=['[$([O][H])]'],
        del_smarts=['[$([H][O])]']
    ),

    'thiol': FGType(
        fg_smarts="[S][H]",
        bonder_smarts=['[$([S][H])]'],
        del_smarts=['[$([H][S])]']
    ),

    'bromine': FGType(
        fg_smarts='*[Br]',
        bonder_smarts=['[$(*[Br])]'],
        del_smarts=['[$([Br]*)]']
    ),

    'iodine': FGType(
        fg_smarts='*[I]',
        bonder_smarts=['[$(*[I])]'],
        del_smarts=['[$([I]*)]']
    ),

    'alkyne': FGType(
        fg_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]([H])#[C])]'],
        del_smarts=['[$([H][C]#[C])]']
    ),

    'terminal_alkene': FGType(
        fg_smarts='[C]=[C]([H])[H]',
        bonder_smarts=['[$([C]=[C]([H])[H])]'],
        del_smarts=(
            ['[$([H][C]([H])=[C])]']*2 +
            ['[$([C](=[C])([H])[H])]']
        )
    ),

    'boronic_acid': FGType(
        fg_smarts='[B]([O][H])[O][H]',
        bonder_smarts=['[$([B]([O][H])[O][H])]'],
        del_smarts=(
            ['[$([O]([H])[B][O][H])]']*2 +
            ['[$([H][O][B][O][H])]']*2
        )
    ),

    # This amine functional group only deletes one of the
    # hydrogen atoms when a bond is formed.
    'amine2': FGType(
        fg_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        del_smarts=['[$([H][N][H])]']
    ),

    'secondary_amine': FGType(
        fg_smarts='[H][N]([#6])[#6]',
        bonder_smarts=[
            '[$([N]([H])([#6])[#6])]'
        ],
        del_smarts=['[$([H][N]([#6])[#6])]']
    ),

    'diol': FGType(
        fg_smarts='[H][O][#6]~[#6][O][H]',
        bonder_smarts=['[$([O]([H])[#6]~[#6][O][H])]']*2,
        del_smarts=['[$([H][O][#6]~[#6][O][H])]']*2
    ),

    'difluorene': FGType(
        fg_smarts='[F][#6]~[#6][F]',
        bonder_smarts=['[$([#6]([F])~[#6][F])]']*2,
        del_smarts=['[$([F][#6]~[#6][F])]']*2
    ),

    'dibromine': FGType(
        fg_smarts='[Br][#6]~[#6][Br]',
        bonder_smarts=['[$([#6]([Br])~[#6][Br])]']*2,
        del_smarts=['[$([Br][#6]~[#6][Br])]']*2
    ),

    'alkyne2': FGType(
        fg_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]#[C][H])]'],
        del_smarts=['[$([H][C]#[C])]', '[$([C](#[C])[H])]']
    ),

    'ring_amine': FGType(
        fg_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
        bonder_smarts=[
            '[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
            '[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
        ],
        del_smarts=(
                ['[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]']*2 +
                ['[$([H][#6](~[#6R1])~[#6][N]([H])[H])]']
        )
    ),

}
