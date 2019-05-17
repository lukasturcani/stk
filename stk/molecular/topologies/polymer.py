import numpy as np
import rdkit.Chem.AllChem as rdkit

from .base import Topology
from ...utilities import dedupe, add_fragment_props, remake


class Linear(Topology):
    """
    A class representing linear polymers.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string showing the repeating unit of the :class:`.Polymer`.
        For example, ``"AB"`` or ``"ABB"``. The building block with
        index ``0`` in :attr:`.MacroMolecule.building_blocks` is
        labelled as ``"A"`` while index ``1`` as ``"B"`` and so on.

    orientation : :class:`tuple` of :class:`float`
        For each character in the repeating unit, a value between ``0``
        and ``1`` (both inclusive) must be given in a :class:`list`. It
        indicates the probability that each monomer will have its
        orientation along the chain flipped.

    n : :class:`int`
        The number of repeating units which are used to make the
        polymer.

    ends : :class:`str`
        The string represents how the end groups of the polymer are
        treated. If ``'h'`` the functional groups at the end of the
        polymer are converted into hydrogem atoms. If ``'fg'`` they are
        kept as the original functional group.

    """

    def __init__(self, repeating_unit, orientation, n, ends='fg'):
        """
        Initializes a :class:`Linear` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string showing the repeating unit of the
            :class:`.Polymer`. For example, ``"AB"`` or ``"ABB"``. The
            building block with index ``0`` in
            :attr:`.MacroMolecule.building_blocks` is labelled as
            ``"A"`` while index ``1`` as ``"B"`` and so on.

        orientation : :class:`tuple` of :class:`float`
            For each character in the repeating unit, a value between
            ``0`` (inclusive) and ``1`` (inclusive) must be given.
            The values give the probability that each monomer is
            flipped by 180 degrees when being added to the chain. If
            ``0`` then the monomer is guaranteed to not flip. If ``1``
            it is guaranteed to flip. This allows the user to create
            head-to-head or head-to-tail chains, as well as chain with
            a preference for head-to-head or head-to-tail if a number
            between ``0`` and ``1`` is chosen.

        n : :class:`int`
            The number of repeating units which are used to make the
            polymer.

        ends : :class:`str`, optional
            The string represents how the end groups of the polymer are
            treated. If ``'h'`` the functional groups at the end of the
            polymer are converted into hydrogem atoms. If ``'fg'`` they
            are kept as the original functional group.

        """

        self.repeating_unit = repeating_unit
        self.orientation = tuple(orientation)
        self.n = n
        self.ends = ends
        super().__init__(track_fgs=False)

    def cleanup(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.ends == 'h':
            self.hygrogen_ends(macro_mol)

    def hygrogen_ends(self, macro_mol):
        """
        Removes all deleter atoms and adds hydrogens.

        In polymers, you want to replace the functional groups at the
        ends with hydrogen atoms.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        deleters = []
        for func_group in self.reactor.func_groups:
            deleters.extend(func_group.deleter_ids)

        emol = rdkit.EditableMol(macro_mol.mol)
        for atom_id in sorted(deleters, reverse=True):
            emol.RemoveAtom(atom_id)
        macro_mol.mol = remake(emol.GetMol())
        macro_mol.mol = rdkit.AddHs(macro_mol.mol, addCoords=True)

    def place_mols(self, macro_mol):
        """
        Places monomers side by side.

        The monomers are placed along the x-axis, so that the vector
        running between the functional groups is placed along the axis.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Make a map from monomer label to object.
        mapping = {}
        # Assign every monomer a label ("A", "B", "C", etc.).
        for label, monomer in zip(dedupe(self.repeating_unit),
                                  macro_mol.building_blocks):
            mapping[label] = monomer

        # Make string representing the entire polymer, not just the
        # repeating unit.
        polymer = self.repeating_unit*self.n
        # Get the direction for each monomer along the entire chain,
        # not just the repeating unit.
        dirs = self.orientation*self.n

        # Go through the repeating unit and place each monomer.
        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            bb = mapping[label]
            macro_mol.bb_counter.update([bb])
            original_position = bb.mol.GetConformer().GetPositions().T

            # Flip or not flip the monomer as given by the probability
            # in `mdir`.
            mdir = np.random.choice([1, -1], p=[mdir, 1-mdir])
            bb.set_orientation2([mdir, 0, 0])

            # The first building block should be placed at 0, the
            # others have positions calculated based on bb size.
            x_coord = self._x_position(macro_mol, bb) if i else 0
            monomer_mol = bb.set_bonder_centroid([x_coord, 0, 0])
            monomer_mol = rdkit.Mol(monomer_mol)

            bb_index = macro_mol.building_blocks.index(bb)
            add_fragment_props(monomer_mol, bb_index, i)

            # Check which functional group is at the back and which
            # one at the front.
            n_fgs = len(bb.func_groups)
            if n_fgs == 2:
                c1, c2 = list(bb.bonder_centroids())
                front = 1 if c1[0] < c2[0] else 0

            # Flag to see if the first functional group in
            # self.reactor.func_groups should be bonded.
            if i == 0:
                self.bond_first = True if n_fgs == 1 else False

            num_atoms = macro_mol.mol.GetNumAtoms()
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                                              monomer_mol)

            for fg in bb.func_groups:
                if n_fgs == 2:
                    id_ = 2*i + 1 if fg.id == front else 2*i
                elif len(self.reactor.func_groups) == 0:
                    id_ = 0
                else:
                    id_ = len(self.reactor.func_groups)

                func_group = fg.shifted_fg(id_, num_atoms)
                self.reactor.func_groups.append(func_group)

            bb.set_position_from_matrix(original_position)

    def bonded_fgs(self, macro_mol):
        """
        Yields functional groups to react.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Yields
        -------
        :class:`tuple` of :class:`int`
            Holds the ids of the functional groups set to react.

        """

        fgs = sorted(self.reactor.func_groups, key=lambda fg: fg.id)

        start = 0 if self.bond_first else 1
        for i in range(start, len(self.reactor.func_groups)-1, 2):
            yield fgs[i], fgs[i+1]

    def _x_position(self, macro_mol, bb):
        """
        Calculates the x coordinate on which to place `bb`.

        Does this by checking the most how for down the x axis
        `macro_mol` stretches and checking the distance between
        the minimum x position of `bb` and its centroid.
        It then tries to place `bb` about 3 A away from `macro_mol`.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        bb : :class:`.StructUnit`
            The building block to be added to `macro_mol`.

        Returns
        -------
        :class:`float`
            The x coordinate on which to place `bb`.

        """

        mm_max_x = max(macro_mol.all_atom_coords(),
                       key=lambda x: x[1][0])[1][0]
        bb_min_x = min(bb.all_atom_coords(),
                       key=lambda x: x[1][0])[1][0]
        bb_len = bb.centroid()[0] - bb_min_x
        return mm_max_x + bb_len + 3
