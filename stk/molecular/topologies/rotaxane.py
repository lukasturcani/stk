"""
Defines rotaxane topologies.

"""


import numpy as np
import rdkit.Chem.AllChem as rdkit

from .base import Topology
from ...utilities import dedupe, add_fragment_props


class NRotaxane(Topology):
    """
    Represents the topology of [n]rotaxanes.

    This class assumes one axle with (n-1) macrocycles threaded on it.
    The macrocycles are spaced evenly along the thread in repeating
    patterns analogous to non-bonded monomers in :class:`.Linear`. The
    orientation of the macrocycle defines the threading direction, thus
    giving access to different mechanical stereoisomers. The axle is
    automatically found in the :class:`list` of
    :attr:`.MacroMolecule.building_blocks` so only the order of the
    macrocycles in the list is important for construction.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string showing the repeating unit of the macrocycles within
        :class:`.Rotaxane`. For example, ``"AB"`` or ``"ABB"``, would
        implied two or three macrocycles threaded, respectively.
        The building block with index ``0`` in
        :attr:`.MacroMolecule.building_blocks` is labelled as ``"A"``
        while index ``1`` as ``"B"`` and so on.

    orientation : :class:`tuple` of :class:`float`
        For each character in the repeating unit, a value between ``0``
        and ``1`` (both inclusive) must be given in a :class:`list`. It
        indicates the probability that each macrocycle will have its
        orientation along the axle flipped.  If
        ``0`` then the macrocycle is guaranteed be aligned with the
        axle. If ``1`` it is guaranteed to be aligned against the
        axle. This allows the user to create stereoisomers.

    n : :class:`int`
        The number of repeating units in the rotaxane.

    """

    def __init__(self, repeating_unit, orientation, n):
        """
        Initialize a :class:`NRotaxane` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string showing the repeating unit of the macrocycles
            within :class:`.Rotaxane`. For example, ``"AB"`` or
            ``"ABB"``, would implied two or three macrocycles threaded,
            respectively. The building block with index ``0`` in
            :attr:`.MacroMolecule.building_blocks` is labelled as
            ``"A"`` while index ``1`` as ``"B"`` and so on.

        orientation : :class:`tuple` of :class:`float`
            For each character in the repeating unit, a value between
            ``0`` (inclusive) and ``1`` (inclusive) must be given.
            The values give the probability that each macrocycle is
            threaded onto the axle in the opposite directions. If
            ``0`` then the macrocycle is guaranteed be aligned with the
            axle. If ``1`` it is guaranteed to be aligned against the
            axle. This allows the user to create stereoisomers.

        n : :class:`int`
            The number of macrocyclerepeating units which are used to
            make rotaxane. Constructs [n*len(repeat_unit)+1]rotaxane.

        """

        self.repeating_unit = repeating_unit
        self.orientation = tuple(orientation)
        self.n = n
        super().__init__(track_fgs=False)

    def place_mols(self, macro_mol):
        """
        Distribute the macrocycles evenly along the axis.

        The axle is positioned along the x axis and the macrocycles
        are distributed evenly and rotates to that the vectors normal
        to their planes lie along the the axle.

        Parameters
        ----------
        macro_mol : :class:`.Rotaxane`
            The rotaxane being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Identify the axle and the macrocycles in the building_blocks.
        bbs = macro_mol.building_blocks
        axle = next(m for m in bbs if not hasattr(m, 'cycle_atoms'))

        cycles = [m for m in bbs if m is not axle]

        # Make a map from monomer label to object.
        mapping = {}
        # Assign every monomer a label ("A", "B", "C", etc.).
        for l, monomer in zip(dedupe(self.repeating_unit), cycles):
            mapping[l] = monomer

        # Make string representing the entire set of macrocycles,
        # not just the repeating unit.
        polycycle = self.repeating_unit*self.n

        # Get the direction for each macrocycle along the axle,
        # not just the repeating unit.
        dirs = self.orientation*self.n

        # Place the axle along the x axis with the centroid at origin.
        axle_dir = axle.direction()
        axle.set_orientation(axle_dir, [1, 0, 0])
        axle.set_position([0, 0, 0])

        macro_mol.bb_counter.update([axle])
        macro_mol.mol = rdkit.CombineMols(macro_mol.mol, axle.mol)

        # Find the limiting x coordinates to space the cycles evenly.
        min_x, max_x = self._minmax_x(axle)
        spacing = (max_x - min_x) / (len(polycycle)+1)

        # Space the macrocycles along the axle and rotate so that the
        # vectors normal to their planes lie along the x axis.

        for i, (label, mdir) in enumerate(zip(polycycle, dirs)):
            cycle = mapping[label]
            macro_mol.bb_counter.update([cycle])
            ring_ids = cycle.cycle_atoms()
            normal = cycle.plane_normal(atom_ids=ring_ids)
            org_pos = cycle.mol.GetConformer().GetPositions().T

            # Rotate the macrocycle towards the x or -x direction
            # as given by the probability in `mdir`.
            mdir = np.random.choice([1, -1], p=[mdir, 1-mdir])
            cycle.set_orientation(normal, [1, 0, 0])
            cycle.rotate(np.pi, [0, 0, 1]) if mdir < 0 else None

            # Position the macrocycle along the axle.
            cycle_x = min_x + (i + 1) * spacing
            mono_cyc = cycle.shift([cycle_x, 0, 0] -
                                   cycle.atom_centroid(ring_ids))

            cycle_index = macro_mol.building_blocks.index(cycle)
            add_fragment_props(mono_cyc, cycle_index, i)

            macro_mol.mol = rdkit.CombineMols(mono_cyc, macro_mol.mol)
            cycle.set_position_from_matrix(org_pos)

    def bonded_fgs(self, macro_mol):
        """
        Yield functional groups to react.

        Parameters
        ----------
        macro_mol : :class:`.Rotaxane`
            The polymer being assembled.

        Yields
        -------
        :class:`tuple` of :class:`int`
            Holds the ids of the functional groups set to react.

        """

        return iter(())

    def _minmax_x(self, axle, exclude_ids=None):
        """
        Calculate the maximum and minimum x coordinate along the axle.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        axle : :class:`.StructUnit`
            The axle of the rotaxane.

        Returns
        -------
        :class:`float`
            The minium x coordinate along the axle.

        :class:`float`
            The minium x coordinate along the axle.

        """

        conf = axle.mol.GetConformer()
        xyz = np.array(conf.GetPositions())

        if exclude_ids is not None:
            xyz = np.delete(xyz, exclude_ids, axis=0)

        max_x = np.amax(xyz, axis=0)[0]
        min_x = np.amin(xyz, axis=0)[0]

        return min_x, max_x
