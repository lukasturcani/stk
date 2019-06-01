"""
Defines macrocycle topologies.

"""

import numpy as np
import rdkit.Chem.AllChem as rdkit

from .base import Topology
from ...utilities import dedupe, add_fragment_props


class Cyclic(Topology):
    """
    Represents cyclic polymers forming a macrocycle.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string showing the repeating unit of the :class:`.Cyclic`.
        For example, ``"AB"`` or ``"ABB"``. The building block with
        index ``0`` in :attr:`.MacroMolecule.building_blocks` is
        labelled as ``"A"`` while index ``1`` as ``"B"`` and so on.

    orientation : :class:`tuple` of :class:`float`
        For each character in the repeating unit, a value between ``0``
        and ``1`` (both inclusive) must be given in a :class:`list`. It
        indicates the probability that each monomer will have its
        orientation along the chain flipped. If
        ``0`` then the monomer is guaranteed to not flip. If ``1``
        it is guaranteed to flip. This allows the user to create
        head-to-head or head-to-tail chains, as well as chain with
        a preference for head-to-head or head-to-tail if a number
        between ``0`` and ``1`` is chosen.

    n : :class:`int`
        The number of repeating units which in the macrocycle.

    """

    def __init__(self, repeating_unit, orientation, n):
        """
        Initialize a :class:`Cyclic` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string showing the repeating unit of the
            :class:`.Cyclic`. For example, ``"AB"`` or ``"ABB"``. The
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
            macrocycle.

        """

        self.repeating_unit = repeating_unit
        self.orientation = tuple(orientation)
        self.n = n
        super().__init__(track_fgs=False)

    def place_mols(self, macro_mol):
        """
        Place monomers on a large circle.

        The monomers are placed around a large circle centred at the
        origin, so that the vector running between the functional
        groups is placed in a direction tangent to the circle.

        Parameters
        ----------
        macro_mol : :class:`.Macrocycle`
            The macrocycle being assembled.

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

        # Calculate the radius of the circle so that the monomers
        # do not clash.
        md = max(x.max_diameter()[0]
                 for x in macro_mol.building_blocks)
        radius = (md + 2)/(2*np.sin(np.pi / len(polymer)))

        # Go through the repeating unit and place each monomer.
        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            bb = mapping[label]
            macro_mol.bb_counter.update([bb])
            original_position = bb.mol.GetConformer().GetPositions().T

            # Flip or not flip the monomer as given by the probability
            # in `mdir`.
            mdir = np.random.choice([1, -1], p=[mdir, 1-mdir])
            bb.set_orientation2([mdir, 0, 0])

            # Find the xyz coordinates of the monomer on the circle.
            theta = 2 * np.pi * i / len(polymer)
            x_coord = radius * np.cos(theta)
            y_coord = radius * np.sin(theta)

            # Check which functional group is at the back and which
            # one at the front.

            bb.set_orientation2([1, 0, 0])

            c1, c2 = list(bb.bonder_centroids())
            front = 1 if c1[0] < c2[0] else 0

            # Place the functional groups on the tangent to the circle.

            bb.rotate(theta + np.pi/2, [0, 0, 1])
            monomer_mol = bb.set_position([x_coord, y_coord, 0])
            monomer_mol = rdkit.Mol(monomer_mol)

            bb_index = macro_mol.building_blocks.index(bb)
            add_fragment_props(monomer_mol, bb_index, i)

            num_atoms = macro_mol.mol.GetNumAtoms()
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                                              monomer_mol)

            # Set ids for self.reactor.func_groups.
            for fg in bb.func_groups:
                id_ = 2*i + 1 if fg.id == front else 2*i

                func_group = fg.shifted_fg(id_, num_atoms)
                self.reactor.func_groups.append(func_group)

            bb.set_position_from_matrix(original_position)

    def bonded_fgs(self, macro_mol):
        """
        Yield functional groups to react.

        Parameters
        ----------
        macro_mol : :class:`.Macrocycle`
            The polymer being assembled.

        Yields
        -------
        :class:`tuple` of :class:`int`
            Holds the ids of the functional groups set to react.

        """

        fgs = sorted(self.reactor.func_groups, key=lambda fg: fg.id)

        for i in range(1, len(self.reactor.func_groups)-1, 2):
            yield fgs[i], fgs[i+1]

        yield fgs[0], fgs[-1]
