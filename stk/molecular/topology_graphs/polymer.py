"""
Defines :class:`.Polymer` topologies.

"""

import numpy as np
import logging

from .topology_graph import TopologyGraph, Vertex, Edge


logger = logging.getLogger(__name__)


class LinearVertex(Vertex):
    def __init__(self, x, y, z, degree, direction):
        self._direction = np.array([direction, 0, 0])
        super().__init__(x, y, z, degree)

    def place_building_block(self, bb, conformer_id):
        if len(bb.func_groups) > 2:
            logger.warning(
                'You are placing a building block which has more than '
                'two functional groups along the backbone of '
                'a Linear topology. You can remove extra functional '
                'groups from the func_groups attribute to remove '
                'this message.'
            )

        bb.set_centroid(
            position=self._coord,
            atom_ids=bb.get_bonder_ids(fg_ids=(0, 1)),
            conformer_id=conformer_id
        )
        bonder_vector = next(bb.get_bonder_direction_vectors(
            fg_ids=(0, 1),
            conformer_id=conformer_id
        ))
        bb.apply_rotation_between_vectors(
            start=bonder_vector,
            end=self._direction,
            origin=self._coord,
            conformer_id=conformer_id
        )
        return bb.get_position_matrix(conformer_id=conformer_id)

    def _assign_vertex_positions(self, bb, conformer_id):
        fg1, fg2 = sorted(
            iterable=bb.func_groups,
            key=lambda fg: bb.get_centroid(
                atom_ids=fg.bonder_ids,
                conformer_id=conformer_id
            )[0]
        )
        self._positions[0].func_group = fg1
        self._positions[1].func_group = fg2


class TerminalVertex(LinearVertex):

    def __init__(self, x, y, z, direction):
        super().__init__(x, y, z, 1, direction)

    def place_building_block(self, bb, conformer_id):
        if len(bb.func_groups) != 1:
            return super().place_building_block(bb, conformer_id)

        bb.set_centroid(
            position=self._coord,
            atom_ids=bb.get_bonder_ids(fg_ids=(0, )),
            conformer_id=conformer_id
        )
        centroid_vector = next(
            bb.get_centroid_centroid_direction_vector(
                conformer_id=conformer_id
            )
        )
        bb.apply_rotation_between_vectors(
            start=centroid_vector,
            end=self._cap_direction,
            origin=self._coord
        )
        return bb.get_position_matrix(conformer_id=conformer_id)

    def _assign_vertex_positions(self, bb, conformer_id):
        if len(bb.func_groups) != 1:
            fg2, fg1 = sorted(
                iterable=bb.func_groups,
                key=lambda fg: bb.get_centroid(
                    atom_ids=fg.bonder_ids,
                    conformer_id=conformer_id
                )[0]
            )
        else:
            fg1 = bb.func_groups[0]

        self._positions[0].func_group = fg1


class HeadVertex(TerminalVertex):
    _cap_direction = [-1, 0, 0]


class TailVertex(TerminalVertex):
    _cap_direction = [1, 0, 0]


class Linear(TopologyGraph):
    """
    Represents linear polymer topology graphs.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string showing the repeating unit of the :class:`.Polymer`.
        For example, ``"AB"`` or ``"ABB"``. The building block with
        index ``0`` in :attr:`.ConstructedMolecule.building_blocks` is
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
            :attr:`.ConstructedMolecule.building_blocks` is labelled as
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

        head, *body, tail = orientation*n
        vertices = [_HeadVertex(0, 0, 0, head)]
        edges = []
        for i, direction in enumerate(body, 1):
            v = _LinearVertex(
                x=i, y=0, z=0, degree=2, direction=direction
            )
            vertices.append(v)
            p1 = vertices[i-1].positions[-1]
            p2 = vertices[i].positions[0]
            edges.append(Edge(p1, p2))

        vertices.append(_TailVertex(len(vertices), 0, 0, tail))
        p1 = vertices[-2].positions[-1]
        p2 = vertices[-1].positions[0]
        edges.append(Edge(p1, p2))

        super().__init__(vertices, edges)

    def _clean_up(self, mol):
        """
        Deletes the atoms which are lost during construction.

        Parameters
        ----------
        mol : :class:`.Polymer`
            The polymer being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.ends == 'h':
            self._hygrogen_ends(mol)

        super()._clean_up(mol)

    def _hygrogen_ends(self, mol):
        """
        Removes all deleter atoms and adds hydrogens.

        In polymers, you want to replace the functional groups at the
        ends with hydrogen atoms.

        Parameters
        ----------
        mol : :class:`.Polymer`
            The polymer being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def _get_bb_map(self, mol):
        raise NotImplementedError()

    def _get_conformer_map(self, mol):
        raise NotImplementedError()

    def _get_scale(self, mol, bb_map, conformer_map):
        maximum_diameter = 0
        for vertex, bb in bb_map.items():
            conformer_id = conformer_map[vertex]
            d = bb.get_maximum_diameter(conformer_id=conformer_id)
            if d > maximum_diameter:
                maximum_diameter = d
        return maximum_diameter
