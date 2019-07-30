"""
Defines topologies of polymers.

"""

import logging
import re
import numpy as np

from .topology_graph import TopologyGraph, Vertex, Edge


logger = logging.getLogger(__name__)


class _LinearVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def __init__(self, x, y, z, orientation):
        """
        Initialize a :class:`.LinearVertex`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        orientation : :class:`float`
            Can be any number from ``0`` to ``1``, both inclusive. It
            specifies the probability the building block placed on the
            vertex will have its orientation along the chain flipped.

        """

        self._orientation = orientation
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            If ``True`` the :attr:`edges` attribute of the clone will
            be empty.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone._orientation = self._orientation
        return clone

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        `building_block` is placed such that its bonder-bonder
        direction vector is either parallel or anti-parallel to the
        polymer chain.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        if len(building_block.func_groups) > 2:
            logger.warning(
                'You are placing a building block which has more than '
                'two functional groups along the backbone of '
                'a Linear topology. You can remove extra functional '
                'groups from the func_groups attribute to remove '
                'this message.'
            )

        building_block.set_centroid(
            position=self._coord,
            atom_ids=building_block.get_bonder_ids(fg_ids=(0, 1))
        )
        bonder_vector = next(
            building_block.get_bonder_direction_vectors(
                fg_ids=(0, 1)
            )
        )[-1]

        p = [1-self._orientation, self._orientation]
        direction = np.random.choice([1, -1], p=p)
        building_block.apply_rotation_between_vectors(
            start=bonder_vector,
            target=[direction, 0, 0],
            origin=self._coord
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        """
        Assign functional groups to edges.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

        fg_map : :class:`dict`
            A mapping from :class:`.FunctionalGroup` instances in
            `building_block` to the equivalent
            :class:`.FunctionalGroup` instances in the molecule being
            constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        fg1, fg2 = sorted(
            building_block.func_groups,
            key=lambda fg: building_block.get_centroid(
                atom_ids=fg.get_bonder_ids()
            )[0]
        )
        self.edges[0].assign_func_group(fg_map[fg1])
        self.edges[1].assign_func_group(fg_map[fg2])

    def __repr__(self):
        x, y, z = self._coord
        cls_name = (
            f'{__name__}.{self.__class__.__name__}'
        )
        # Make sure that the name has all the topology_graph submodule
        # names.
        p = re.compile(r'.*?topology_graphs\.(.*)', re.DOTALL)
        cls_name = p.findall(cls_name)[0]
        return (
            f'{cls_name}({x}, {y}, {z}, '
            f'orientation={self._orientation})'
        )


class _TerminalVertex(_LinearVertex):
    """
    Represents a vertex at the end of a polymer chain.

    Do not instantiate this class directly, use :class:`.HeadVertex` or
    :class:`.TailVertex` instead.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        `building_block` is placed such that the centroid-centroid
        direction vector is always pointing away from the center of the
        chain, when `building_block` has only one
        :class:`.FunctionalGroup`.

        If the `building_block` has more than one
        :class:`.FunctionalGroup` then this method behaves in the same
        was as for :class:`.LinearVertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed on the :class:`.Vertex`.
        """

        if len(building_block.func_groups) != 1:
            return super().place_building_block(building_block)

        building_block.set_centroid(
            position=self._coord,
            atom_ids=building_block.get_bonder_ids(fg_ids=(0, )),
        )
        centroid_vector = next(
            building_block.get_centroid_centroid_direction_vector()
        )
        building_block.apply_rotation_between_vectors(
            start=centroid_vector,
            # _cap_direction is defined by a subclass.
            target=[self._cap_direction, 0, 0],
            origin=self._coord
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        """
        Assign functional groups to edges.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

        fg_map : :class:`dict`
            A mapping from :class:`.FunctionalGroup` instances in
            `building_block` to the equivalent
            :class:`.FunctionalGroup` instances in the molecule being
            constructed.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`ValueError`
            If `building_block` does not have one or two functional
            groups.

        """

        if len(building_block.func_groups) == 2:
            fgs = sorted(
                building_block.func_groups,
                key=lambda fg: building_block.get_centroid(
                    atom_ids=fg.get_bonder_ids()
                )[0]
            )
            fg_index = 0 if self._cap_direction == 1 else -1
            fg = fgs[fg_index]

        elif len(building_block.func_groups) == 1:
            fg = building_block.func_groups[0]

        else:
            raise ValueError(
                'The building block of a polymer '
                'must have 1 or 2 functional groups.'
            )

        self.edges[0].assign_func_group(fg_map[fg])


class _HeadVertex(_TerminalVertex):
    """
    Represents a vertex at the head of a polymer chain.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = -1


class _TailVertex(_TerminalVertex):
    """
    Represents a vertex at the tail of a polymer chain.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = 1


class Linear(TopologyGraph):
    """
    Represents linear polymer topology graphs.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string specifying the repeating unit of the polymer.
        For example, ``"AB"`` or ``"ABB"``. Letters are assigned to
        building block molecules in the order they are passed to
        :meth:`.ConstructedMolecule.__init__`.

    orientations : :class:`tuple` of :class:`float`
        For each character in the repeating unit, a value between ``0``
        and ``1`` (both inclusive) must be given in a :class:`list`. It
        indicates the probability that each monomer will have its
        orientation along the chain flipped. If ``0`` then the
        monomer is guaranteed to not flip. If ``1`` it is
        guaranteed to flip. This allows the user to create
        head-to-head or head-to-tail chains, as well as chain with
        a preference for head-to-head or head-to-tail if a number
        between ``0`` and ``1`` is chosen.

    n : :class:`int`
        The number of repeating units which are used to make the
        polymer.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    def __init__(self, repeating_unit, orientations, n, processes=1):
        """
        Initialize a :class:`Linear` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string specifying the repeating unit of the polymer.
            For example, ``'AB'`` or ``'ABB'``. Letters are assigned to
            building block molecules in the order they are passed to
            :meth:`.ConstructedMolecule.__init__`.

        orientations : :class:`tuple` of :class:`float`
            For each character in the repeating unit, a value
            between ``0`` and ``1`` (both inclusive) must be given in
            a :class:`list`. It indicates the probability that each
            monomer will have its orientation along the chain flipped.
            If ``0`` then the monomer is guaranteed to not flip. If
            ``1`` it is guaranteed to flip. This allows the user to
            create head-to-head or head-to-tail chains, as well as
            chain with a preference for head-to-head or head-to-tail if
            a number between ``0`` and ``1`` is chosen.

        n : :class:`int`
            The number of repeating units which are used to make the
            polymer.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self.repeating_unit = repeating_unit
        self.orientations = tuple(orientations)
        self.n = n

        head, *body, tail = orientations*n
        vertices = [_HeadVertex(0, 0, 0, head)]
        edges = []
        for i, orientation in enumerate(body, 1):
            v = _LinearVertex(
                x=i, y=0, z=0, orientation=orientation
            )
            vertices.append(v)
            edges.append(Edge(vertices[i-1], vertices[i]))

        vertices.append(_TailVertex(len(vertices), 0, 0, tail))
        edges.append(Edge(vertices[-2], vertices[-1]))

        super().__init__(tuple(vertices), tuple(edges), processes)

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Assignment is done by modifying
        :attr:`.ConstructedMolecule.building_block_vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance being
            constructed.

        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        polymer = self.repeating_unit*self.n
        bb_map = {
            letter: bb for letter, bb in zip(polymer, building_blocks)
        }
        for letter, vertex in zip(polymer, self.vertices):
            bb = bb_map[letter]
            mol.building_block_vertices[bb].append(vertex)

    def _get_scale(self, mol):
        """
        Get the scale used for the positions of :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        """

        return max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )

    def __repr__(self):
        return (
            f'polymer.Linear({self.repeating_unit!r}, '
            f'{self.orientations!r}, '
            f'{self.n!r})'
        )
