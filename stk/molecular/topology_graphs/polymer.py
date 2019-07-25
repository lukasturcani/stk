"""
Defines :class:`.Polymer` topologies.

"""

import logging
import re
from collections import defaultdict

from .topology_graph import TopologyGraph, Vertex, Edge


logger = logging.getLogger(__name__)


class LinearVertex(Vertex):
    """
    Represents a vertex in the middle of the chain.

    Attributes
    ----------
    _direction : :class:`int`
        Can be ``1``or ``-1`` to signify if the building block placed
        on the vertex should be placed parallel or anti-parallel
        to the chain.

    """

    def __init__(self, x, y, z, direction):
        self._direction = direction
        super().__init__(x, y, z)

    def clone(self):
        """
        Create a clone of the instance.

        Returns
        -------
        :class:`LinearVertex`
            A clone with the same position and direction and connected
            to the same :class:`.Edge` objects.

        """

        clone = super().clone()
        clone._direction = self._direction
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
        *_, bonder_vector = next(
            building_block.get_bonder_direction_vectors(
                fg_ids=(0, 1)
            )
        )
        building_block.apply_rotation_between_vectors(
            start=bonder_vector,
            target=[self._direction, 0, 0],
            origin=self._coord
        )
        return building_block.get_position_matrix()

    def _assign_func_groups_to_edges(self, building_block):
        fg1, fg2 = sorted(
            building_block.func_groups,
            key=lambda fg: building_block.get_centroid(
                atom_ids=fg.get_bonder_ids()
            )[0]
        )
        self.edges[0].assign_func_group(fg1)
        self.edges[1].assign_func_group(fg2)

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
            f'{cls_name}({x}, {y}, {z}, direction={self._direction})'
        )


class TerminalVertex(LinearVertex):
    """
    Represents a :class:`.Vertex` on the end of a polymer chain.

    Do not instantiate this class directly, use :class:`.HeadVertex` or
    :class:`.TailVertex` instead.

    Attributes
    ----------
    _cap_direction : :class:`int`
        The direction to use if the building block placed on the
        vertex only has 1 :class:`.FunctionalGroup`.

    """

    def __init__(self, x, y, z, direction):
        super().__init__(x, y, z, direction)

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
            target=[self._cap_direction, 0, 0],
            origin=self._coord
        )
        return building_block.get_position_matrix()

    def _assign_func_groups_to_edges(self, building_block):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be assigned to one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

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
            fg2, fg1 = sorted(
                building_block.func_groups,
                key=lambda fg: building_block.get_centroid(
                    atom_ids=fg.get_bonder_ids()
                )[0]
            )
        elif len(building_block.func_groups) == 1:
            fg1 = building_block.func_groups[0]
        else:
            raise ValueError(
                'The building block of a polymer '
                'must have 1 or 2 functional groups.'
            )

        self.edges[0].assign_func_group(fg1)


class HeadVertex(TerminalVertex):
    """
    Represents a vertex at the head of the chain.

    """

    _cap_direction = -1


class TailVertex(TerminalVertex):
    """
    Represents a vertex at the tail of the chain.

    """

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

    orientation : :class:`tuple` of :class:`float`
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

    ends : :class:`str`
        The string represents how the end groups of the polymer are
        treated. If ``'h'`` the functional groups at the end of the
        polymer are converted into hydrogem atoms. If ``'fg'`` they are
        kept as the original functional group.

    """

    def __init__(
        self,
        repeating_unit,
        orientation,
        n,
        ends='fg',
        processes=1
    ):
        """
        Initialize a :class:`Linear` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string specifying the repeating unit of the polymer.
            For example, ``'AB'`` or ``'ABB'``. Letters are assigned to
            building block molecules in the order they are passed to
            :meth:`.ConstructedMolecule.__init__`.

        orientation : :class:`tuple` of :class:`float`
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

        ends : :class:`str`, optional
            The string represents how the end groups of the polymer are
            treated. If ``'h'`` the functional groups at the end of the
            polymer are converted into hydrogem atoms. If ``'fg'`` they
            are kept as the original functional group.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self.repeating_unit = repeating_unit
        self.orientation = tuple(orientation)
        self.n = n
        self.ends = ends

        head, *body, tail = orientation*n
        vertices = [HeadVertex(0, 0, 0, head)]
        edges = []
        for i, direction in enumerate(body, 1):
            v = LinearVertex(
                x=i, y=0, z=0, direction=direction
            )
            vertices.append(v)
            edges.append(Edge(vertices[i-1], vertices[i]))

        vertices.append(TailVertex(len(vertices), 0, 0, tail))
        edges.append(Edge(vertices[-2], vertices[-1]))

        super().__init__(vertices, edges, processes)

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Note
        ----
        This method will modify
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

        mol.building_block_vertices = defaultdict(list)
        polymer = self.repeating_unit*self.n
        bb_map = {
            letter: bb for letter, bb in zip(polymer, building_blocks)
        }
        for letter, vertex in zip(polymer, self.vertices):
            bb = bb_map[letter]
            mol.building_block_vertices[bb].append(vertex)

    def _clean_up(self, mol):
        """
        Delete the atoms which are lost during construction.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.ends == 'h':
            self._hygrogen_ends(mol)

        super()._clean_up(mol)

    def _hygrogen_ends(self, mol):
        """
        Remove all deleter atoms and add hydrogens.

        In polymers, you may want to replace the functional groups at
        the ends with hydrogen atoms.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError('TODO')

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
