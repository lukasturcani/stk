"""
Metal Complex
=============

#. :class:`.SquarePlanar`

"""

import logging

from .topology_graph import TopologyGraph, Vertex, Edge


logger = logging.getLogger(__name__)


class MetalComplex(TopologyGraph):
    """
    Represents single-molecule metal complex topology graphs.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------


    """

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
            for bb in mol.building_block_vertices[1:]
        )


class SquarePlanar(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    See :class:`.MetalComplex` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertices = (
        Vertex(0, 0, 0),
        Vertex(0, 1, 0),
        Vertex(0, 0, 1),
        Vertex(0, -1, 0),
        Vertex(0, 0, -1),
    )

    edges = (
        Edge(vertices[0], vertices[1]),
        Edge(vertices[0], vertices[2]),
        Edge(vertices[0], vertices[3]),
        Edge(vertices[0], vertices[4]),
    )
