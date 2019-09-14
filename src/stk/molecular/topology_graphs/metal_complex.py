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
