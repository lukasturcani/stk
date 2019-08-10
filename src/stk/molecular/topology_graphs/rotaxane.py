"""
Defines rotaxane topologies.

"""


import numpy as np
import rdkit.Chem.AllChem as rdkit

from .topology_graph import TopologyGraph, Vertex


class _AxleVertex(Vertex):
    ...


class _CycleVertex(Vertex):
    ...


class NRotaxane(TopologyGraph):
    """
    Represents [n]rotaxane topology graphs.

    This class assumes one axle with (n-1) macrocycles threaded on it.
    The macrocycles are spaced evenly along the thread in repeating
    patterns. The threaded macrocycles can be described analagously
    to monomers in linear polymers, in terms of a repeating unit,
    except that no bonds are formed between them.

    The axle must be provided first to the `building_blocks` in
    :class:`.ConstructedMolecule.__init__`.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------
    .. code-block:: python

        import stk

    """

    def __init__(self, repeating_unit, orientation, n):
        """
        Initialize a :class:`NRotaxane` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string specifying the repeating unit of the macrocycles.
            For example, ``'AB'`` or ``'ABB'``. Letters are assigned to
            building block molecules in the order they are passed to
            :meth:`.ConstructedMolecule.__init__`.

        orientations : :class:`tuple` of :class:`float`
            For each character in the repeating unit, a value
            between ``0`` and ``1`` (both inclusive) must be given in
            a :class:`tuple`. It indicates the probability that each
            macrocycle will have its orientation along the axle
            flipped. If ``0`` then the macrocycle is guaranteed not to
            flip. If ``1`` it is guaranteed to flip. This allows the
            user to create head-to-head or head-to-tail chains, as well
            as chain with a preference for head-to-head or head-to-tail
            if a number between ``0`` and ``1`` is chosen.

        n : :class:`int`
            The number of repeating units threaded along the axle.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.


        """

        self._repeating_unit = repeating_unit
        self._orientations = orientation
        self._n = n
