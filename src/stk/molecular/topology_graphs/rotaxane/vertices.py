"""
Rotaxane Vertices
=================

"""

import rdkit.Chem.AllChem as rdkit

from ..topology_graph import Vertex


class AxleVertex(Vertex):
    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {}


class CycleVertex(Vertex):
    """
    Places the cycles in a :class:`NRotaxane`.

    """

    def __init__(self, id, position, flip):
        """
        Initialize a :class:`.CycleVertex` instance.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        flip : :class:`bool`
            If ``True``, the orientation of building blocks placed by
            the vertex will be flipped.

        """

        self._flip = flip
        super().__init__(id, position)

    def get_flip(self):
        """
        Return ``True`` if the vertex flips building blocks it places.

        Returns
        -------
        :class:`bool`
            ``True`` if the vertex flips building blocks it places.

        """

        return self._flip

    def clone(self):
        clone = super().clone()
        clone._flip = self._flip
        return clone

    def place_building_block(self, building_block, edges):
        rdkit_mol = building_block.to_rdkit_mol()
        macrocycle = max(rdkit.GetSymmSSSR(rdkit_mol), key=len)
        return building_block.with_centroid(
            position=self._position,
            atom_ids=macrocycle,
        ).with_rotation_between_vectors(
            start=building_block.get_plane_normal(macrocycle),
            target=[-1 if self._flip else 1, 0, 0],
            origin=self._position
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {}

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'flip={self._flip})'
        )
