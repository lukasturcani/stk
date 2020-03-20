"""
Constructed Molecule
====================

"""

import logging
import numpy as np

from .molecule import Molecule

logger = logging.getLogger(__name__)


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    Examples
    --------
    *Initialization*

    A :class:`.ConstructedMolecule` is initialized from a
    :class:`.TopologyGraph`, which is initialized from some
    :class:`.BuildingBlock` instances.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)CC=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        tetrahedron = stk.cage.FourPlusSix((bb1, bb2))
        cage1 = stk.ConstructedMolecule(tetrahedron)

    Depending on the :class:`.TopologyGraph`, a
    :class:`ConstructedMolecule` may be used to construct other
    :class:`ConstructedMolecule` instances.

    .. code-block:: python

        benzene = stk.BuildingBlock('c1ccccc1')
        cage_complex = stk.host_guest.Complex(
            host=cage1,
            guest=benzene,
        )
        cage_complex = stk.ConstructedMolecule(host_guest_complex)

    Obviously, the initialization of the :class:`.ConstructedMolecule`
    depends mostly on the specifics of the :class:`.TopologyGraph`
    used, and the documentation of those classes should be examined
    for more examples.

    """

    def __init__(self, topology_graph):
        """
        Initialize a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        topology_graph : :class:`.TopologyGraph`
            The topology graph of the constructed molecule.

        """

        construction_result = topology_graph.construct()
        super().__init__(
            atoms=construction_result.get_atoms(),
            bonds=construction_result.get_bonds(),
            position_matrix=construction_result.get_position_matrix(),
        )
        self._atom_infos = construction_result.get_atom_infos()
        self._bond_infos = construction_result.get_bond_infos()
        self._topology_graph = topology_graph

    def clone(self):
        clone = super().clone()
        clone._atom_infos = self._atom_infos
        clone._bond_infos = self._bond_infos
        clone._topology_graph = self._topology_graph
        return clone

    def get_building_blocks(self):
        """
        Yield the building blocks of the constructed molecule.

        Building blocks are yielded in an order based on their
        position in the constructed molecule. For two topologically
        equivalent constructed molecules, but with different building
        blocks, equivalently positioned building blocks will be
        yielded at the same time.

        Yields
        ------
        :class:`.BuildingBlock`
            A building block of the constructed molecule.

        """

        yield from self._topology_graph.get_building_blocks()

    def get_num_building_block(self, building_block):
        """
        Get the number of times `building_block` is present.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block whose frequency in the constructed
            molecule is desired.

        Returns
        -------
        :class:`int`
            The number of times `building_block` was used in the
            construction of the constructed molecule.

        """

        return (
            self._topology_graph.get_num_building_block(building_block)
        )

    def with_building_blocks(self, building_block_map):
        """
        Return a made of different building blocks.

        Parameters
        ----------
        building_block_map : :class:`dict`
            Maps a building block in the current constructed molecule
            to the building block which should replace it in the clone.
            If a building block should be not replaced in the clone, it
            can be omitted from the map.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The clone. Has the same type as the original constructed
            molecule.

        Examples
        --------
        .. code-block:: python

            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
            bb2 = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

            linear = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(bb1, bb2),
                    repeating_unit='A',
                    num_repeating_units=15,
                ),
            )

            bb3 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
            # All bb1 instances are replaced by bb3, but bb2 remains
            # in place.
            clone = linear.with_building_blocks({
                bb1: bb3,
            })

        """

        return self.clone()._with_building_blocks(building_block_map)

    def _with_building_blocks(self, building_block_map):
        """
        Modify the instance.

        """

        topology_graph = self._topology_graph.with_building_blocks(
            building_block_map=building_block_map,
        )
        self._topology_graph = topology_graph

        construction_result = topology_graph.construct()
        self._atoms = construction_result.get_atoms()
        self._bonds = construction_result.get_bonds()
        self._position_matrix = np.array(
            construction_result.get_position_matrix()
        )
        self._atom_infos = construction_result.get_atom_infos()
        self._bond_infos = construction_result.get_bond_infos()
        return self

    def get_atom_infos(self, atom_ids=None):
        """
        Yield data about atoms in the molecule.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms whose data is desired. If ``None``,
            data on all atoms will be yielded. Can be a single
            :class:`int`, if data on a single atom is desired.

        Yields
        ------
        :class:`.AtomInfo`
            Data about an atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        for atom_id in atom_ids:
            yield self._atom_infos[atom_id]

    def get_bond_infos(self):
        """
        Yield data about bonds in the molecule.

        Yields
        ------
        :class:`.BondInfo`
            Data about a bond.

        """

        yield from self._bond_infos

    def __str__(self):
        return f'{self.__class__.__name__}({self._topology_graph!r})'

    def __repr__(self):
        return str(self)
