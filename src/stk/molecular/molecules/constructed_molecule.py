"""
Constructed Molecule
====================

"""

import logging
import rdkit.Chem.AllChem as rdkit

from .molecule import Molecule
from ..atoms import AtomInfo
from ..bonds import BondInfo

logger = logging.getLogger(__name__)


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    Examples
    --------
    *Initialization*

    A :class:`.ConstructedMolecule` is initialized from a
    :class:`.TopologyGraph`, which is typically initialized from some
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

    A :class:`ConstructedMolecule` may be used to construct other
    :class:`ConstructedMolecule` instances, though you will probably
    have to convert it to a :class:`.BuildingBlock` first

    .. code-block:: python

        benzene = stk.BuildingBlock('c1ccccc1')
        cage_complex = stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(cage1),
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
        self._num_building_blocks = {
            building_block:
                topology_graph.get_num_building_block(building_block)
            for building_block in topology_graph.get_building_blocks()
        }

    @classmethod
    def init(
        cls,
        atoms,
        bonds,
        position_matrix,
        atom_infos,
        bond_infos,
        num_building_blocks,
    ):
        """
        Initialize a :class:`.ConstructedMolecule` from its components.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms of the molecule.

        bond : :class:`tuple` of :class:`.Bond`
            The bonds of the molecule.

        position_matrix : :class:`numpy.ndarray`
            A ``(n, 3)`` position matrix of the molecule.

        atom_infos : :class:`tuple` of :class:`.AtomInfo`
            The atom infos of the molecule.

        bond_infos : :class:`tuple` of :class:`.BondInfo`
            The bond infos of the molecule.

        num_building_blocks : :class:`dict`
            Maps each building block of the constructed molecule to
            the number of times it is present in it.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The constructed molecule.

        """

        molecule = cls.__new__(cls)
        Molecule.__init__(molecule, atoms, bonds, position_matrix)
        molecule._atom_infos = atom_infos
        molecule._bond_infos = bond_infos
        molecule._num_building_blocks = dict(num_building_blocks)
        return molecule

    def clone(self):
        clone = super().clone()
        clone._atom_infos = self._atom_infos
        clone._bond_infos = self._bond_infos
        clone._num_building_blocks = dict(self._num_building_blocks)
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
        :class:`.Molecule`
            A building block of the constructed molecule.

        """

        yield from self._num_building_blocks

    def get_num_building_block(self, building_block):
        """
        Get the number of times `building_block` is present.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block whose frequency in the constructed
            molecule is desired.

        Returns
        -------
        :class:`int`
            The number of times `building_block` was used in the
            construction of the constructed molecule.

        """

        return self._num_building_blocks[building_block]

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

    def _with_canonical_atom_ordering(self):
        ordering = rdkit.CanonicalRankAtoms(self.to_rdkit_mol())
        id_map = {
            new_id: atom.get_id()
            for new_id, atom in zip(ordering, self._atoms)
        }
        super()._with_canonical_atom_ordering()
        atom_map = {
            old_id: self._atoms[new_id]
            for old_id, new_id in enumerate(ordering)
        }
        old_atom_infos = self._atom_infos

        def get_atom_info(atom):
            info = old_atom_infos[id_map[atom.get_id()]]
            return AtomInfo(
                atom=atom,
                building_block=info.get_building_block(),
                building_block_id=info.get_building_block_id(),
            )

        def get_bond_info(info):
            return BondInfo(
                bond=info.get_bond().with_atoms(atom_map),
                building_block=info.get_building_block(),
                building_block_id=info.get_building_block_id(),
            )

        self._atom_infos = tuple(map(get_atom_info, self._atoms))
        self._bond_infos = tuple(map(get_bond_info, self._bond_infos))
        return self
