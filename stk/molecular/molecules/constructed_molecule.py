"""
Defines :class:`ConstructedMolecule`.

"""

import logging
import rdkit.Chem.AllChem as rdkit
from collections import Counter

from . import elements
from .molecule import Molecule
from .. import topology_graphs
from ..functional_groups import FunctionalGroup
from ...utilities import remake

logger = logging.getLogger(__name__)


class ConstructionError(Exception):
    ...


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    A :class:`ConstructedMolecule` requires at least 2 basic pieces of
    information: which building block molecules are used to construct
    the molecule and what the :class:`.TopologyGraph` of the
    constructe molecule is. The construction of the molecular structure
    is performed by :meth:`.TopologyGraph.construct`. This method does
    not have to be called explicitly by the user, it will be called
    automatically during initialization.

    The building block molecules used for construction can be either
    :class:`.BuildingBlock` instances or other
    :class:`.ConstructedMolecule` instances, or a combination both.

    Each :class:`.TopologyGraph` subclass may add additional attributes
    to the :class:`ConstructedMolecule`, which will be described within
    its documentation.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        Extends :class:`.Molecule.atoms`. Each :class:`.Atom`
        instance is guaranteed to have two additional attributes. The
        first is :attr:`building_block`, which holds the building
        block :class:`.Molecule` from which that
        :class:`.Atom` came. If the :class:`.Atom` did not come from a
        building block, but was added during contruction, the value
        of this attribute will be ``None``.

        The second attribute is :attr:`building_block_id`.
        Every time a building block is added to the
        :class:`ConstructedMolecule`, all added atoms will have the
        same :attr:`building_block_id`. This means that if you use the
        same building block twice during construction, each
        :class:`.Atom` in the :class:`.ConstructedMolecule`
        will have a :attr:`building_block_id` of either
        ``0``or ``1``.

    building_block_vertices : :class:`dict`
        Maps the :class:`.Molecule` instances used for construction,
        which can be either :class:`.BuildingBlock` or
        :class:`.ConstructedMolecule`, to the
        :class:`~.topologies.base.Vertex` objects they are placed on
        during construction. The :class:`dict` has the form

        .. code-block:: python

            building_block_vertices = {
                BuildingBlock(...): [Vertex(...), Vertex(...)],
                BuildingBlock(...): [
                    Vertex(...),
                    Vertex(...),
                    Vertex(...),
                ]
                ConstructedMolecule(...): [Vertex(...)]
            }

        Each :class:`.BuildingBlock` and :class:`ConstructedMolecule`
        can be mapped to multiple :class:`~.topologies.base.Vertex`
        objects.

    building_block_counter : :class:`collections.Counter`
        A counter keeping track of how many times each building block
        molecule appears in the :class:`ConstructedMolecule`.

    topology_graph : :class:`.TopologyGraph`
        Defines the topology graph of :class:`ConstructedMolecule` and
        is responsible for constructing it.

    bonds_made : :class:`int`
        The net number of bonds added during construction.

    func_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The remnants of building block functional groups present in the
        molecule. They track which atoms belonged to functional groups
        in the building block molecules. The id of each
        :class:`.FunctionalGroup` should match its index in
        :attr:`func_groups`.

    Methods
    -------
    :meth:`__init__`
    :meth:`to_dict`

    Examples
    --------

    """

    def __init__(
        self,
        building_blocks,
        topology_graph,
        building_block_vertices=None,
        use_cache=False
    ):
        """
        Initialize a :class:`ConstructedMolecule` instance.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        topolog_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`, optional
            Maps the :class:`.Molecule` in  `building_blocks` to the
            :class:`~.topologies.base.Vertex` in `topology_graph`.
            Each :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` can be mapped to multiple
            :class:`~.topologies.base.Vertex` objects. See the
            examples section in the :class:`.ConstructedMolecule`
            class docstring to help understand how this parameter
            is used. If ``None``, building block molecules will be
            assigned to vertices at random.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.ConstructedMolecule` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`ConstructedMolecule` does
            not yet exist the created one will be added to the cache.

        """

        # If building_block_vertices is "None" then
        # topology_graph.construct() is responsible for creating the
        # dict.
        self.building_block_vertices = building_block_vertices
        self.topology_graph = topology_graph
        self.atoms = []
        self.bonds = []
        self.bonds_made = 0
        self.func_groups = []
        self.building_block_counter = Counter()
        self._position_matrix = []

        try:
            topology_graph.construct(self, building_blocks)

        except Exception as ex:
            errormsg = (
                'Construction failure.\n'
                '\n'
                'topology_graph\n'
                '--------------\n'
                f'{topology_graph}\n'
                '\n'
                'building blocks\n'
                '---------------\n'
            )

            bb_blocks = []
            for i, bb in enumerate(building_blocks):
                bb_blocks.append(
                    f'{bb}\n\n'
                    'MDL MOL BLOCK\n'
                    '-------------\n'
                    f'{bb._to_mdl_mol_block()}'
                )

            errormsg += '\n'.join(bb_blocks)
            raise ConstructionError(errormsg) from ex

        self.atoms = tuple(self.atoms)
        self.bonds = tuple(self.bonds)
        self.func_groups = tuple(self.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(self.func_groups):
            func_group.id = id_

    def to_dict(self, include_attrs=None, ignore_missing_attrs=False):
        """
        Return a :class:`dict` representation of the molecule.

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of additional attributes of the molecule to be
            added to the :class:`dict`. Each attribute is saved as a
            string using :func:`repr`. These attributes are also
            passed down recursively to the building block molecules.

        ignore_missing_attrs : :class:`bool`, optional
            If ``False`` and an attribute in `include_attrs` is not
            held by the :class:`ConstructedMolecule`, an error will be
            raised.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        bb_counter = []
        building_blocks = []
        building_block_vertices = []
        for i, bb in enumerate(self.building_block_vertices):
            bb_counter.append((i, self.bb_counter[bb]))
            building_blocks.append(bb.to_dict(include_attrs, True))
            building_block_vertices[i] = [
                repr(v) for v in self.building_block_vertices[bb]
            ]

        d = {
            'building_blocks': building_blocks,
            'building_block_vertices': building_block_vertices,
            'building_block_counter': bb_counter,
            'bonds_made': self.bonds_made,
            'class': self.__class__.__name__,
            'mol_block': self._to_mdl_mol_block(),
            'topology_graph': repr(self.topology_graph),
            'func_groups': repr(self.func_groups),
            'atoms': repr(self.atoms),
        }

        d.update(
            {attr: repr(getattr(self, attr)) for attr in include_attrs}
        )
        return d

    @classmethod
    def _init_from_dict(cls, mol_dict, use_cache):
        """
        Intialize from a :class:`dict` representation.

        Parameters
        ----------
        mol_dict : :class:`dict`
            A :class:`dict` representation of a molecule generated
            by :meth:`to_dict`.

        use_cache : :class:`bool`
            If ``True``, a new instance will not be made if a cached
            and identical one already exists, the one which already
            exists will be returned. If ``True`` and a cached,
            identical instance does not yet exist the created one will
            be added to the cache.

        Returns
        -------
        :class:`ConstructedMolecule`
            The molecule described by `mol_dict`.

        """

        d = dict(mol_dict)
        tops = vars(topology_graphs)
        topology_graph = eval(d.pop('topology_graph'), tops)

        bbs = [
            Molecule.init_from_dict(bb_dict, use_cache)
            for bb_dict in d.pop('building_blocks')
        ]
        key = cls._get_key(cls, bbs, topology_graph, use_cache)
        if key in cls._cache and use_cache:
            return cls.cache[key]

        obj = cls.__new__(cls)

        mol = remake(rdkit.MolFromMolBlock(
            molBlock=d.pop('mol_block'),
            sanitize=False,
            removeHs=False
        ))

        obj.building_block_counter = Counter()
        obj.building_block_vertices = {}

        counter = d.pop('bb_counter')
        vertices = d.pop('building_block_vertices')
        for i, bb in enumerate(bbs):
            obj.building_block_counter[bb] = counter[i]
            obj.building_block_vertices[bb] = [
                eval(v, tops) for v in vertices[i]
            ]

        obj.topology_graph = topology_graph
        obj.bonds_made = d.pop('bonds_made')
        obj._key = key
        obj._position_matrix = mol.GetConformer().GetPositions().T
        obj.atoms = eval(d.pop('atoms'), vars(elements))

        obj.bonds = tuple(
            elements.Bond(
                obj.atoms[b.GetBeginAtomIdx()],
                obj.atoms[b.GetEndAtomIdx()],
                b.GetBondTypeAsDouble()
            )
            for b in mol.GetBonds()
        )

        # Globals for eval.
        g = {'FunctionalGroup': FunctionalGroup}
        obj.func_groups = tuple(eval(d.pop('func_groups'), g))
        if use_cache:
            cls._cache[key] = obj

        d.pop('class')
        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        return obj

    @staticmethod
    def _get_key(
        self,
        building_blocks,
        topology_graph,
        building_block_vertices=None,
        use_cache=False
    ):
        """
        Get the key used for caching the molecule.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        topolog_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`, optional
            Maps the :class:`.Molecule` in  `building_blocks` to the
            :class:`~.topologies.base.Vertex` in `topology_graph`.
            Each :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` can be mapped to multiple
            :class:`~.topologies.base.Vertex` objects. See the
            examples section in the :class:`.ConstructedMolecule`
            class docstring to help understand how this parameter
            is used. If ``None``, building block molecules will be
            assigned to vertices at random.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.ConstructedMolecule` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`ConstructedMolecule` does
            not yet exist the created one will be added to the cache.

        """

        bb_keys = frozenset(x._key for x in building_blocks)
        return bb_keys, repr(topology_graph)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            '(building_blocks='
            f'{[str(x) for x in self.building_blocks]}, '
            f'topology_graph={self.topology_graph!r})'
        )

    def __repr__(self):
        return str(self)
