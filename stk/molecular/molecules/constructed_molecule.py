"""
Defines :class:`ConstructedMolecule`.

.. _`molecular construction`:

A more detailed description of molecular construction.
------------------------------------------------------

This is a step-by-step guide of how :class:`.ConstructedMolecule`
instances are constructed.

First, you create :class:`.BuildingBlock` instances of the building blocks
which make up the :class:`ConstructedMolecule`:

.. code-block:: python

    bb = BuildingBlock('/path/to/struct/file.mol2', ['amine'])

The :class:`.BuildingBlock` instances are initialized using paths to
molecular structure files or :mod:`rdkit` molecules or with SMILES
strings. Initializing a :class:`.BuildingBlock` automatically completes
steps 1 to 4.

    1. Place an :mod:`rdkit` instance of the molecule into
       :attr:`.BuildingBlock.mol`, i.e.

       .. code-block:: python

           bb.mol  # <rdkit.Chem.rdchem.Mol at 0x7f961a8f1f80>

    2. Scan the path of the structure file for the names of functional
       groups. (Alternatively the names of functional groups can be
       supplied to the initializer). Find the :class:`.FGInfo` instance
       for each functional group.


Which functional groups are recognized by ``stk``?

The module :mod:`.functional_groups` defines the class :class:`.FGInfo`
and a :class:`tuple` of instances of this class called
:data:`functional_groups`. If you put an :class:`.FGInfo` instance into
:data:`functional_groups`, the functional group will be recognized.

    3. Using :class:`.FGInfo` create :class:`.FunctionalGroup`
       instances, which determine the bonder and deleter atoms in the
       molecule. These identify which atoms form bonds during
       the construction of a :class:`ConstructedMolecule` and which
       ones are deleted. Place the :class:`.FunctionalGroup` instances
       into :attr:`.BuildingBlock.func_groups`.

       .. code-block:: python

           bb.func_groups
           # (
           #     FunctionalGroup(
           #         id_=0,
           #         atom_ids=(45, 21, 0),
           #         bonder_ids=(21, ),
           #         deleter_ids=(0, 45),
           #         info=FGInfo('amine')
           #     ),
           #     FunctionalGroup(
           #         id_=1,
           #         atom_ids=(47, 23, 15),
           #         bonder_ids=(47, ),
           #         deleter_ids=(23, 15),
           #         info=FGInfo('amine')
           #     )
           # )

    5. Initialize an instance of :class:`.ConstructedMolecule`.

       .. code-block:: python

           mol = ConstructedMolecule([bb1, bb2], Topology())

       Normally, :class:`.ConstructedMolecule` and :class:`.Topology`
       will not be used directly. Instead, classes derived from these
       will be used. For example,

           .. code-block:: python

               polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

    6. Run :meth:`.Topology.construct` inside
       :meth:`ConstructedMolecule.__init__`.

    7. The details of :meth:`.Topology.construct` will vary depending
       on the :class:`.Topology` class used. However, the basic
       structure is the same (steps 8 - 10).

    8. Use :meth:`.Topology.place_mols` to combine the :mod:`rdkit`
       molecules of all building blocks into a single :mod:`rdkit`
       instance. The combined :mod:`rdkit` instance is placed into the
       :attr:`.ConstructedMolecule.mol`. attribute.
       :meth:`.Topology.place_mols` also usually
       keeps track of each functional groups in the constructed
       molecule. If a buliding block is placed during the construction
       of a molecule, the atom ids have to be shifted upward by some
       amount. :meth:`.FunctionalGroup.shifted_fg` performs this
       operation.

    9. Use :meth:`.Topology.prepare` to run any additional operations
       before joining up the building blocks and deleting extra
       atoms, this method may do nothing.

    10. Use :meth:`.Topology.bonded_fgs` to yield the functional groups
        which react. The :class:`.FunctionalGroup` instances in the
        constructed molecule are passed to :meth:`.Reactor.react`,
        which performs the reaction. See the documentation of
        :class:`.Reactor` for information on how reactions are carried
        out.

    11. Run :meth:`.Topology.cleanup` to perform any final operations
        on the constructed molecule. Can be nothing.

After all this you should have a :mod:`rdkit` instance of the
constructed molecule, which should be placed into
:attr:`ConstructedMolecule.mol`.

.. _`adding constructed molecules`:

Extending stk: Adding new types of constructed molecules.
---------------------------------------------------------

To add new constructed molecules create a new class which inherits
:class:`.ConstructedMolecule`.

If you're adding a new class of constructed molecules, it quite likely
you want to add a new :class:`.Topology` class. See the
:mod:`.topologies.base` for guidance on adding these. The topology
class does the construction of the molecule from the building blocks.

"""

import logging
import rdkit.Chem.AllChem as rdkit
from collections import Counter

from .molecule import Molecule
from .. import topologies
from ..functional_groups import FunctionalGroup

logger = logging.getLogger(__name__)


class ConstructionError(Exception):
    ...


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    A :class:`ConstructedMolecule` requires at least 2 basic pieces of
    information: which building block molecules are used to construct
    the molecule and what the :class:`.Topology` of the constructed
    molecule is. The construction of the molecular structure is
    performed by :meth:`.Topology.construct`. This method does not
    have to be called explicitly by the user, it will be called
    automatically during initialization.

    Each :class:`.Topology` may add additional attributes to the
    :class:`ConstructedMolecule`, which will be documented by that
    class.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        Extends :class:`.Molecule.atoms`. Each :class:`.Atom`
        instance has two additional attributes. The
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

    building_blocks : :class:`dict`
        Maps the building blocks used for construction, which can be
        either :class:`.BuildingBlock` or
        :class:`.ConstructedMolecule`, to the
        :class:`~.topologies.base.Vertex` objects they are placed on
        during construction. The :class:`dict` has the form

        .. code-block:: python

            building_blocks = {
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

    building_block_conformers : :class:`list` of :class:`dict`
        

    bb_counter : :class:`collections.Counter`
        A counter keeping track of how many times each building block
        molecule appears in the :class:`ConstructedMolecule`.

    topology : :class:`.Topology`
        Defines the topology of :class:`ConstructedMolecule` and
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
    :meth:`add_conformer`

    Examples
    --------

    """

    def __init__(
        self,
        building_blocks,
        topology,
        bb_conformers=None,
        use_cache=False
    ):
        """
        Initialize a :class:`ConstructedMolecule` instance.

        Parameters
        ---------
        building_blocks : :class:`list` of :class:`.BuildingBlock`
            The :class:`.BuildingBlock` instances which
            represent the building block molecules of the
            :class:`ConstructedMolecule`. Only one
            :class:`.BuildingBlock` instance is present per building
            block, even if multiples of that building block join up to
            form the :class:`ConstructedMolecule`.

        topology : :class:`.Topology`
            Defines the topology of the :class:`ConstructedMolecule`
            and constructs it.

        bb_conformers : :class:`list` of :class:`int`, optional
            The ids of the building block conformers to be used. Must
            be equal in length to `building_blocks` and orders must
            correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.ConstructedMolecule` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`ConstructedMolecule` does
            not yet exist the created one will be added to the cache.

        """

        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(building_blocks))]

        self.building_blocks = building_blocks
        self.topology = topology

        try:
            # Ask the ``Topology`` instance to construct the
            # molecule. This creates the `mol`, `bonds_made` and
            # `func_groups` attributes.
            topology.construct(self, bb_conformers)

        except Exception as ex:
            errormsg = (
                'Construction failure.\n'
                '\n'
                'topology\n'
                '--------\n'
                f'{topology}\n'
                '\n'
                'building blocks\n'
                '---------------\n'
            )

            bb_blocks = []
            for i, bb in enumerate(building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append(
                    f'{bb.__class__.__name__} '
                    f'{[info.name for info in bb.func_group_infos]}\n'
                    f'{bb.mdl_mol_block(bb_conf)}'
                )

            errormsg += '\n'.join(bb_blocks)
            raise ConstructionError(errormsg) from ex

        self.func_groups = tuple(self.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(self.func_groups):
            func_group.id = id_

        super().__init__()

    def add_conformer(self, bb_conformers):
        """
        Construct a new conformer.

        Parameters
        ----------
        bb_conformers : :class:`list` of :class:`int`
            The ids of the building block conformers to be used. Must
            be equal in length to :attr:`building_blocks` and the
            orders must correspond. If ``None``, then ``-1`` is used
            for all building blocks.

        Returns
        -------
        :class:`int`
            The id of the new conformer.

        """

        # Save the original rdkit molecule.
        original_mol = self._mol
        # Construct a new molecule.
        try:
            # Ask the ``Topology`` instance to construct the
            # molecule. This creates the `mol`, `bonds_made`
            # and `func_groups` attributes.
            self.topology.construct(self, bb_conformers)

        except Exception as ex:
            self._mol = original_mol
            errormsg = (
                'Construction failure.\n'
                '\n'
                'topology\n'
                '--------\n'
                f'{self.topology}\n'
                '\n'
                'building blocks\n'
                '---------------\n'
            )

            bb_blocks = []
            for i, bb in enumerate(self.building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append(
                    f'{bb.__class__.__name__} '
                    f'{[info.name for info in bb.func_group_infos]}\n'
                    f'{bb.mdl_mol_block(bb_conf)}'
                )

            errormsg += '\n'.join(bb_blocks)
            raise ConstructionError(errormsg) from ex

        self.func_groups = tuple(self.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(self.func_groups):
            func_group.id = id_

        # Get the new conformer.
        new_conf = rdkit.Conformer(self._mol.GetConformer())
        # Add it to the original molecule.
        new_id = original_mol.AddConformer(new_conf, True)
        self._mol = original_mol
        return new_id

    def _to_dict(self, include_attrs=None):
        """
        Returns a :class:`dict` representation of the molecule.

        The representation has the form

        .. code-block:: python

            {
                'class' : 'Polymer',
                'mol_block' : '''A string holding the V3000 mol
                                 block of the molecule.''',
                'building_blocks' : {bb1._to_dict(), bb2._to_dict()},
                'topology' : 'Copolymer(repeating_unit="AB")',
            }

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecule to be added to
            the :class:`dict`. Each attribute is saved as a string
            using :func:`repr`.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        conformers = [
            (
                conf.GetId(),
                self.to_mdl_mol_block(conformer=conf.GetId())
            )
            for conf in self._mol.GetConformers()
        ]

        d = {
            'bb_counter': [
                (key._to_dict(), val)
                for key, val in self.bb_counter.items()
            ],
            'bonds_made': self.bonds_made,
            'class': self.__class__.__name__,
            'conformers': conformers,
            'building_blocks': [
                x._to_dict() for x in self.building_blocks
            ],
            'topology': repr(self.topology),
            'func_groups': repr(self.func_groups)

        }

        d.update(
            {attr: repr(getattr(self, attr)) for attr in include_attrs}
        )
        return d

    @classmethod
    def _init_from_dict(cls, mol_dict, use_cache):
        """
        Initialize from a :class:`dict` representation.

        Parameters
        ----------
        mol_dict : :class:`dict`
            A dictionary holding the attribute data of the molecule.

        use_cache : :class:`bool`
            If ``True``, a new :class:`.ConstructedMolecule` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`ConstructedMolecule` does
            not yet exist the created one will be added to the cache.

        Returns
        -------
        :class:`ConstructedMolecule`
            The molecule.

        """

        d = dict(mol_dict)
        d.pop('building_blocks')
        d.pop('class')

        bb_counter = Counter({
            Molecule.from_dict(key): val
            for key, val in d.pop('bb_counter')
        })
        bbs = list(bb_counter)
        topology = eval(d.pop('topology'),  topologies.__dict__)

        key = cls._generate_key(bbs, topology)
        if key in cls._cache and use_cache:
            return cls.cache[key]

        obj = cls.__new__(cls)

        (conf_id, mol_block), *confs = d.pop('conformers')
        obj._mol = rdkit.MolFromMolBlock(
            molBlock=mol_block,
            sanitize=False,
            removeHs=False
        )
        obj._mol.GetConformer().SetId(conf_id)

        for conf_id, mol_block in confs:
            conf_mol = rdkit.MolFromMolBlock(
                molBlock=mol_block,
                sanitize=False,
                removeHs=False
            )
            conf = conf_mol.GetConformer()
            conf.SetId(conf_id)
            obj._mol.AddConformer(conf)

        obj.topology = topology
        obj.bb_counter = bb_counter
        obj.bonds_made = d.pop('bonds_made')
        obj._key = key
        obj.building_blocks = bbs

        # Globals for eval.
        g = {'FunctionalGroup': FunctionalGroup}
        obj.func_groups = tuple(eval(d.pop('func_groups'), g))
        if use_cache:
            cls._cache[key] = obj

        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        return obj

    @classmethod
    def _generate_key(
        cls,
        building_blocks,
        topology,
        bb_conformers,
        use_cache
    ):
        """
        Generates the key used for caching the molecule.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.BuildingBlock`
            The :class:`.BuildingBlock` instances which
            represent the building block molecules of the
            :class:`ConstructedMolecule`. Only one
            :class:`.BuildingBlock` instance is present per building
            block, even if multiples of that building block join up to
            form the :class:`ConstructedMolecule`.

        topology : :class:`.Topology`
            Defines the topology of the :class:`ConstructedMolecule`
            and constructs it.

        bb_conformers : :class:`list` of :class:`int`
            The ids of the building block conformers to be used. Must
            be equal in length to `building_blocks` and orders must
            correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        use_cache : :class:`bool`
            This argument is ignored but included to be maintain
            compatiblity the the :meth:`__init__` signature.

        """

        bb_keys = frozenset(x.key for x in building_blocks)
        return bb_keys, repr(topology)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            '(building_blocks='
            f'{[str(x) for x in self.building_blocks]}, '
            f'topology={self.topology!r})'
        )

    def __repr__(self):
        return str(self)
