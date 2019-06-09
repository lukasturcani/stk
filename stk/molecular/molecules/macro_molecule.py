"""
Defines :class:`MacroMolecule`.

The :class:`MacroMolecule` class represents an assembled macromolecule.
It requires at least 2 basic pieces of information: which monomers are
used to assemble the macromolecule and what the topology/structure
of the macromolecule is.

:attr:`MacroMolecule.building_blocks` holds a :class:`list` of
:class:`.StructUnit` instances. These represent the monomers which
make up the macromolecule. Only one :class:`.StructUnit` instance per
monomer type is held. So if 4 of one type of monomer and 2 of another
type of monomer form a macromolecule, only 2 :class:`.StructUnit`
instances are in :attr:`MacroMolecule.building_blocks`.

:attr:`MacroMolecule.topology` holds a :class:`.Topology` instance.
This instance is responsible for assembling the macromolecule from the
building blocks. The building should happen in
:meth:`MacroMolecule.__init__` via :meth:`.Topology.build`. The
:meth:`~.Topology.build` method places the assembled macromolecule in
:attr:`MacroMolecule.mol` as an ``rdkit`` molecule.

.. _`macromolecular assembly`:

A more detailed description of macromolecular assembly.
-------------------------------------------------------

This is a step-by-step guide of how macromolecular assembly is carried
out and what the classes do.

First you create :class:`.StructUnit` instances of the building blocks
which make up the macromolecule:

.. code-block:: python

    bb = StructUnit('/path/to/struct/file.mol2', ['amine'])

The :class:`.StructUnit` instances are initialized using paths to
molecular structure files. (Initializing a :class:`.StructUnit`
automatically completes steps 1 to 4.)

    1. Place an ``rdkit`` instance of the molecule into
       :attr:`StructUnit.mol`, i.e.

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
       macromolecular assembly and which ones are deleted. Place the
       :class:`.FunctionalGroup` instances into
       :attr:`StructUnit.func_groups`.

       .. code-block:: python

           bb.func_groups
           # (FunctionalGroup(id=0,
           #                  atom_ids=(45, 21, 0),
           #                  bonder_ids=(21, ),
           #                  deleter_ids=(0, 45),
           #                  info=FGInfo('amine')),
           #  FunctionalGroup(id=1,
           #                  atom_ids=(47, 23, 15),
           #                  bonder_ids=(47, ),
           #                  deleter_ids=(23, 15),
           #                  info=FGInfo('amine')))

    5. Give the :class:`.StructUnit` and :class:`.Topology` instances
       to the macromolecule's initializer.

       .. code-block:: python

           macro_mol = MacroMolecule([bb1, bb2], Topology())

       Normally, :class:`MacroMolecule` and :class:`.Topology` will
       not be used directly. Instead, classes derived from these
       will be used. For example,

           .. code-block:: python

               macro_mol = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

    6. Run :meth:`.Topology.build` inside
       :meth:`MacroMolecule.__init__`.

    7. The details of :meth:`.Topology.build` will vary depending on
       the :class:`.Topology` class used. However, the basic structure
       is the same (steps 8 - 10).

    8. Use :meth:`.Topology.place_mols` to combine the ``rdkit``
       molecules of all the building blocks into a single ``rdkit``
       instance. The combined ``rdkit`` instance is placed into
       ``macro_mol.mol``. :meth:`.Topology.place_mols` also usually
       keeps track of each functional group in the macromolecule.
       If a buliding block is placed in a macromolecule, the atom
       ids have to be shifted upward by some amount.
       :meth:`.FunctionalGroup.shifted_fg` performs this operation.

    9. Use :meth:`.Topology.prepare` to run any additional operations
       before joining up the building blocks and deleting extra
       atoms, this method may do nothing.

    10. Use :meth:`.Topology.bonded_fgs` to yield the functional groups
        which react. The :class:`.FunctionalGroup`s in the
        macromolecule are passed to :meth:`.Reactor.react`, which
        performs the reaction. See the documentation of
        :class:`Reactor` for information on how reactions are carried
        out.

    11. Run :meth:`.Topology.cleanup` to perform any final operations
        on the assembled molecule. Can be nothing.

After all this you should have a ``rdkit`` instance of the
macromolecule which should be placed into :attr:`MacroMolecule.mol`.

.. _`adding macromolecules`:

Extending stk: Adding new macromolecules.
-----------------------------------------

To add new macromolecules create a new class which inherits
:class:`MacroMolecule`.

If you're adding a new class of macromolecules, it quite likely you
want to add a new :class:`.Topology` class. See the
:mod:`.topologies.base` for guidance on adding these. The topology
class does the assembly of the macromolecule from the building blocks.

"""

import logging
import rdkit.Chem.AllChem as rdkit
from collections import Counter, defaultdict
from inspect import signature

from .molecule import Molecule
from .. import topologies
from ..functional_groups import FunctionalGroup
from ...utilities import OPTIONS

logger = logging.getLogger(__name__)


class Cached(type):
    """
    A metaclass for creating classes which create cached instances.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = dict()

    def __call__(self, *args, **kwargs):
        sig = signature(self.__init__)
        sig = sig.bind_partial(self, *args, **kwargs)
        sig.apply_defaults()
        sig = sig.arguments
        key = self.gen_key(sig['building_blocks'], sig['topology'])

        if key in self.cache and OPTIONS['cache']:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            if OPTIONS['cache']:
                self.cache[key] = obj
            return obj


class MacroMoleculeBuildError(Exception):
    ...


class MacroMolecule(Molecule, metaclass=Cached):
    """
    A representing assembled macromolecules.

    Because of the computational cost associated with macromolecule
    assembly, instances of this class are cached. This means that
    providing the same arguments to the initializer will not build a
    different instance with the same attribute values. It will yield
    the original instance, retrieved from memory.

    This class is not intended to be used directly but should be
    inherited by subclasses representing specific macromolecules. The
    :class:`Cage` and :class:`Polymer` classes are examples of this.
    Any information or methods that apply generally to all
    macromolecule should be defined within this class while specific
    non-general data should be included in the derived classes.

    Attributes
    ----------
    building_blocks : :class:`list` of :class:`.StructUnit`
        This attribute holds :class:`.StructUnit` instances which
        represent the monomers forming the macromolecule. Only one
        :class:`.StructUnit` instance is needed per building block,
        even if multiples of that molecule join up to form the
        macromolecule.

    bb_counter : :class:`collections.Counter`
        A counter keeping track of how much of each building block is
        used to form the macromolecule. Added by
        :func:`.Topology.build`.

    topology : :class:`.Topology`
        Defines the shape of macromolecule and assembles it.

    bonds_made : :class:`int`
        The number of bonds made during assembly. Added by
        :func:`.Topology.build`.

    func_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The remnants of building block functional groups present in the
        molecule. These functional groups track which atoms belonged to
        functional groups in the building block molecules. The id of
        each :class:`.FunctionalGroup` should match its index.

    """

    def __init__(self,
                 building_blocks,
                 topology,
                 name="",
                 note="",
                 bb_conformers=None):
        """
        Initialize a :class:`MacroMolecule` instance.

        Parameters
        ---------
        building_blocks : :class:`list` of :class:`.StructUnit`
            The :class:`.StructUnit` instances of building blocks
            forming the macromolecule.

        topology : :class:`.Topology`
            Defines the shape of macromolecule and assembles it.

        name : :class:`str`, optional
            A name which can be given to the molcule for easy
            identification.

        note : :class:`str`, optional
            A note or comment about the molecule.

        bb_conformers : :class:`list` of :class:`int`, optional
            The ids of the building block conformers to be used. Must
            be equal in length to `building_blocks` and orders must
            correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        """

        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(building_blocks))]

        self.building_blocks = building_blocks
        self.topology = topology

        try:
            # Ask the ``Topology`` instance to assemble/build the
            # macromolecule. This creates the `mol` and `func_groups`
            # attributes.
            topology.build(self, bb_conformers)

        except Exception as ex:
            errormsg = ('Build failure.\n'
                        '\n'
                        'topology\n'
                        '--------\n'
                        '{}\n'
                        '\n'
                        'building blocks\n'
                        '---------------\n').format(topology)

            bb_blocks = []
            for i, bb in enumerate(building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append('{} {}\n{}'.format(
                    bb.__class__.__name__,
                    [info.name for info in bb.func_group_infos],
                    bb.mdl_mol_block(bb_conf)))

            errormsg += '\n'.join(bb_blocks)

            raise MacroMoleculeBuildError(errormsg) from ex

        self.func_groups = tuple(self.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(self.func_groups):
            func_group.id = id_

        super().__init__(name, note)
        self.save_rdkit_atom_props({'mol_index', 'bb_index'})

    def add_conformer(self, bb_conformers):
        """
        Assembles a new conformer.

        Parameters
        ----------
        bb_conformers : :class:`list` of :class:`int`
            The ids of the building block conformers to be used. Must
            be equal in length to :attr:`building_blocks` and orders
            must correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        Returns
        -------
        :class:`int`
            The id of the new conformer.

        """

        # Save the original rdkit molecule.
        original_mol = self.mol
        # Build a new molecule.
        try:
            # Ask the ``Topology`` instance to assemble/build the
            # macromolecule. This creates the `mol` and `func_groups`
            # attributes.
            self.topology.build(self, bb_conformers)

        except Exception as ex:
            self.mol = rdkit.Mol()
            errormsg = ('Build failure.\n'
                        '\n'
                        'topology\n'
                        '--------\n'
                        '{}\n'
                        '\n'
                        'building blocks\n'
                        '---------------\n').format(self.topology)

            bb_blocks = []
            for i, bb in enumerate(self.building_blocks):
                bb_conf = bb_conformers[i]
                bb_blocks.append('{} {}\n{}'.format(
                    bb.__class__.__name__,
                    [info.name for info in bb.functional_group_infos],
                    bb.mdl_mol_block(bb_conf)))

            errormsg += '\n'.join(bb_blocks)

            raise MacroMoleculeBuildError(errormsg) from ex

        self.func_groups = tuple(self.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(self.func_groups):
            func_group.id = id_

        # Get the new conformer.
        new_conf = rdkit.Conformer(self.mol.GetConformer())
        # Add it to the original molecule.
        new_id = original_mol.AddConformer(new_conf, True)
        self.mol = original_mol
        return new_id

    def building_block_cores(self, bb):
        """
        Yields the "cores" of the building block molecules.

        The structure of the yielded cores has the geometry found in
        the macromolecule.

        Parameters
        ----------
        bb : :class:`int`
            The index of a building block molecule within
            :attr:`building_blocks`. The cores of this molecule are
            yielded.

        Yields
        ------
        :class:`rdkit.Mol`
            The core of a building block molecule, as found in the
            macromolecule.

        """

        mols = defaultdict(set)
        for atom_id, props in self.atom_props.items():
            correct_bb = props.get('bb_index', float('nan')) == bb
            if correct_bb and self.is_core_atom(atom_id):
                mols[props['mol_index']].add(atom_id)

        for mol in mols.values():
            core = rdkit.EditableMol(self.mol)
            for atom in reversed(range(self.mol.GetNumAtoms())):
                if atom not in mol:
                    core.RemoveAtom(atom)

            yield core.GetMol()

    def json(self, include_attrs=None):
        """
        Returns a JSON representation of the molecule.

        The representation has the form

        .. code-block:: python

            {
                'class' : 'Polymer',
                'mol_block' : '''A string holding the V3000 mol
                                 block of the molecule.''',
                'building_blocks' : {bb1.json(), bb2.json()},
                'topology' : 'Copolymer(repeating_unit="AB")',
                'note' : 'A nice molecule.',
                'name' : 'Poly-Benzene',
                'atom_props': {0: {'prop1': 1.0,
                                   'prop2': 'value1'}}
            }

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecule to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        conformers = [
            (conf.GetId(), self.mdl_mol_block(conformer=conf.GetId()))
            for conf in self.mol.GetConformers()
        ]

        json = {
            'bb_counter': [(key.json(), val) for key, val in
                           self.bb_counter.items()],
            'bonds_made': self.bonds_made,
            'class': self.__class__.__name__,
            'conformers': conformers,
            'building_blocks': [
                x.json() for x in self.building_blocks
            ],
            'topology': repr(self.topology),
            'note': self.note,
            'name': self.name,
            'atom_props': self.atom_props,
            'func_groups': repr(self.func_groups)

        }

        json.update(
            {attr: repr(getattr(self, attr)) for attr in include_attrs}
        )
        return json

    @classmethod
    def _json_init(cls, json_dict):
        """
        Completes a JSON initialization.

        This function is not to be used. Use :meth:`Molecule.load`
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : :class:`dict`
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        d = dict(json_dict)
        d.pop('building_blocks')
        d.pop('class')

        bb_counter = Counter({Molecule.from_dict(key): val for
                              key, val in d.pop('bb_counter')})
        bbs = list(bb_counter)
        topology = eval(d.pop('topology'),  topologies.__dict__)

        key = cls.gen_key(bbs, topology)
        if key in cls.cache and OPTIONS['cache']:
            return cls.cache[key]

        obj = cls.__new__(cls)

        (conf_id, mol_block), *confs = d.pop('conformers')
        obj.mol = rdkit.MolFromMolBlock(molBlock=mol_block,
                                        sanitize=False,
                                        removeHs=False)
        obj.mol.GetConformer().SetId(conf_id)

        for conf_id, mol_block in confs:
            conf_mol = rdkit.MolFromMolBlock(molBlock=mol_block,
                                             sanitize=False,
                                             removeHs=False)
            conf = conf_mol.GetConformer()
            conf.SetId(conf_id)
            obj.mol.AddConformer(conf)

        obj.topology = topology
        obj.bb_counter = bb_counter
        obj.bonds_made = d.pop('bonds_made')
        obj.note = d.pop('note')
        obj.name = d.pop('name') if d.pop('load_names') else ''
        obj.key = key
        obj.building_blocks = bbs
        obj.atom_props = {int(key): value for key, value in
                          d.pop('atom_props').items()}
        # Globals for eval.
        g = {'FunctionalGroup': FunctionalGroup}
        obj.func_groups = tuple(eval(d.pop('func_groups'), g))
        if OPTIONS['cache']:
            cls.cache[key] = obj

        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        return obj

    @staticmethod
    def gen_key(building_blocks, topology):
        """
        Generates the key used for caching the molecule.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.StructUnit`
            The building blocks used to make the macromolecule.

        topology : :class:`.Topology`
            The topology used to make the macromolecule.

        Returns
        -------
        :class:`tuple`
            The key used for caching the macromolecule.

        """

        bb_keys = frozenset(x.key for x in building_blocks)
        return bb_keys, repr(topology)

    def bb_distortion(self, bb_conformers=None, conformer=-1):
        """
        Rmsd difference of building blocks before and after assembly.

        The function looks at each building block in the macromolecule
        and calculates the rmsd between the "free" version and the one
        present in the macromolecule. The mean of these rmsds is
        returned.

        Atoms which form the functional group of the building blocks
        and hydrogens are excluded from the calculation.

        Parameters
        ----------
        bb_conformers : :class:`list` of :class:`int`
            The ids of building block conformers to use. 1 id for each
            building block, in an order corresponding to
            :attr:`building_blocks`. If ``None``, all conformer ids
            default to ``-1``.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The mean rmsd of the macromole's building blocks to their
            "free" counterparts.

        """

        if bb_conformers is None:
            bb_conformers = [
                -1 for _ in range(len(self.building_blocks))
            ]

        # Go through each of the building blocks. For each building
        # block get the core. Get the corrospending cores in the
        # macromolecules and add the rmsd to the sum. Increment the
        # count to calculate the mean later.
        rmsd = 0
        n = 0
        for i, bb in enumerate(self.building_blocks):
            free = bb.core()
            am = [(x, x) for x in range(free.GetNumAtoms())]
            for frag in self.building_block_cores(i):
                rmsd += rdkit.AlignMol(free,
                                       frag,
                                       bb_conformers[i],
                                       conformer,
                                       atomMap=am)
                n += 1
        return rmsd / n

    def __str__(self):
        return "{}(building_blocks={}, topology={!r})".format(
                        self.__class__.__name__,
                        [str(x) for x in self.building_blocks],
                        self.topology)

    def __repr__(self):
        return str(self)
