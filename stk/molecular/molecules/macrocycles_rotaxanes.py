"""
Defines classes that describe macrocycles and rotaxanes.

There is a family of classes dealing with macrocycles, a cyclic
oligomer topology used to construct macrocycles, and classes related
to rotaxanes.

A base class :class:`MacrocycleBase` contains the methods used by all
macrocycles, no matter whether loaded or constructed using ``stk``.
A daugter :class:`MacrocycleStructUnit` also inherits
:class:`StructUnit` and is used to load macrocycles that were not
constructed in ``stk``, while:class:`Macrocycle` inherits
:class:`MacroMolecule` and is a result of construction within the
``stk``. Either class can construct rotaxanes. The class
:class:`MacrocycleBase` defines a couple of methods that are useful
to generate properties of macrocycle.
:meth:`MacrocycleBase.macro_atoms` returns the coordinates and indices
of the atoms forming the largest ring in the macrocycle (used for
threading). The method relies on the Smallest Set of Symmetric Rings
and hence its result is not unique, but that should not cause any
problem for most applications. A vector normal to the plane of the
macrocycle is returned by :meth:`MacrocycleBase.macrocycle_plane'.

One way to construct macrocycles is to use :meth:`Cyclic`, which
behaves analogously to construction of linear polymers but the monomers
are placed on a circumference of a large circle and the two terminal
monomers are joined together to close the macrocycle. Otherwise it has
attributes akin to :class:`Linear`, i.e. :attr:`Cyclic.repeating_unit`,
:attr:`Cyclic.orientation`, and :attr:`Cyclic.n`.

A new :class:`MacroMolecule` called :class:`Rotaxane` is defined. This
class stores constructed rotaxanes, normally originating from
:class:`NRotaxane` topology. The :class:`NRotaxane` takes the axle and
a :class:`list` of :class:`MacrocycleBase` objects to be threaded.
The macrocycles are spaced evenly along the axle in a direction
specified analogously to polymers in :attr:`NRotaxane.orientation`.
This allows to construct mechanical isomers, with opposite orientations
of the macrocycle relative to the axle caps.

"""

from .macro_molecule import MacroMolecule
from .struct_unit import StructUnit


class MacrocycleBase:
    """
    Used to represent macrocycles.

    Macrocycles are molecules that contain a large cycle. A simple
    example is a polymer with two ends connected together. This base
    class allows for the macrocyles to be initialised as either
    :class:`StructUnit` or :class:`MacroMolecule` but to be equally
    identified as macrocycles.

    """

    def macro_atoms(self, xyz=None, conformer=-1):
        """
        Find the macrocyclic atoms in the molecule.

        Coordinates of the atoms comprising the largest ring in the
        macrocycle are returned and the xyz coordinates file containing
        only those atoms can be saved. The method uses
        :meth:`rdkit.GetSymSSSR()` to identify the Smallest Set of
        Symmetric Rings, so the results are not unique. This should not
        be a problem in most applications.

        Parameters
        ----------
        xyz : :class:`str`
            A path where the xyz file should be saved. If None then
            no file is produced. InChKey is used as the filename.

        Returns
        -------
        :class:`list` of :class:`list` of :class:`float`
            Coordinates of the atoms in the largest ring in the format
            [atomic_number, *xyz].

        :class:`list` of :class:`int`
            Atom ids of the atoms comprising the largest ring.

        """
        ssr = rdkit.GetSymmSSSR(self.mol)
        conf = self.mol.GetConformer(conformer)
        ring_atom_ids = list(max(ssr, key=len))
        macrocycle = (self.mol.GetAtomWithIdx(i)
                      for i in ring_atom_ids)
        macro_coords = [[atom.GetAtomicNum(),
                         *conf.GetAtomPosition(atom.GetIdx())]
                        for atom in macrocycle]

        if xyz is not None:
            name = rdkit.MolToInchiKey(self.mol)
            xyz_file = f'{len(macro_coords)}\n\n'

            for anum, *coords in macro_coords:
                xyz_file += f'{anum} {coords[0]} {coords[1]} '
                xyz_file += f'{coords[2]}\n'

            if not os.path.exists(xyz):
                os.makedirs(xyz)

            with open(f'{xyz}/{name}.xyz', 'w') as f:
                f.write(xyz_file)

        return macro_coords, ring_atom_ids


class MacrocycleStructUnit(MacrocycleBase, StructUnit):
    """
    Used to represent macrocyles loaded as :class:`StructUnit`s.

    """
    pass


class Macrocycle(MacrocycleBase, MacroMolecule):
    """
    Used to represent macrocycles constructed by ``stk``.

    """
    pass
