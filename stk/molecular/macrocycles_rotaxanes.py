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

The module also introduces extensions to :class:`Molecule`, namely
:meth:`Molecule.direction` which calculates the best linear direction
of the molecule, and :meth:`Molecule.plane_normal` which returns a
vector normal to the plane of the molecule (or a subset of its atoms).

"""
