"""
Defines classes that describe rotaxanes.

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

class Rotaxane(MacroMolecule):
    """
    Used to represent rotaxanes.

    """
    pass
