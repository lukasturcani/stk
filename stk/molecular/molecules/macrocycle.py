"""
Defines classes that describe macrocycles.

There is a family of classes dealing with macrocycles, a cyclic
oligomer topology used to construct macrocycles.

"""

from .macro_molecule import MacroMolecule
from .struct_unit import StructUnit
import rdkit.Chem.AllChem as rdkit


class MacrocycleBase:
    """
    Used to represent macrocycles.

    Macrocycles are molecules that contain a large cycle. A simple
    example is a polymer with two ends connected together. This base
    class allows for the macrocyles to be initialised as either
    :class:`.StructUnit` or :class:`.MacroMolecule` but to be equally
    identified as macrocycles.

    """

    def cycle_atoms(self, conformer=-1):
        """
        Find the macrocyclic atoms in the molecule.

        Notes
        -----
        The approach identifies the Smallest Set of Symmetric Rings
        and as a result one of multiple rings of the same size can be
        chose arbitrarily, making the results not unique. This should
        not be a problem in most applications.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`list` of :class:`int`
            Atom ids of the atoms comprising the largest ring.

        """

        ssr = rdkit.GetSymmSSSR(self.mol)
        return list(max(ssr, key=len))


class MacrocycleStructUnit(MacrocycleBase, StructUnit):
    """
    Used to represent macrocyles loaded as :class:`.StructUnit`.

    """
    pass


class Macrocycle(MacrocycleBase, MacroMolecule):
    """
    Used to represent macrocycles constructed by ``stk``.

    """
    pass
