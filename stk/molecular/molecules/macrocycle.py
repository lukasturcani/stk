"""
Defines classes that describe macrocycles.

There is a family of classes dealing with macrocycles, a cyclic
oligomer topology used to construct macrocycles.

"""

from .macro_molecule import MacroMolecule
from .struct_unit import StructUnit
import rdkit.Chem.AllChem as rdkit
import os


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

        Returns
        -------
        :class:`.list` of :class:`.int`
            Atom ids of the atoms comprising the largest ring.

        """

        ssr = rdkit.GetSymmSSSR(self.mol)
        ring_atom_ids = list(max(ssr, key=len))

        return ring_atom_ids

    def cycle_coords(self, path=None, conformer=-1):
        """
        Find the coordinates of the macrocyclic atoms in the molecule.

        Coordinates of the atoms comprising the largest ring in the
        macrocycle are found and the xyz coordinates file containing
        only those atoms can be saved.

        Notes
        -----
        The approach identifies the Smallest Set of Symmetric Rings
        and as a result one of multiple rings of the same size can be
        chose arbitrarily, making the results not unique. This should
        not be a problem in most applications.

        Parameters
        ----------
        path : :class:`.str`
            A path where the xyz file should be saved. If ``None`` then
            no file is produced. InChKey is used as the filename.

        Returns
        -------
        :class:`list` of :class:`list` of :class:`float`
            Coordinates of the atoms in the largest ring in the format
            ``[atom_index, x, y, z]``.

        """

        ssr = rdkit.GetSymmSSSR(self.mol)
        conf = self.mol.GetConformer(conformer)
        macrocycle = (self.mol.GetAtomWithIdx(i)
                      for i in max(ssr, key=len))
        macro_coords = [[atom.GetIdx(),
                         *conf.GetAtomPosition(atom.GetIdx())]
                        for atom in macrocycle]

        if path is not None:
            name = rdkit.MolToInchiKey(self.mol)
            xyz_file = f'{len(macro_coords)}\n\n'

            for anum, *coords in macro_coords:
                xyz_file += f'{anum} {coords[0]} {coords[1]} '
                xyz_file += f'{coords[2]}\n'

            if not os.path.exists(path):
                os.makedirs(path)

            with open(f'{path}/{name}.xyz', 'w') as f:
                f.write(xyz_file)

        return macro_coords


class MacrocycleStructUnit(MacrocycleBase, StructUnit):
    """
    Used to represent macrocyles loaded as :class:`.StructUnit`s.

    """
    pass


class Macrocycle(MacrocycleBase, MacroMolecule):
    """
    Used to represent macrocycles constructed by ``stk``.

    """
    pass
