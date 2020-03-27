
class MoleculePresent(Terminator):
    """
    Stops the EA if a specific molecule has been found.

    """

    def __init__(self, mol):
        """
        Initialize a :class:`MoleculePresent` instance.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule which if present in any of the EA's generations
            causes it to stop running.

        """

        self._mol = mol
        self._is_same_molecule = partial(
            self._is_same_molecule,
            mol.to_rdkit_mol()
        )

    def terminate(self, progress):
        """
        Return ``True`` if :attr:`mol` is in `progress`.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if :attr:`mol` in `progress`, ``False`` otherwise.

        """

        # Check for the presence of the molecule, starting with the
        # newest generation first.
        for pop in progress.subpopulations:
            if any(
                self._is_same_molecule(mol.to_rdkit_mol())
                for mol in pop
            ):
                return True
        return False

    @staticmethod
    def _is_same_molecule(mol1, mol2):
        return rdkit.MolToInchi(mol1) == rdkit.MolToInchi(mol2)


