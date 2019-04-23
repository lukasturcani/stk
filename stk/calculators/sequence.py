from .optimizers import Optimizer


class Sequence(Optimizer):
    """


    """

    def __init__(self, *calculators):
        """
        Initializes a :class:`Sequence` instance.

        Parameters
        ----------
        *calculators : :class:`Calculator`
            The calculators which get used in sequence.

        """

        self.calculators = calculators
        # skip_optimized is False because it is the job of the
        # optimizers in calculators to toggle skipping for themselves.
        super().__init__(skip_optimized=False)

    def optimize(self, mol, conformer=-1):
        """
        Chains multiple :class:`Optimizer` instances together.

        Attributes
        ----------
        optimizers : :class:`tuple` of :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order in this :class:`tuple`.

        Examples
        --------
        Let's say we want to embed a molecule with ETKDG first and then
        minimize it with the MMFF force field.

        .. code-block:: python

            import rdkit.Chem.AllChem as rdkit
            mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
            etkdg = RDKitEmbedder(rdkit.ETKDG())
            mmff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
            optimizer = OptimizerPipeline(etkdg, mmff)
            optimizer.optimize(mol)



        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self.optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol.name}".')
            optimizer.optimize(mol, conformer)
