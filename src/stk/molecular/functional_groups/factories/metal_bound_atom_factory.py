"""
Metal Bound Atom Factory
========================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import MetalBoundAtom


class MetalBoundAtomFactory(FunctionalGroupFactory):
    """
    Creates :class:`.MetalBoundAtom` instances.

    Creates functional groups from substructures, which match the
    ``[metal]~[atom]`` functional group string.

    Examples
    --------
    * Creating MetalBoundAtom Functional Groups *

    You want to create a building block which has
    :class:`.MetalBoundAtom` functional groups for each :attr:`smiles`
    in the molecule.

    .. code-block:: python

        import stk

        # Build single atoms to place on metal centre topology.
        # Metal atom with 6 functional groups.
        atom = rdkit.MolFromSmiles('[Fe+2]')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        metal_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
        atom_0 = list(metal_atom.get_atoms())[0]
        metal_atom = metal_atom.with_functional_groups(
            (stk.SingleAtom(atom_0) for i in range(6))
        )


        # Nitrogen atom.
        atom = rdkit.MolFromSmiles('N')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        binding_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
        atom_0 = list(binding_atom.get_atoms())[0]
        binding_atom = binding_atom.with_functional_groups(
            (stk.SingleAtom(atom_0) for i in range(1))
        )


        # Build an Fe atom with octahedrally coordinated N atoms.
        metal_centre = stk.ConstructedMolecule(
            topology_graph=stk.metal_centre.Octahedral(
                building_blocks={
                    metal_atom: 0,
                    binding_atom: range(1, 7)
                }
            )
        )

        # Define a building block with 6 new MetalBoundAtom functional
        # groups.
        metal_centre = stk.BuildingBlock.init_from_molecule(
            metal_centre,
            functional_groups=[
                stk.MetalBoundAtomFactory(
                    atom_smarts='[N]'
                    metal_smarts='[Fe]',
                )
            ]
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, atom_smarts, metal_smarts):
        """
        Initialize a :class:`.MetalBoundAtomFactory` instance.

        Parameters
        ----------
        atom_smarts : :class:`str`
            SMARTS defining the *bonder* atom.

        metal_smarts : :class:`str`
            SMARTS defining the *metal* atom.

        """

        self._smarts = f'[{metal_smarts}]~[{atom_smarts}]'

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids(self._smarts, molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield MetalBoundAtom(
                atom=atoms[1],
                metal=atoms[0],
            )
