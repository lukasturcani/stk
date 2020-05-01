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
    ``[metal_smiles][atom_smiles]`` functional group string.

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
        metal_atom = stk.BuildingBlock.init_from_rdkit_mol(
            atom,
            functional_groups=[stk.SingleAtomFactory(
                num_functional_groups=6
            )]
        )

        # Nitrogen atom.
        atom = rdkit.MolFromSmiles('N')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        binding_atom = stk.BuildingBlock.init_from_rdkit_mol(
            atom,
            functional_groups=[stk.SingleAtomFactory(
                num_functional_groups=1
            )]
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
                    atom_smiles='[N]'
                    metal_smiles='[Fe]',
                )
            ]
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(self, atom_smiles, metal_smiles):
        """
        Initialize a :class:`.MetalBoundAtomFactory` instance.

        Parameters
        ----------
        atom_smiles : :class:`str`
            SMILES defining the *bonder* atom.

        metal_smiles : :class:`str`
            SMILES defining the *metal* atom.

        """

        self._smiles = f'[{metal_smiles}][{atom_smiles}]'

    def get_functional_groups(self, molecule):
        for atom_ids in _get_atom_ids(self._smiles, molecule):
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield MetalBoundAtom(
                atom=atoms[1],
                metal=atoms[0],
            )
