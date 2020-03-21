import rdkit.Chem.AllChem as rdkit


def get_inchi_key(molecule):
    """
    Get the InChIKey of `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose InChIKey is needed.

    Returns
    -------
    :class:`str`
        The InChIKey.

    """

    return rdkit.MolToInchiKey(molecule.to_rdkit_mol())


class GetBuildingBlockKey:
    """

    """

    def __init__(self, molecule_key, functional_group_key):
        """

        """

        self._molecule_key = molecule_key
        self._functional_group_key = functional_group_key

    def __call__(self, building_block):
        """

        """

        functional_group_keys = '-'.join(
            self._functional_group_key(functional_group)
            for functional_group
            in building_block.get_functional_groups()
        )
        placer_ids = ''.join(
            str(id_) for id_ in building_block.get_placer_ids()
        )
        return (
            f'{self._molecule_key(building_block)}-'
            f'{functional_group_keys}-'
            f'{placer_ids}'
        )
