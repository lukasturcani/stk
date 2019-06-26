"""
Defines electronic property calculators.

See :mod:`.electronic_property`.

"""


class ElectronicPropertyCalculator:
    def dipole_moment(self, mol, conformer=-1):
        raise NotImplementedError()

    def electron_affinity(self, mol, conformer=-1):
        raise NotImplementedError()

    def ionization_potential(self, mol, conformer=-1):
        raise NotImplementedError()
