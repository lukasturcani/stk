import stk


def test_mmff(amine2):
    mmff = stk.MMFFEnergy()
    assert mmff.energy(amine2) < mmff.energy(amine2, 1)

    # To test that the cache is not being used, put a random object
    # into it, and test that it was not returned.
    obj = object()
    mmff.cache[(amine2, 1)] = obj
    assert mmff.energy(amine2, 1) is not obj

    # Test that the cache is being used by putting a random object into
    # it and making sure it gets returned.
    cached_mmff = stk.MMFFEnergy(use_cache=True)
    cached_mmff.cache[(amine2, 1)] = obj
    assert mmff.energy(amine2, 1) is obj


def test_uff(amine2):
    uff = stk.UFFEnergy()
    assert uff.energy(amine2) < uff.energy(amine2, 1)

    # To test that the cache is not being used, put a random object
    # into it, and test that it was not returned.
    obj = object()
    uff.cache[(amine2, 1)] = obj
    assert uff.energy(amine2, 1) is not obj

    # Test that the cache is being used by putting a random object into
    # it and making sure it gets returned.
    cached_uff = stk.UFFEnergy(use_cache=True)
    cached_uff.cache[(amine2, 1)] = obj
    assert uff.energy(amine2, 1) is obj


def test_formation(polymer, amine2):
    mmff = stk.MMFFEnergy(use_cache=True)

    water = stk.StructUnit.smiles_init('[H]O[H]')
    products = [water]*3
    formation = stk.FormationEnergy(
                        energy_calculator=mmff,
                        reactants=polymer.building_blocks,
                        products=products)

    reactant_energy = sum(
        mmff.energy(bb) for bb in polymer.building_blocks
    )
    product_energy = mmff.energy(water)*3 + mmff.energy(polymer)
    formation_energy = product_energy - reactant_energy
    assert formation.energy(polymer) - formation_energy < 1e-4
