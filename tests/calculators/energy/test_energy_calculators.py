import stk


def test_mmff(amine2, amine2_conf1):
    mmff = stk.MMFFEnergy()
    assert mmff.get_energy(amine2) < mmff.get_energy(amine2_conf1)


def test_uff(amine2, amine2_conf1):
    uff = stk.UFFEnergy()
    assert uff.get_energy(amine2) < uff.get_energy(amine2_conf1)


def test_cache_use(amine2):
    mmff = stk.MMFFEnergy()
    mmff.get_energy(amine2)
    # Since use_cache is False the cache should be empty.
    assert not mmff._cache

    # To test that the cache is not being used, put a random object
    # into it, and test that it was not returned.
    obj = object()
    mmff._cache[amine2] = obj
    assert mmff.get_energy(amine2) is not obj

    # Test that the cache is being filled when use_cache is True.
    cached_mmff = stk.MMFFEnergy(use_cache=True)
    assert not cached_mmff._cache
    cached_mmff.get_energy(amine2)
    assert cached_mmff._cache

    # Test that the cache is being used by putting a random object into
    # it and making sure it gets returned.
    cached_mmff._cache[amine2] = obj
    assert cached_mmff.get_energy(amine2) is obj


def test_formation(polymer, amine2, water):
    mmff = stk.MMFFEnergy(use_cache=True)
    building_blocks = list(polymer.building_block_vertices.keys())

    formation = stk.FormationEnergy(
        energy_calculator=mmff,
        reactants=building_blocks,
        products=[water]*3
    )
    reactant_energy = sum(
        mmff.get_energy(bb) for bb in building_blocks
    )
    product_energy = (
        mmff.get_energy(water)*3 + mmff.get_energy(polymer)
    )
    formation_energy = product_energy - reactant_energy
    assert formation.get_energy(polymer) - formation_energy < 1e-4
