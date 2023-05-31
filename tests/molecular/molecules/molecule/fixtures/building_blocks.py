import stk


def get_pd_atom() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=([0, 0, 0],),
    )


def get_fe_atom() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[Fe+2]",
        functional_groups=(
            stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(6)
        ),
        position_matrix=([0, 0, 0],),
    )


def get_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=C([H])C(C3=C([H])"
            "C([H])=NC([H])=C3[H])=C2[H])=C1[H]"
        ),
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


def get_other_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=NC([H])=C2[H])" "=C1[H]"
        ),
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


def get_palladium_bi_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[H]N([H])C([H])([H])C([H])([H])N([H])[H]",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#7]~[#6]",
                bonders=(0,),
                deleters=(),
            ),
        ],
    )


def get_palladium_cispbi_sqpl() -> stk.BuildingBlock:
    molecule = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
            metals={get_pd_atom(): 0},
            ligands={get_palladium_bi_1(): 0},
            reaction_factory=stk.DativeReactionFactory(
                reaction_factory=stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    },
                ),
            ),
        ),
    )

    return stk.BuildingBlock.init_from_molecule(
        molecule=molecule,
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[Pd]~[#7]",
                bonders=(0,),
                deleters=(),
                placers=(0, 1),
            ),
        ],
    )


def get_ditopic_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "[H]C1=C([H])C(C2=C([H])C([H])=C(Br)C([H])=C2[H])=C([H])"
            "C([H])=C1Br"
        ),
        functional_groups=[stk.BromoFactory()],
    )


def get_tritopic_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "[H]C1=C([H])C(N(C2=C([H])C([H])=C(Br)C([H])=C2[H])"
            "C2=C([H])C("
            "[H])=C(Br)C([H])=C2[H])=C([H])C([H])=C1Br"
        ),
        functional_groups=[stk.BromoFactory()],
    )


def get_tetratopic_linker() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles=(
            "[H]C1=C([H])C(C(C2=C([H])C([H])=C(Br)C([H])=C2[H])C"
            "(C2=C([H])C([H])=C(Br)C([H])=C2[H])C2=C([H])C([H])=C"
            "(Br)C([H])=C2[H])=C([H])C([H])=C1Br"
        ),
        functional_groups=[stk.BromoFactory()],
    )


def get_complex_ligand() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="[H]C1=NC(C([H])=NBr)=C([H])C([H])=C1[H]",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


def get_iron_complex() -> stk.BuildingBlock:
    iron_complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.OctahedralDelta(
            metals={get_fe_atom(): 0},
            ligands={get_complex_ligand(): (0, 1, 2)},
            reaction_factory=stk.DativeReactionFactory(
                reaction_factory=stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    }
                )
            ),
        )
    )
    return stk.BuildingBlock.init_from_molecule(
        molecule=iron_complex, functional_groups=[stk.BromoFactory()]
    )


def get_iron_bi_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="BrN=Cc1ccccn1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )


def get_mo_1() -> stk.BuildingBlock:
    return stk.BuildingBlock(
        smiles="c1cc2c(cn1)CCCCC2",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )
