Functional Groups
=================

.. toctree::
  :maxdepth: 1

  Overview and Examples <_autosummary/stk.FunctionalGroup>
  Alcohol <_autosummary/stk.Alcohol>
  Aldehyde <_autosummary/stk.Aldehyde>
  Alkene <_autosummary/stk.Alkene>
  Alkyne <_autosummary/stk.Alkyne>
  Amide <_autosummary/stk.Amide>
  Boronic Acid <_autosummary/stk.BoronicAcid>
  Bromo <_autosummary/stk.Bromo>
  Carboxylic Acid <_autosummary/stk.CarboxylicAcid>
  Dibromo <_autosummary/stk.Dibromo>
  Difluoro <_autosummary/stk.Difluoro>
  Diol <_autosummary/stk.Diol>
  Fluoro <_autosummary/stk.Fluoro>
  Generic Functional Group <_autosummary/stk.GenericFunctionalGroup>
  Iodo <_autosummary/stk.Iodo>
  Primary Amino <_autosummary/stk.PrimaryAmino>
  Ring Amine <_autosummary/stk.RingAmine>
  Secondary Amino <_autosummary/stk.SecondaryAmino>
  Single Atom <_autosummary/stk.SingleAtom>
  Thioacid <_autosummary/stk.Thioacid>
  Thiol <_autosummary/stk.Thiol>

Functional groups define which atoms of a :class:`.BuildingBlock` are
modified during :class:`.ConstructedMolecule` construction, and which
are used to position it.
The class of a :class:`.FunctionalGroup`
affects which :class:`.Reaction` can be used with it.
See the abstract base class :class:`.FunctionalGroup` for more
information.
