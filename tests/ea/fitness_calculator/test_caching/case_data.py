import typing
from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    fitness_calculator: stk.FitnessCalculator
    record: stk.MoleculeRecord
    fitness_value: typing.Any
