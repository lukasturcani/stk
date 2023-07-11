from dataclasses import dataclass
from typing import Any

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    fitness_normalizer: stk.FitnessNormalizer
    fitness_values: dict[stk.MoleculeRecord[Any], Any]
    normalized: dict[stk.MoleculeRecord[Any], Any]
