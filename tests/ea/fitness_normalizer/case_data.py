from dataclasses import dataclass
from typing import Any

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    fitness_normalizer: stk.FitnessNormalizer[stk.MoleculeRecord[Any]]
    fitness_values: dict[stk.MoleculeRecord[Any], Any]
    normalized: dict[stk.MoleculeRecord[Any], Any]

    @staticmethod
    def new(
        fitness_normalizer: stk.FitnessNormalizer[stk.MoleculeRecord[Any]],
        fitness_values: dict[stk.MoleculeRecord[Any], tuple[Any, Any]],
    ) -> "CaseData":
        return CaseData(
            fitness_normalizer=fitness_normalizer,
            fitness_values={
                record: fitness_value[0]
                for record, fitness_value in fitness_values.items()
            },
            normalized={
                record: fitness_value[1]
                for record, fitness_value in fitness_values.items()
            },
        )
