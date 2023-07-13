import typing

T = typing.TypeVar("T")


def get_rank_fitness(population: dict[T, float]) -> dict[T, float]:
    sorted_population = sorted(
        zip(population.values(), population.keys()),
        reverse=True,
    )
    return {
        record: 1 / rank
        for rank, (_, record) in enumerate(sorted_population, 1)
    }
