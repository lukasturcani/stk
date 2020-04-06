def get_rank_fitness(population):
    population = sorted(
        population,
        key=lambda record: record.get_fitness_value(),
    )
    return {
        record: 1/rank
        for rank, record in enumerate(population, 1)
    }
