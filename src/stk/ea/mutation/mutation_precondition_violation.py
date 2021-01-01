class MutationPreconditionViolation(Exception):
    """
    Raised when a mutation fails due to a precondition violation.

    A precondition violation occurs when a molecule cannot undergo
    mutation, because it does not satisfy the necessary conditions
    for that mutation to take place. For example, if you try to
    mutate the hydrogen atoms of a molecule with no hydrogen atoms,
    that could be treated as a precondition violation. Note that
    in this example, it would also be reasonable to just do nothing,
    and return the unmutated molecule as its own mutant. It is up to
    the implementor of the mutation algorithm to decide what makes the
    most sense.

    Note that because an implementation of a
    :class:`.EvolutionaryAlgorithm` might not remove duplicates of a
    molecule in the population, returning a molecule as its own mutant
    may be undesirable, because it could lead to it being
    overrepresented in the population. Raising
    :class:`MutationPreconditionViolation` could be an effective
    way of avoiding such a pitfall.

    However, the default implementation of the
    :class:`.EvolutionaryAlgorithm` does actually remove duplicate
    molecules from the population, so you would avoid such an pitfall
    regardless of the approach taken.

    """

    pass
