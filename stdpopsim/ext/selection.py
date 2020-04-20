import attr


@attr.s(kw_only=True)
class MutationType(object):
    dominance_coeff = attr.ib(default=0.5, type=float)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(factory=lambda: [0])
    convert_to_substitution = attr.ib(default=True, type=bool)
    # MutationTypes with non-zero weight will be simulated by SLiM,
    # using rates obtained from the relative weights of the types.
    # I.e. the weights will be used in the ``proportion`` parameter
    # to SLiM's :func:`initializeGenomicElementType()`.
    weight = attr.ib(default=0, type=float)


def slim_mutation_frac(mutation_types):
    """
    The fraction of mutations that should be added by SLiM.
    The remainder are added by msprime.mutate() once the SLiM part
    of the simulation is complete.
    """
    weighted = any(mut_type.weight > 0 for mut_type in mutation_types)
    return 1 if weighted else 0


def KimDFE():
    """
    Return neutral and negative MutationType()s representing a human DFE.
    Kim et al. (2018), p.23, http://doi.org/10.1371/journal.pgen.1007741
    """
    neutral = MutationType(weight=1.0)
    gamma_shape = 0.186  # shape, alpha
    gamma_mean = -0.01314833  # expected value, alpha*beta
    h = 0.5/(1-7071.07*gamma_mean)
    negative = MutationType(
            weight=2.31, dominance_coeff=h, distribution_type="g",
            distribution_args=[gamma_mean, gamma_shape])
    return [neutral, negative]


@attr.s
class GenerationAfter(object):
    time = attr.ib(type=float)


class ExtendedEvent(object):
    pass


@attr.s(kw_only=True)
class DrawMutation(ExtendedEvent):
    time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    coordinate = attr.ib(type=int)
    save = attr.ib(type=bool, default=False)


@attr.s(kw_only=True)
class ChangeMutationFitness(ExtendedEvent):
    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    selection_coeff = attr.ib(type=float)
    dominance_coeff = attr.ib(type=float)


@attr.s(kw_only=True)
class ConditionOnAlleleFrequency(ExtendedEvent):
    start_time = attr.ib(type=float)
    end_time = attr.ib(type=float)
    mutation_type_id = attr.ib(type=int)
    population_id = attr.ib(type=int)
    op = attr.ib(type=str, default=None)
    allele_frequency = attr.ib(type=float, default=None)
    save = attr.ib(type=bool, default=False)

    op_types = ("<", "<=", ">", ">=")

    def __attrs_post_init__(self):
        if self.op not in self.op_types:
            raise ValueError(f"Invalid conditioning op `{self.op}`")
        if not (0 <= self.allele_frequency <= 1):
            raise ValueError("Must have 0 <= allele_frequency <= 1.")
        if ((self.op == "<" and self.allele_frequency == 0) or
           (self.op == ">" and self.allele_frequency == 1)):
            raise ValueError(
                    f"allele_frequency {self.op} {self.allele_frequency} "
                    f"would require the allele frequency to be out of range.")

    @classmethod
    def op_id(cls, op):
        for i, _op in enumerate(cls.op_types):
            if op == _op:
                return i
        return -1
