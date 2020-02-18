import attr


@attr.s(kw_only=True)
class MutationType(object):
    dominance_coeff = attr.ib(default=0.5, type=float)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(factory=lambda: [0, ])
    convert_to_substitution = attr.ib(default=True, type=bool)
    weight = attr.ib(default=0, type=float)


def neutral_mutation_frac(mutation_types):
    """
    Return the fraction of mutations in mutation_types that are neutral.
    """
    neutral_weight = 0
    total_weight = 0
    for mut_type in mutation_types:
        assert isinstance(mut_type, MutationType)
        if mut_type.distribution_type == "f" and mut_type.distribution_args[0] == 0:
            # we define neutral mutations as "f" mutations with s=0
            neutral_weight += mut_type.weight
        total_weight += mut_type.weight
    if total_weight == 0:
        neutral_frac = 1
    else:
        neutral_frac = neutral_weight / total_weight
    return neutral_frac


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

    @op.validator
    def _op_validator(self, _attr, op):
        if op not in self.op_types:
            raise ValueError(f"Invalid conditioning op `{op}`")

    @classmethod
    def op_id(cls, op):
        for i, _op in enumerate(cls.op_types):
            if op == _op:
                return i
        return -1
