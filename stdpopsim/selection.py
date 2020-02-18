import attr
import msprime


@attr.s(frozen=True, kw_only=True)
class MutationType(object):
    dominance_coeff = attr.ib(default=0.5, type=float)
    distribution_type = attr.ib(default="f", type=str)
    distribution_args = attr.ib(factory=lambda: [0, ])
    convert_to_substitution = attr.ib(default=True, type=bool)
    weight = attr.ib(default=0, type=float)


def KimDFE():
    """
    Kim et al. (2018), p.23, http://doi.org/10.1371/journal.pgen.1007741
    """
    neutral = MutationType(
            weight=1.0, dominance_coeff=0.5, distribution_type="f",
            distribution_args=[0, ])
    gamma_shape = 0.186  # shape, alpha
    gamma_mean = -0.01314833  # expected value, alpha*beta
    h = 0.5/(1-7071.07*gamma_mean)
    negative = MutationType(
            weight=2.31, dominance_coeff=h, distribution_type="g",
            distribution_args=[gamma_mean, gamma_shape])
    return [neutral, negative]


class DrawMutation(msprime.DemographicEvent):
    def __init__(self, time, mutation_type_id, population_id, coordinate):
        super().__init__("draw_mutation", time)
        self.mutation_type_id = mutation_type_id
        self.population_id = population_id
        self.coordinate = round(coordinate)


class ChangeMutationFitness(msprime.DemographicEvent):
    def __init__(
            self, time, mutation_type_id, population_id, selection_coeff,
            dominance_coeff):
        super().__init__("change_mutation_fitness", time)
        self.mutation_type_id = mutation_type_id
        self.population_id = population_id
        self.selection_coeff = selection_coeff
        self.dominance_coeff = dominance_coeff
