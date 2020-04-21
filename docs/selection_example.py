import random
import logging

import stdpopsim

logger = logging.getLogger(__name__)


def adaptive_introgression(seed):
    """
    Adaptive introgression in HomSap/PapuansOutOfAfrica_10J19.

    A neutral mutation is drawn in Denisovans, transmitted to Papuans via a
    migration pulse, and is then positively selected in the Papuan population.
    The time of mutation introduction, the time of selection onset, and the
    selection coefficient, are each random variables.
    """
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("PapuansOutOfAfrica_10J19")
    contig = species.get_contig("chr1", length_multiplier=0.001)
    samples = model.get_samples(
            100,  # YRI
            0,    # CEU
            0,    # CHB
            100,  # Papuan
            2,    # DenA
            2)    # NeaA

    # One mutation type, which we'll use for the positively selected mutation.
    # Neutral mutations will be added by the SLiM engine as usual, after the
    # SLiM phase of the simulation has completed.
    positive = stdpopsim.ext.MutationType(convert_to_substitution=False)
    mutation_types = [positive]
    mut_id = len(mutation_types)

    # We need some demographic model parameters to set bounds on the timing
    # of random variables and extended_events (below).
    # These values were copied from the PapuansOutOfAfrica_10J19 model
    # implementation, but can also be found in the catalog documentation.
    T_Den_Nea_split = 15090
    T_DenA_Den1_split = 9750
    T_Den1_Papuan_mig = 29.8e3 / model.generation_time
    # The drawn mutation is transmitted via Den1.
    T_Den_split = T_DenA_Den1_split
    T_mig = T_Den1_Papuan_mig

    # Draw random variables.
    rng = random.Random(seed)
    # Time of mutation introduction.
    T_mut = rng.uniform(T_Den_split, T_Den_Nea_split)
    # Time of selection onset. We use a lower bound of 1000 years ago, so that
    # weaker selection has time to act.
    T_sel = rng.uniform(1000 / model.generation_time, T_mig)
    # Selection coefficient.
    s = rng.uniform(0.001, 0.1)
    logger.info(f"Parameters: T_mut={T_mut:.3f}, T_sel={T_sel:.3f}, s={s:.3g}")

    # Place the drawn mutation in the middle of the contig.
    coordinate = round(contig.recombination_map.get_length() / 2)

    pop = {p.id: i for i, p in enumerate(model.populations)}

    # Thinking forwards in time, we define a number of extended events that
    # correspond to drawing the mutation, conditioning on the new allele not
    # being lost, and it's selection. The extended events will be time-sorted
    # by the SLiM engine, so the ordering here is only to aid clarity.
    # Like msprime.DemographicEvents, all the times used here have units of
    # generations before present.
    extended_events = [
        # Draw mutation in DenA.
        stdpopsim.ext.DrawMutation(
                time=T_mut, mutation_type_id=mut_id, population_id=pop["DenA"],
                coordinate=coordinate,
                # Save state before the mutation is introduced.
                save=True),
        # Because the drawn mutation is neutral at the time of introduction,
        # it's likely to be lost due to drift. To avoid this, we condition on
        # the mutation having AF > 0 in DenA. If this condition is false at any
        # time between start_time and end_time, the simulation will be
        # restored to the most recent save point.
        # Conditioning should start one generation after T_mut (not at T_mut!),
        # to avoid checking for the mutation before SLiM can introduced it.
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_mut),
                end_time=T_Den_split,
                mutation_type_id=mut_id, population_id=pop["DenA"],
                op=">", allele_frequency=0),
        # Denisovans split into DenA and Den1 at time T_Den_split,
        # so now we condition on having AF > 0 in Den1.
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_Den_split),
                end_time=T_mig,
                mutation_type_id=mut_id, population_id=pop["Den1"],
                op=">", allele_frequency=0,
                # Update save point at start_time (if the condition is met).
                save=True),
        # The Den1 lineage has migrants entering the Papaun lineage at T_mig,
        # so condition on AF > 0 in Papuans.
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_mig), end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                op=">", allele_frequency=0,
                # Update save point at start_time (if the condition is met).
                # If the Den1 migrants didn't carry the mutation, we will
                # instead restore to the previous save point.
                save=True),
        # The mutation is positively selected in Papuans at T_sel.
        # Note that this will have no effect, unless/until a mutation with the
        # specified mutation_type_id is found in the population.
        stdpopsim.ext.ChangeMutationFitness(
                start_time=T_sel, end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                selection_coeff=s, dominance_coeff=0.5),
        # Condition on AF > 0.05 in Papuans at the end of the simulation.
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0, end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                op=">", allele_frequency=0.05),
        ]

    # Simulate.
    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
            model, contig, samples,
            seed=rng.randrange(1, 2**32),
            mutation_types=mutation_types,
            extended_events=extended_events,
            slim_scaling_factor=10,
            slim_burn_in=0.1,
            # Set slim_script=True to print the script instead of running it.
            # slim_script=True,
            )

    return ts, T_mut, T_sel, s


if __name__ == "__main__":
    import sys
    import stdpopsim.cli
    import collections

    if len(sys.argv) == 2:
        seed = int(sys.argv[1])
    else:
        seed = 1234

    # Setup logging at verbosity level 2.
    args = collections.namedtuple("_", ["quiet", "verbosity"])(False, 2)
    stdpopsim.cli.setup_logging(args)

    ts, T_mut, T_sel, s = adaptive_introgression(seed)
