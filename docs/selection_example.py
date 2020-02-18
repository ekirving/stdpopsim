import sys
import random

import stdpopsim


def adaptive_introgression(seed, verbosity=2):
    """
    Adaptive introgression in HomSap/PapuansOutOfAfrica_10J19.

    A neutral mutation is drawn in Denisovans, transmitted to Papuans via a
    migration pulse, and is then positively selected in the Papuan population.
    The time of mutation introduction, the time of selection onset, and the
    selection coefficient are each random variables.
    """
    rng = random.Random(seed)

    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("PapuansOutOfAfrica_10J19")
    contig = species.get_contig("chr1", length_multiplier=0.001)
    samples = model.get_samples(100, 0,  0, 100, 2, 2)

    mutation_types = stdpopsim.ext.KimDFE()
    positive = stdpopsim.ext.MutationType(convert_to_substitution=False)
    mutation_types.append(positive)
    mut_id = len(mutation_types)

    # TODO: get these from the model itself
    T_Den_Nea_split = 15090
    T_DenA_Den1_split = 9750
    # T_DenA_Den2_split = 12500
    T_Den1_Papuan_mig = 29.8e3 / model.generation_time
    # T_Den2_Papuan_mig = 45.7e3 / model.generation_time

    # The drawn mutation is transmitted via Den1.
    T_Den_split = T_DenA_Den1_split
    T_mig = T_Den1_Papuan_mig

    T_mut = rng.uniform(T_Den_split, T_Den_Nea_split)
    T_sel = rng.uniform(1e3 / model.generation_time, T_mig)
    s = rng.uniform(0.001, 0.1)

    pop = {p.id: i for i, p in enumerate(model.populations)}
    coordinate = round(contig.recombination_map.get_length() / 2)

    extended_events = [
        # Draw mutation in Denisovans.
        stdpopsim.ext.DrawMutation(
                time=T_mut, mutation_type_id=mut_id, population_id=pop["DenA"],
                coordinate=coordinate,
                # Save state before the mutation is introduced
                save=True),
        # Mutation is positively selected in Papuans
        stdpopsim.ext.ChangeMutationFitness(
                start_time=T_sel, end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                selection_coeff=s, dominance_coeff=0.5),
        # Allele frequency conditioning. If the condition is not met, we
        # restore to the most recent save point.
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_mut),
                end_time=T_Den_split,
                mutation_type_id=mut_id, population_id=pop["DenA"],
                op=">", allele_frequency=0),
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_Den_split),
                end_time=T_mig,
                mutation_type_id=mut_id, population_id=pop["Den1"],
                op=">", allele_frequency=0,
                # Update save point at start_time.
                save=True),
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=stdpopsim.ext.GenerationAfter(T_mig), end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                op=">", allele_frequency=0,
                # Update save point at start_time.
                save=True),
        stdpopsim.ext.ConditionOnAlleleFrequency(
                start_time=0, end_time=0,
                mutation_type_id=mut_id, population_id=pop["Papuan"],
                op=">", allele_frequency=0.1),
        ]

    engine = stdpopsim.get_engine("slim")
    ts = engine.simulate(
            model, contig, samples,
            seed=rng.randrange(1, 2**32),
            verbosity=verbosity,
            mutation_types=mutation_types,
            extended_events=extended_events,
            slim_scaling_factor=10,
            # slim_script=True,
            )

    return ts, T_mut, T_sel, s


if __name__ == "__main__":
    if len(sys.argv) == 2:
        seed = int(sys.argv[1])
    else:
        seed = 1234
    ts, T_mut, T_sel, s = adaptive_introgression(seed)
    print(f"Params: T_mut={T_mut}, T_sel={T_sel}, s={s}", file=sys.stderr)
