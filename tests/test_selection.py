import numpy as np

import stdpopsim


def ai():
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("PapuansOutOfAfrica_10J19")
    contig = species.get_chunk(length=100*1000, genetic_map="HapMapII_GRCh37")
    samples = model.get_samples(200, 0,  0, 200, 2, 2)

    mutation_types = stdpopsim.KimDFE()
    positive = stdpopsim.MutationType(convert_to_substitution=False)
    mutation_types.append(positive)
    mut_id = len(mutation_types)

    # TODO: get these from the model itself
    T_Den_Nea_split = 15090
    T_Den1_Papuan_mig = 29.8e3 / model.generation_time
    # T_Den2_Papuan_mig = 45.7e3 / model.generation_time
    T_mig = T_Den1_Papuan_mig  # should flip a coin to choose Den1 or Den2.

    T_mut = np.random.uniform(T_mig, T_Den_Nea_split)
    T_sel = np.random.uniform(1e3 / model.generation_time, T_mig)
    s = np.random.uniform(0.001, 0.1)

    pop = {p.id: i for i, p in enumerate(model.populations)}
    coordinate = round(contig.recombination_map.get_length() / 2)
    draw_mut = stdpopsim.DrawMutation(T_mut, mut_id, pop["DenA"], coordinate)
    selection_onset = stdpopsim.ChangeMutationFitness(
            T_sel, mut_id, pop["Papuan"], s, 0.5)
    extended_events = [draw_mut, selection_onset]

    engine = stdpopsim.get_engine("slim")
    engine.simulate(
            model, contig, samples,
            verbosity=2,
            mutation_types=mutation_types,
            extended_events=extended_events,
            slim_script=True,
            # slim_no_recapitation=True,
            # slim_no_burnin=True,
            )


if __name__ == "__main__":
    ai()
