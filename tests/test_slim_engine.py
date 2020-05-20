"""
Tests for SLiM simulation engine.
"""
import os
import re
import io
import sys
import itertools
import unittest
import tempfile
import math
from unittest import mock

import tskit
import pyslim
import msprime

import stdpopsim
import stdpopsim.cli
from . test_cli import capture_output

IS_WINDOWS = sys.platform.startswith("win")


def slim_simulate_no_recap(seed=1234, **kwargs):
    """
    Return the tree sequence produced by SLiM, without recapitation, etc.
    """
    kwargs.update(slim_script=True)
    engine = stdpopsim.get_engine("slim")
    out, _ = capture_output(engine.simulate, **kwargs)

    # Find the name of the temporary trees_file in the script.
    match = re.search(r'"trees_file",\s*"([^"]*)"', out)
    assert match is not None
    tmp_trees_file = match.group(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        script_file = os.path.join(tmpdir, "script.slim")
        trees_file = os.path.join(tmpdir, "out.trees")
        # Write out the script with a new location for the trees_file.
        out = out.replace(tmp_trees_file, trees_file)
        with open(script_file, "w") as f:
            f.write(out)
        engine._run_slim(script_file, seed=seed)
        ts = pyslim.load(trees_file)
    return ts


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestAPI(unittest.TestCase):

    def test_bad_params(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)

        for scaling_factor in (0, -1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=scaling_factor,
                    dry_run=True)

        for burn_in in (-1, -1e-6):
            with self.assertRaises(ValueError):
                engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_burn_in=burn_in,
                    dry_run=True)

    def test_script_generation(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")

        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AncientEurasia_9K19")
        samples = model.get_samples(10, 20, 30, 40, 50, 60, 70)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

        model = species.get_demographic_model("AmericanAdmixture_4B11")
        samples = model.get_samples(10, 10, 10)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                slim_script=True)
        self.assertTrue("sim.registerLateEvent" in out)

    def test_recombination_map(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh37")
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=model, contig=contig, samples=samples,
                dry_run=True)

    def test_simulate(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("AraTha")
        contig = species.get_contig("chr5", length_multiplier=0.001)
        model = stdpopsim.PiecewiseConstantSize(species.population_size)
        samples = model.get_samples(10)
        ts = engine.simulate(
                demographic_model=model, contig=contig, samples=samples,
                slim_scaling_factor=10, slim_burn_in=0)
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

    def test_recap_and_rescale(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(10, 10, 10)

        for weight, seed in zip((0, 1), (1234, 2345)):
            if weight:
                mutation_types = [stdpopsim.ext.MutationType(weight=1)]
                extended_events = None
            else:
                mutation_types = None
                extended_events = []
            ts1 = engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    mutation_types=mutation_types, extended_events=extended_events,
                    slim_scaling_factor=10, slim_burn_in=0, seed=seed)
            ts2_headless = slim_simulate_no_recap(
                    demographic_model=model, contig=contig, samples=samples,
                    mutation_types=mutation_types, extended_events=extended_events,
                    slim_scaling_factor=10, slim_burn_in=0,
                    seed=seed)
            ts2 = engine.recap_and_rescale(
                    ts2_headless,
                    demographic_model=model, contig=contig, samples=samples,
                    mutation_types=mutation_types, extended_events=extended_events,
                    slim_scaling_factor=10, seed=seed)

            tables1 = ts1.dump_tables()
            tables2 = ts2.dump_tables()

            self.assertEqual(tables1.nodes, tables2.nodes)
            self.assertEqual(tables1.edges, tables2.edges)
            self.assertEqual(tables1.mutations, tables2.mutations)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestCLI(unittest.TestCase):

    def docmd(self, _cmd):
        cmd = ("-e slim --slim-scaling-factor 20 --slim-burn-in 0 "
               f"{_cmd} -l 0.001 -c chr1 -s 1234 -q 10").split()
        return capture_output(stdpopsim.cli.stdpopsim_main, cmd)

    def test_script_generation(self):
        out, _ = self.docmd("--slim-script HomSap")
        self.assertTrue("sim.registerLateEvent" in out)

        # msprime.MassMigration demographic events, with proportion<1.0
        # low level migration
        out, _ = self.docmd("--slim-script HomSap -d AncientEurasia_9K19")
        self.assertTrue("sim.registerLateEvent" in out)
        # simultaneous mass migrations, with proportions summing to 1.0
        out, _ = self.docmd("--slim-script HomSap -d AmericanAdmixture_4B11")
        self.assertTrue("sim.registerLateEvent" in out)

    def test_simulate(self):
        saved_slim_env = os.environ.get("SLIM")
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"--slim-path slim HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap -o {f.name}")
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_samples, 10)

        # verify sample counts for a multipopulation demographic model
        with tempfile.NamedTemporaryFile(mode="w") as f:
            cmd = ("-e slim --slim-scaling-factor 20 --slim-burn-in 0 "
                   f"HomSap -o {f.name} -l 0.001 -c chr1 -s 1234 -q "
                   "-d OutOfAfrica_3G09 0 0 8").split()
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
            ts = tskit.load(f.name)
        self.assertEqual(ts.num_populations, 3)
        observed_counts = [0, 0, 0]
        for sample in ts.samples():
            observed_counts[ts.get_population(sample)] += 1
        self.assertEqual(observed_counts[0], 0)
        self.assertEqual(observed_counts[1], 0)
        self.assertEqual(observed_counts[2], 8)
        self.assertTrue(all(tree.num_roots == 1 for tree in ts.trees()))

    def test_dry_run(self):
        # --dry-run should run slim, but not create an output file.
        with mock.patch("subprocess.Popen", autospec=True) as mocked_popen:
            # Popen is used as a context manager, so we frob the return value
            # of the returned context manager, rather than the Popen mock itself.
            proc = mocked_popen.return_value.__enter__.return_value
            proc.returncode = 0
            proc.stdout = io.StringIO()
            proc.stderr = io.StringIO()
            with tempfile.NamedTemporaryFile(mode="w") as f:
                self.docmd(f"HomSap --dry-run -o {f.name}")
        mocked_popen.assert_called_once()
        self.assertTrue("slim" in mocked_popen.call_args[0][0])
        with tempfile.NamedTemporaryFile(mode="w") as f:
            self.docmd(f"HomSap --dry-run -o {f.name}")
            self.assertEqual(os.stat(f.name).st_size, 0)

    def test_bad_slim_environ_var(self):
        saved_slim_env = os.environ.get("SLIM")

        os.environ["SLIM"] = "nonexistent"
        with self.assertRaises(FileNotFoundError):
            self.docmd("HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env

    def test_bad_slim_path(self):
        saved_slim_env = os.environ.get("SLIM")

        with self.assertRaises(FileNotFoundError):
            self.docmd("--slim-path nonexistent HomSap")

        if saved_slim_env is None:
            del os.environ["SLIM"]
        else:
            os.environ["SLIM"] = saved_slim_env


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestWarningsAndErrors(unittest.TestCase):
    """
    Checks that warning messages are printed when appropriate.
    """
    def test_odd_sample_warning(self):
        cmd = "-e slim --slim-script HomSap -d OutOfAfrica_2T12 4 6 -q".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 0)

        cmd = "-e slim --slim-script HomSap -d OutOfAfrica_2T12 4 5 -q".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 1)

        cmd = "-e slim --slim-script HomSap -d OutOfAfrica_2T12 3 5 -q".split()
        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            capture_output(stdpopsim.cli.stdpopsim_main, cmd)
        self.assertEqual(mock_warning.call_count, 2)

    def triplet(self):
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22", length_multiplier=0.001)
        return engine, species, contig

    def test_bad_population_size_addSubPop(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

    def test_no_populations_in_generation_1(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100)
        samples = model.get_samples(2)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_bad_population_size_addSubpopSplit(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.IsolationWithMigration(
                NA=1000, N1=100, N2=1000, T=1000, M12=0, M21=0)
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_bad_population_size_setSubpopulationSize(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(100, (1000, 1000))
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=True)

    def test_sample_size_too_big(self):
        engine, species, contig = self.triplet()
        model = stdpopsim.PiecewiseConstantSize(1000)
        samples = model.get_samples(300)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)

    def exp_decline(self, N0=100, N1=1000, T=1000):
        """
        One population model with exponential decline in population size.
        Used for testing that growth rates are handled appropriately.
        """
        r = math.log(N0 / N1) / T
        return stdpopsim.DemographicModel(
                id="exp_decline",
                description="exp_decline",
                long_description="exp_decline",
                populations=[stdpopsim.models._pop0],
                generation_time=1,
                population_configurations=[
                    msprime.PopulationConfiguration(
                        initial_size=N0, growth_rate=r,
                        metadata=stdpopsim.models._pop0.asdict())
                    ],
                demographic_events=[
                    msprime.PopulationParametersChange(
                        time=T, initial_size=N1, growth_rate=0, population_id=0),
                    ],
                )

    def test_bad_population_size_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(2)

        with mock.patch("warnings.warn", autospec=True) as mock_warning:
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)
        mock_warning.assert_called_once()

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=200, dry_run=False)

    def test_sample_size_too_big_exp_decline(self):
        engine, species, contig = self.triplet()
        model = self.exp_decline()
        samples = model.get_samples(30)

        with self.assertRaises(stdpopsim.SLiMException):
            engine.simulate(
                    demographic_model=model, contig=contig, samples=samples,
                    slim_scaling_factor=10, dry_run=True)


class TestSlimAvailable(unittest.TestCase):
    """
    Checks whether SLiM is available or not on platforms that support it.
    """
    def test_parser_has_options(self):
        parser = stdpopsim.cli.stdpopsim_cli_parser()
        with mock.patch("sys.exit", autospec=True):
            _, stderr = capture_output(parser.parse_args, ["--help"])
            # On windows we should have no "slim" options
            self.assertEqual(IS_WINDOWS, "slim" not in stderr)

    def test_engine_available(self):
        all_engines = [engine.id for engine in stdpopsim.all_engines()]
        self.assertEqual(IS_WINDOWS, "slim" not in all_engines)


class PiecewiseConstantSizeMixin(object):
    """
    Mixin that sets up a simple demographic model used by multiple unit tests.
    """
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig("chr22", length_multiplier=0.001)  # ~50 kb

    N0 = 1000  # size in the present
    N1 = 500  # ancestral size
    T = 500  # generations since size change occurred
    T_mut = 300  # introduce a mutation at this generation
    model = stdpopsim.PiecewiseConstantSize(N0, (T, N1))
    model.generation_time = 1
    samples = model.get_samples(100)
    mutation_types = [stdpopsim.ext.MutationType(convert_to_substitution=False)]
    mut_id = len(mutation_types)

    def allele_frequency(self, ts):
        """
        Get the allele frequency of the drawn mutation.
        """
        # surely there's a simpler way!
        assert ts.num_mutations == 1
        alive = list(itertools.chain.from_iterable(
                    ts.individual(i).nodes for i in ts.individuals_alive_at(0)))
        mut = next(ts.mutations())
        tree = ts.at(ts.site(mut.site).position)
        have_mut = [u for u in alive if tree.is_descendant(u, mut.node)]
        af = len(have_mut) / len(alive)
        return af


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestMutationTypes(unittest.TestCase, PiecewiseConstantSizeMixin):
    def test_single_mutation_type_in_script(self):
        engine = stdpopsim.get_engine("slim")
        out, _ = capture_output(
                engine.simulate,
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, slim_script=True)
        self.assertEqual(out.count("initializeMutationType"), 1)

        mutation_types = [stdpopsim.ext.MutationType(weight=1)]
        out, _ = capture_output(
                engine.simulate,
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_script=True)
        self.assertEqual(out.count("initializeMutationType"), 1)

    def test_multiple_mutation_types_in_script(self):
        engine = stdpopsim.get_engine("slim")
        mutation_types = [stdpopsim.ext.MutationType(weight=1),
                          stdpopsim.ext.MutationType(weight=2)]
        out, _ = capture_output(
                engine.simulate,
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_script=True)
        self.assertEqual(out.count("initializeMutationType"), 2)

        mutation_types = stdpopsim.ext.KimDFE()
        positive = stdpopsim.ext.MutationType(convert_to_substitution=False)
        mutation_types.append(positive)
        out, _ = capture_output(
                engine.simulate,
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_script=True)
        self.assertEqual(out.count("initializeMutationType"), 3)

    def test_unweighted_mutations_are_not_simulated_by_slim(self):
        mutation_types = [stdpopsim.ext.MutationType(convert_to_substitution=True),
                          stdpopsim.ext.MutationType(convert_to_substitution=False)]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_scaling_factor=10, slim_burn_in=0.1)
        self.assertEqual(ts.num_sites, 0)

        mutation_types = [stdpopsim.ext.MutationType(),
                          stdpopsim.ext.MutationType(
                              weight=0, distribution_type="g",
                              distribution_args=[-0.01, 0.2])]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_scaling_factor=10, slim_burn_in=0.1)
        self.assertEqual(ts.num_sites, 0)

    def test_weighted_mutations_are_simulated_by_slim(self):
        mutation_types = [stdpopsim.ext.MutationType(
            weight=1, convert_to_substitution=True)]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_scaling_factor=10, slim_burn_in=0.1)
        self.assertTrue(ts.num_sites > 0)

        mutation_types = [stdpopsim.ext.MutationType(
            weight=1, convert_to_substitution=False)]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_scaling_factor=10, slim_burn_in=0.1)
        self.assertTrue(ts.num_sites > 0)

        mutation_types = stdpopsim.ext.KimDFE()
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=mutation_types,
                slim_scaling_factor=10, slim_burn_in=0.1)
        self.assertTrue(ts.num_sites > 0)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestAlleleFrequencyConditioning(unittest.TestCase, PiecewiseConstantSizeMixin):

    def test_drawn_mutation_not_lost(self):
        extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut, mutation_type_id=self.mut_id,
                    population_id=0, coordinate=100, save=True),
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0, end_time=0,
                    mutation_type_id=self.mut_id, population_id=0,
                    op=">", allele_frequency=0),
                ]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=self.mutation_types,
                extended_events=extended_events, slim_scaling_factor=10,
                slim_burn_in=0.1)
        self.assertEqual(ts.num_mutations, 1)

    def test_drawn_mutation_is_lost(self):
        extended_events = [
                stdpopsim.ext.DrawMutation(
                    time=self.T_mut, mutation_type_id=self.mut_id,
                    population_id=0, coordinate=100, save=True),
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0, end_time=0,
                    mutation_type_id=self.mut_id, population_id=0,
                    op="<=", allele_frequency=0),
                ]
        ts = slim_simulate_no_recap(
                demographic_model=self.model, contig=self.contig,
                samples=self.samples, mutation_types=self.mutation_types,
                extended_events=extended_events, slim_scaling_factor=10,
                slim_burn_in=0.1)
        self.assertEqual(ts.num_mutations, 0)

    def test_drawn_mutation_meets_AF_threshold(self):
        for af_threshold, seed in zip((0.01, 0.1, 0.2), (1, 2, 3)):
            extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut, mutation_type_id=self.mut_id,
                        population_id=0, coordinate=100, save=True),
                    # Condition on desired AF at end of simulation.
                    stdpopsim.ext.ConditionOnAlleleFrequency(
                        start_time=0, end_time=0,
                        mutation_type_id=self.mut_id, population_id=0,
                        op=">=", allele_frequency=af_threshold),
                    ]
            ts = slim_simulate_no_recap(
                    demographic_model=self.model, contig=self.contig,
                    samples=self.samples, mutation_types=self.mutation_types,
                    extended_events=extended_events, slim_scaling_factor=10,
                    slim_burn_in=0.1, seed=seed)
            self.assertEqual(ts.num_mutations, 1)
            self.assertTrue(self.allele_frequency(ts) >= af_threshold)

    def test_bad_AF_conditioning_parameters(self):
        for op, af in [
                # bad ops
                ("=", .5), ("==", 0.5), ("!=", 0.5), ({}, 0.5), ("", 0.5),
                # bad allele frequencies
                ("<", -1), ("<=", 2.0),
                # bad combinations
                ("<", 0), (">", 1)]:
            with self.assertRaises(ValueError):
                stdpopsim.ext.ConditionOnAlleleFrequency(
                    start_time=0, end_time=0, mutation_type_id=self.mut_id,
                    population_id=0, op=op, allele_frequency=af)

    def test_op_id(self):
        op_types = stdpopsim.ext.ConditionOnAlleleFrequency.op_types
        for op in op_types:
            id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id(op)
            self.assertTrue(0 <= id <= len(op_types))
        for op in ("==", "=", "!=", {}, ""):
            id = stdpopsim.ext.ConditionOnAlleleFrequency.op_id(op)
            self.assertEqual(id, -1)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestChangeMutationFitness(unittest.TestCase, PiecewiseConstantSizeMixin):
    # Testing stdpopsim.ext.ChangeMutationFitness is challenging, because
    # the side-effects are not deterministic. But if we condition on fixation
    # of a drawn mutation, such a simulation will be very slow without strong
    # positive selection (because we effectively do rejection sampling on the
    # simulation until we get one that meets the allele frequency condition).
    # So if this test takes more than a few seconds to run, that's a good
    # indication that selection is not acting.
    def test_positive_mutation_meets_AF_threshold(self):
        for af_threshold, seed in zip((0.5, 1), (1, 2)):
            extended_events = [
                    stdpopsim.ext.DrawMutation(
                        time=self.T_mut, mutation_type_id=self.mut_id,
                        population_id=0, coordinate=100, save=True),
                    stdpopsim.ext.ChangeMutationFitness(
                        start_time=stdpopsim.ext.GenerationAfter(self.T_mut),
                        end_time=0,
                        mutation_type_id=self.mut_id, population_id=0,
                        selection_coeff=0.1, dominance_coeff=0.5),
                    # Condition on AF > 0, to restore() immediately if the
                    # allele is lost.
                    stdpopsim.ext.ConditionOnAlleleFrequency(
                        start_time=0, end_time=0,
                        mutation_type_id=self.mut_id, population_id=0,
                        op=">", allele_frequency=0),
                    # Condition on desired AF at end of simulation.
                    stdpopsim.ext.ConditionOnAlleleFrequency(
                        start_time=0, end_time=0,
                        mutation_type_id=self.mut_id, population_id=0,
                        op=">=", allele_frequency=af_threshold),
                    ]
            ts = slim_simulate_no_recap(
                    demographic_model=self.model, contig=self.contig,
                    samples=self.samples, mutation_types=self.mutation_types,
                    extended_events=extended_events, slim_scaling_factor=10,
                    slim_burn_in=0.1, seed=seed, verbosity=2)
            self.assertEqual(ts.num_mutations, 1)
            self.assertTrue(self.allele_frequency(ts) >= af_threshold)


@unittest.skipIf(IS_WINDOWS, "SLiM not available on windows")
class TestExtendedEvents(unittest.TestCase, PiecewiseConstantSizeMixin):

    def test_bad_extended_events(self):
        engine = stdpopsim.get_engine("slim")
        for bad_ee in [
                msprime.PopulationParametersChange(time=0, initial_size=100),
                None, {}, "",
                ]:
            with self.assertRaises(ValueError):
                engine.simulate(
                        demographic_model=self.model, contig=self.contig,
                        samples=self.samples, extended_events=[bad_ee],
                        dry_run=True)
