"""
The command line interface for stdpopsim. Provides standard simulations
at the command line and methods to manage resources used by stdpopsim.
"""
import argparse
import json
import logging
import platform
import sys
import textwrap
import tempfile
import pathlib
import shutil
import functools

import msprime
import tskit
import humanize

import stdpopsim

# resource is from the standard library, but it's not available on
# Windows. We break from the usual import grouping conventions here
# to avoid lots of pep8 complaints about mixing imports and code.
_resource_module_available = False
try:
    import resource
    _resource_module_available = True
except ImportError:
    pass


logger = logging.getLogger(__name__)


def exit(message):
    """
    Exit with the specified error message, setting error status.
    """
    sys.exit(f"{sys.argv[0]}: {message}")


LOG_FORMAT = "%(asctime)s [%(process)d] %(levelname)s %(name)s: %(message)s"


def setup_logging(args):
    log_level = "WARN"
    if args.verbosity > 0:
        log_level = "INFO"
    if args.verbosity > 1:
        log_level = "DEBUG"
    logging.basicConfig(format=LOG_FORMAT, level=log_level)


def get_species_wrapper(species_id):
    try:
        return stdpopsim.get_species(species_id)
    except ValueError as ve:
        exit(str(ve))


def get_model_wrapper(species, model_id):
    try:
        return species.get_demographic_model(model_id)
    except ValueError as ve:
        exit(str(ve))


def get_genetic_map_wrapper(species, genetic_map_id):
    try:
        return species.get_genetic_map(genetic_map_id)
    except ValueError as ve:
        exit(str(ve))


def get_models_help(species_id, model_id):
    """
    Generate help text for the specified species. If model_id is None, generate
    help for all models. Otherwise, it must be a string with a valid model ID.
    """
    species = stdpopsim.get_species(species_id)
    if model_id is None:
        models_text = f"\nAll simulation models for {species.name}\n\n"
        models = [model.id for model in species.demographic_models]
    else:
        models = [model_id]
        models_text = f"\nModel description\n\n"

    # TODO improve this text formatting.
    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for model_id in models:
        model = get_model_wrapper(species, model_id)
        models_text += f"{model.id}: {model.name}\n"
        models_text += wrapper.fill(textwrap.dedent(model.description))
        models_text += "\n\n"

        models_text += indent + "Populations:\n"

        for population in model.populations:
            if population.allow_samples:
                models_text += indent * 2
                models_text += f"{population.name}: {population.description}\n"
        models_text += "\n"

    return models_text


class HelpModels(argparse.Action):
    """
    Action used to produce model help text.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_models_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


def get_genetic_maps_help(species_id, genetic_map_id):
    """
    Generate help text for the given genetic map. If map_id is None, generate
    help for all genetic maps. Otherwise, it must be a string with a valid map
    ID.
    """
    species = stdpopsim.get_species(species_id)
    if genetic_map_id is None:
        maps_text = f"\nAll genetic maps for {species.name}\n\n"
        maps = [genetic_map.name for genetic_map in species.genetic_maps]
    else:
        maps = [genetic_map_id]
        maps_text = f"\nGenetic map description\n\n"

    indent = " " * 4
    wrapper = textwrap.TextWrapper(initial_indent=indent, subsequent_indent=indent)
    for map_id in maps:
        map = get_genetic_map_wrapper(species, map_id)
        maps_text += f"{map.name}\n"
        maps_text += wrapper.fill(textwrap.dedent(map.description))
        maps_text += "\n\n"

    return maps_text


class HelpGeneticMaps(argparse.Action):
    """
    Action used to produce genetic map help text.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        help_text = get_genetic_maps_help(namespace.species, values)
        print(help_text, file=sys.stderr)
        parser.exit()


def get_species_help(species_id):
    """
    Generate help text for the given species with some of the species attributes
    that are not covered by the other helps.
    """
    species = stdpopsim.get_species(species_id)
    species_text = f"\nDefault population parameters for {species.name}:\n"

    species_text += f"Generation time: {species.generation_time}\n"
    species_text += f"Population size: {species.population_size}\n"
    species_text += f"Mutation rate: {species.genome.mean_mutation_rate}\n"
    species_text += f"Recombination rate: {species.genome.mean_recombination_rate}\n"
    return species_text


def get_environment():
    """
    Returns a dictionary describing the environment in which stdpopsim
    is currently running.
    """
    env = {
        "os": {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version(),
        },
        "libraries": {
            "msprime": {"version": msprime.__version__},
            "tskit": {"version": tskit.__version__},
        }
    }
    return env


def get_provenance_dict():
    """
    Returns a dictionary encoding an execution of stdpopsim conforming to the
    tskit provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "stdpopsim",
            "version": stdpopsim.__version__
        },
        "parameters": {
            "command": sys.argv[0],
            "args": sys.argv[1:]
        },
        "environment": get_environment()
    }
    return document


def write_output(ts, args):
    """
    Adds provenance information to the specified tree sequence (ensuring that the
    output is reproducible) and write the resulting tree sequence to output.
    """
    tables = ts.dump_tables()
    logger.debug("Updating provenance")
    provenance = get_provenance_dict()
    tables.provenances.add_row(json.dumps(provenance))
    ts = tables.tree_sequence()
    if args.output is None:
        # There's no way to get tskit to write directly to stdout, so we write
        # to a tempfile first.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = pathlib.Path(tmpdir) / "tmp.trees"
            ts.dump(tmpfile)
            with open(tmpfile, "rb") as f:
                shutil.copyfileobj(f, sys.stdout.buffer)
    else:
        logger.info(f"Writing to {args.output}")
        ts.dump(args.output)


def write_bibtex(engine, model, contig, bibtex_file):
    """
    Write bibtex for available citations to a file."""
    for citation in engine.citations:
        bibtex_file.write(citation.fetch_bibtex())
    if contig.genetic_map is not None:
        for citation in contig.genetic_map.citations:
            bibtex_file.write(citation.fetch_bibtex())
    for citation in model.citations:
        bibtex_file.write(citation.fetch_bibtex())
    bibtex_file.close()


def write_citations(engine, model, contig):
    """
    Write out citation information so that the user knows what papers to cite
    for the simulation engine, the model and the mutation/recombination rate
    information.
    """
    printerr = functools.partial(print, file=sys.stderr)
    printerr("If you use this simulation in published work, please cite:")
    printerr(f"Simulation engine: {engine.name}")
    for citation in engine.citations:
        printerr(f"  {citation}")
    if contig.genetic_map is not None:
        printerr(f"Genetic map: {contig.genetic_map.name}")
        for citation in contig.genetic_map.citations:
            printerr(f"  {citation}")
    if len(model.citations) > 0:
        printerr(f"Simulation model: {model.name}")
        for citation in model.citations:
            printerr(f"  {citation}")


def summarise_usage():
    # Don't report usage on Windows as the resource module is not available.
    #  We could do this using the psutil external library, if demand exists.
    if _resource_module_available:
        rusage = resource.getrusage(resource.RUSAGE_SELF)
        user_time = humanize.naturaldelta(rusage.ru_utime)
        sys_time = rusage.ru_stime
        max_mem = rusage.ru_maxrss
        if sys.platform != 'darwin':
            max_mem *= 1024  # Linux and other OSs (e.g. freeBSD) report maxrss in kb
        max_mem_str = humanize.naturalsize(max_mem, binary=True)
        logger.info("rusage: user={}; sys={:.2f}s; max_rss={}".format(
            user_time, sys_time, max_mem_str))


def add_simulate_species_parser(parser, species):
    header = (
        f"Run simulations for {species.name} using up-to-date genome information, "
        "genetic maps and simulation models from the literature. "
        "NOTE: By default, the tskit '.trees' binary file is written to stdout,"
        "so you should either redirect this to a file or use the '--output' "
        "option to specify a filename."
    )

    description_text = textwrap.fill(header) + "\n" + get_species_help(species.id)

    species_parser = parser.add_parser(
        f"{species.id}",
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=f"Run simulations for {species.name}.")
    species_parser.set_defaults(species=species.id)
    species_parser.set_defaults(genetic_map=None)
    species_parser.set_defaults(chromosome=None)
    species_parser.add_argument(
        "--help-models", action=HelpModels, nargs="?",
        help=(
            "Print descriptions of simulation models and exit. If a model ID "
            "is provided as an argument show help for this model; otherwise "
            "show help for all available models"))
    species_parser.add_argument(
        "-b", "--bibtex-file",
        type=argparse.FileType('w'),
        help="Write citations to a given bib file. This will overwrite the file.",
        default=None,
        action='store')

    # Set metavar="" to prevent help text from writing out the explicit list
    # of options, which can be too long and ugly.
    choices = [gm.name for gm in species.genetic_maps]
    if len(species.genetic_maps) > 0:
        species_parser.add_argument(
            "--help-genetic-maps", action=HelpGeneticMaps, nargs="?",
            help=(
                "Print list of genetic maps and exit. If a genetic map ID is "
                "given as an argument, show help for this map. Otherwise show "
                "help for all available genetic maps"))

    species_parser.add_argument(
        "-q", "--quiet", action='store_true',
        help="Do not write out citation information")
    species_parser.add_argument(
        "-D", "--dry-run", action='store_true', default=False,
        help="Do not run actual simulation")

    if len(species.genetic_maps) > 0:
        species_parser.add_argument(
            "-g", "--genetic-map",
            choices=choices, metavar="", default=None,
            help=(
                "Specify a particular genetic map. By default, a chromosome-specific "
                "uniform recombination rate is used. These default rates are listed in "
                "the catalog: <https://stdpopsim.readthedocs.io/en/latest/catalog.html> "
                "Available maps: "
                f"{', '.join(choices)}. "))

    if len(species.genome.chromosomes) > 1:
        choices = [chrom.id for chrom in species.genome.chromosomes]
        species_parser.add_argument(
            "-c", "--chromosome", choices=choices, metavar="", default=choices[0],
            help=(
                f"Simulate a specific chromosome. "
                f"Options: {', '.join(choices)}. "
                f"Default={choices[0]}."))
    species_parser.add_argument(
        "-l", "--length-multiplier", default=1, type=float,
        help="Simulate a sequence of length l times the named chromosome's length, "
             "using the named chromosome's mutation and recombination rates.")
    species_parser.add_argument(
        "-s", "--seed", default=None, type=int,
        help=(
            "The random seed to use for simulations. If not specified a "
            "high-quality random seed will be generated automatically. "
            "For msprime, seeds must be > 0 and < 2^32."))

    model_help = (
        "Specify a simulation model. If no model is specified, a single population"
        "constant size model is used. Available models:"
        f"{', '.join(model.id for model in species.demographic_models)}"
        ". Please see --help-models for details of these models.")
    species_parser.add_argument(
        "-m", "--model", default=None, metavar="",
        choices=[model.id for model in species.demographic_models],
        help=model_help)
    species_parser.add_argument(
        "-o", "--output",
        help=(
            "Where to write the output tree sequence file. Defaults to "
            "stdout if not specified"))

    species_parser.add_argument(
        "samples", type=int, nargs="+",
        help=(
            "The number of samples to draw from each population. At least "
            "two samples must be specified. The number of arguments that "
            "will be accepted depends on the simulation model that is "
            "specified: for a model that has n populations, we can specify "
            "the number of samples to draw from each of these populations."
            "We do not need to provide sample numbers of each of the "
            "populations; those that are omitted are set to zero."))

    def run_simulation(args):
        if args.model is None:
            model = stdpopsim.PiecewiseConstantSize(species.population_size)
            model.generation_time = species.generation_time
            model.citations.extend(species.population_size_citations)
            model.citations.extend(species.generation_time_citations)
        else:
            model = get_model_wrapper(species, args.model)
        if len(args.samples) > model.num_sampling_populations:
            exit(
                f"Cannot sample from more than {model.num_sampling_populations} "
                "populations")
        samples = model.get_samples(*args.samples)

        contig = species.get_contig(
            args.chromosome, genetic_map=args.genetic_map,
            length_multiplier=args.length_multiplier)
        engine = stdpopsim.get_engine(args.engine)
        logger.info(
            f"Running simulation model {model.name} for {species.name} on "
            f"{contig} with {len(samples)} samples using {engine.name}.")

        kwargs = vars(args)
        kwargs.update(demographic_model=model, contig=contig, samples=samples)
        if not args.quiet:
            write_simulation_summary(engine=engine, model=model, contig=contig,
                                     samples=samples, seed=args.seed)
        if not args.dry_run:
            ts = engine.simulate(**kwargs)
            summarise_usage()
            write_output(ts, args)
        if not args.quiet:
            write_citations(engine, model, contig)
        if args.bibtex_file is not None:
            write_bibtex(engine, model, contig, args.bibtex_file)

    species_parser.set_defaults(runner=run_simulation)


def write_simulation_summary(engine, model, contig, samples, seed=None):
    # Setup writer
    printerr = functools.partial(print, file=sys.stderr)
    indent = " " * 4
    # Header
    dry_run_text = "Simulation information:\n"
    # Get information about the engine
    dry_run_text += f"{indent}Engine: {engine.id} ({engine.get_version()})\n"
    # Get information about model
    dry_run_text += f"{indent}Model id: {model.id}\n"
    dry_run_text += f"{indent}Model name: {model.name}\n"
    # Add seed information if extant
    if seed is not None:
        dry_run_text += f"{indent}Seed: {seed}\n"
    # Get information about the number of samples
    dry_run_text += f"{indent}Population: number_samples (sampling_time_generations):\n"
    sample_counts = [0] * model.num_sampling_populations
    for s in samples:
        sample_counts[s.population] += 1
    for p in range(0, model.num_sampling_populations):
        pop_name = model.populations[p].name
        sample_time = model.populations[p].sampling_time
        dry_run_text += f"{indent * 2}{pop_name}: "
        dry_run_text += f"{sample_counts[p]} ({sample_time})\n"
    # Get information about relevant contig
    gmap = "None" if contig.genetic_map is None else contig.genetic_map.name
    mean_recomb_rate = contig.recombination_map.mean_recombination_rate
    mut_rate = contig.mutation_rate
    contig_len = contig.recombination_map.get_length()
    dry_run_text += f"Contig Description:\n"
    dry_run_text += f"{indent}Contig length: {contig_len}\n"
    dry_run_text += f"{indent}Mean recombination rate: {mean_recomb_rate}\n"
    dry_run_text += f"{indent}Mean mutation rate: {mut_rate}\n"
    dry_run_text += f"{indent}Genetic map: {gmap}\n"
    printerr(dry_run_text)


def run_download_genetic_maps(args):
    species_names = [args.species]
    if args.species is None:
        species_names = [species.id for species in stdpopsim.all_species()]
    for species_id in species_names:
        species = get_species_wrapper(species_id)
        if len(args.genetic_maps) == 0:
            genetic_maps = [gmap.name for gmap in species.genetic_maps]
        else:
            genetic_maps = args.genetic_maps
        for genetic_map_id in genetic_maps:
            genetic_map = get_genetic_map_wrapper(species, genetic_map_id)
            genetic_map.download()


def stdpopsim_cli_parser():

    # TODO the CLI defined by this hierarchical and clumsy, but it's the best
    # I could figure out. It can definitely be improved!
    top_parser = argparse.ArgumentParser(
        description="Command line interface for stdpopsim.")
    top_parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(stdpopsim.__version__))
    top_parser.add_argument(
        "-v", "--verbosity", action='count', default=0,
        help="Increase the verbosity")
    top_parser.add_argument(
        "-c", "--cache-dir", type=str, default=None,
        help=(
            "Set the cache directory to the specified value. "
            "Note that this can also be set using the environment variable "
            "STDPOPSIM_CACHE. If both the environment variable and this "
            "option are set, the option takes precedence. "
            f"Default: {stdpopsim.get_cache_dir()}"))

    top_parser.add_argument(
        "-e", "--engine",
        default=stdpopsim.get_default_engine().id,
        choices=[e.id for e in stdpopsim.all_engines()],
        help="Specify a simulation engine.")

    for engine in stdpopsim.all_engines():
        group = top_parser.add_argument_group(
                f"{engine.name} specific parameters")
        engine.add_arguments(group)

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    for species in stdpopsim.all_species():
        add_simulate_species_parser(subparsers, species)

    download_maps_parser = subparsers.add_parser(
        "download-genetic-maps",
        help="Download genetic maps",
        description=(
            "Download genetic maps and store them in the cache directory. "
            "Maps are downloaded regardless of whether they are already "
            "in the cache or not. Please use the --cache-dir option to "
            "download maps to a specific directory. "))
    download_maps_parser.add_argument(
        "species", nargs="?",
        help=(
            "Download genetic maps for this species. If not specified "
            "download all known genetic maps."))
    download_maps_parser.add_argument(
        "genetic_maps", type=str, nargs="*",
        help=(
            "If specified, download these genetic maps. If no maps "
            "are provided, download all maps for this species."))

    download_maps_parser.set_defaults(runner=run_download_genetic_maps)

    return top_parser


# This function only exists to make mocking out the actual running of
# the program easier.
def run(args):
    args.runner(args)


def stdpopsim_main(arg_list=None):
    parser = stdpopsim_cli_parser()
    args = parser.parse_args(arg_list)
    setup_logging(args)
    if args.cache_dir is not None:
        stdpopsim.set_cache_dir(args.cache_dir)
    run(args)
