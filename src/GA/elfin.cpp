#include <string>
#include <regex>

#include "util.h"
#include "input/SpecParser.hpp"
#include "input/CSVParser.hpp"
#include "input/JSONParser.hpp"
#include "core/EvolutionSolver.hpp"

class Options
{
public:
    Options() {};
    virtual ~Options() {
    };
    std::string xDBFile = "./xDB.json";
    std::string inputFile = "";
    enum InputType { Unknown, CSV, JSON };
    InputType inputType = Unknown;
    std::string settingsFile = "./settings.json";
    std::string outputFile = "./output.json";
    // use signed primitives because we want to check validity
    long popSize = 10000;
    long nIters = 1000;
};

Options options;

DECL_ARG_CALLBACK(helpAndExit); // defined later due to need of bundle size
DECL_ARG_CALLBACK(setSettingsFile) { options.settingsFile = arg_in; }
DECL_ARG_CALLBACK(setInputFile) { options.inputFile = arg_in; }
DECL_ARG_CALLBACK(setXDB) { options.xDBFile = arg_in; }
DECL_ARG_CALLBACK(setOutput) { options.outputFile = arg_in; }
DECL_ARG_CALLBACK(setPopSize) { options.popSize = parse_long(arg_in); }
DECL_ARG_CALLBACK(setNIters) { options.nIters = parse_long(arg_in); }
DECL_ARG_CALLBACK(setLogLevel) { set_log_level((Log_Level) parse_long(arg_in)); }

const argument_bundle argb[] = {
    {"-h", "--help", "Print this help text and exit", false, helpAndExit},
    {"-s", "--settingsFile", "Set settings file (default ./settings.json)", true, setSettingsFile},
    {"-i", "--inputFile", "Set input file", true, setInputFile},
    {"-x", "--xDBFile", "Set xDB file (default ./xDB.json)", true, setXDB},
    {"-o", "--outputFile", "Set output file (default ./output.json)", true, setOutput},
    {"-p", "--popSize", "Set population size (default 10000)", true, setPopSize},
    {"-n", "--nIters", "Set number of iterations (default 1000)", true, setNIters},
    {"-lg", "--logLevel", "Set log level", true, setLogLevel}
};
#define ARG_BUND_SIZE (sizeof(argb) / sizeof(argb[0]))

DECL_ARG_CALLBACK(helpAndExit)
{
    raw("elfin: Repeat Modular Protein Designer GA\n");
    raw("Usage: ./elfin [OPTIONS]\n");
    raw("Note: Setting an option will override the same setting in the json file\n");

    const argument_bundle * ptr = argb;
    print_arg_bundles(&ptr, 1);
    print_arg_bundles(&ptr, ARG_BUND_SIZE - 1);

    exit(1);
}

DECL_ARG_IN_FAIL_CALLBACK(argParseFail)
{
    printf("Argument parsing failed on string: \"%s\"\n", arg_in);
    helpAndExit("");
    exit(1);
}

void parseSettings()
{
    elfin::JSON j = elfin::JSONParser().parse(options.settingsFile);

    if (j["inputFile"] != NULL)
        options.inputFile = j["inputFile"].get<std::string>();

    if (j["xDBFile"] != NULL)
        options.xDBFile = j["xDBFile"].get<std::string>();

    if (j["outputFile"] != NULL)
        options.outputFile = j["inputFile"].get<std::string>();

    // parse as signed primitives because we want to check validity
    if (j["popSize"] != NULL)
        options.popSize = j["popSize"].get<long>();

    if (j["nIters"] != NULL)
        options.nIters = j["nIters"].get<long>();
}

void checkOptions()
{
    // Do basic checks for each option

    // Files
    panicIf(options.xDBFile == "",
            "No xDB file given. Check your settings.json\n");

    panicIf(options.inputFile == "",
            "No input spec file given. Check help using -h\n");

    panicIf(options.settingsFile == "",
            "No settings file file given. Check help using -h\n");

    // Extensions
    if (std::regex_match(
                options.inputFile,
                std::regex("(.*)(\\.csv)$", std::regex::icase)))
    {
        msg("Using CSV input\n");
        options.inputType = Options::InputType::CSV;
    }
    else if (std::regex_match(
                 options.inputFile,
                 std::regex("(.*)(\\.json$)", std::regex::icase)))
    {
        msg("Using JSON input\n");
        options.inputType = Options::InputType::JSON;
    }
    else {
        die("Unrecognized input file type\n");
    }

    // Settings

    panicIf(options.popSize < 0, "Population size cannot be < 0\n");

    panicIf(options.nIters < 0, "Number of iterations cannot be < 0\n");

}

elfin::Points3f parseInput()
{
    using namespace elfin;
    switch (options.inputType)
    {
    case Options::InputType::CSV:
        return CSVParser().parseSpec(options.inputFile);
    case Options::InputType::JSON:
        return JSONParser().parseSpec(options.inputFile);
    default:
        die("Unknown input format\n");
    }
}

/*
 * The elfin design process:
 *   Input:
 *      A vector of Centre of Mass shape specification
 *   Algorithm:
 *      GA with a variety of inheritance and also a desturctive
 *      gene (shape candidate) operator
 *   Output:
 *      A vector of module (node) names suitable for
 *      use by Synth.py to produce full PDB
 */
int main(int argc, const char ** argv)
{
    // Get values from settings json first
    parseSettings();

    // Parse user arguments
    parse_args(argc, argv, ARG_BUND_SIZE, argb, argParseFail);

    checkOptions();

    using namespace elfin;
    RelaMat rm;
    NameIdMap nim;
    JSONParser().parseDB(options.xDBFile, rm , nim);

    Points3f spec = parseInput();

    EvolutionSolver es = EvolutionSolver(rm, spec);
    es.run(options.popSize, options.nIters);

    wrn("TODO: Output a json containing \"nodes\" solution data, and score\n");

    return 0;
}