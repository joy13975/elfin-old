#include <string>
#include <regex>
#include <sstream>

#include "data/TypeDefs.hpp"
#include "util.h"
#include "input/SpecParser.hpp"
#include "input/CSVParser.hpp"
#include "input/JSONParser.hpp"
#include "core/EvolutionSolver.hpp"

namespace elfin
{

static OptionPack options;

DECL_ARG_CALLBACK(helpAndExit); // defined later due to need of bundle size
DECL_ARG_CALLBACK(setSettingsFile) { options.settingsFile = arg_in; }
DECL_ARG_CALLBACK(setInputFile) { options.inputFile = arg_in; }
DECL_ARG_CALLBACK(setXDB) { options.xDBFile = arg_in; }
DECL_ARG_CALLBACK(setOutput) { options.outputFile = arg_in; }

DECL_ARG_CALLBACK(setChromoLenDev) { options.chromoLenDev = parse_float(arg_in); }
DECL_ARG_CALLBACK(setAvgPairDist) { options.avgPairDist = parse_float(arg_in); }
DECL_ARG_CALLBACK(setRandSeed) { options.randSeed = parse_long(arg_in); }

DECL_ARG_CALLBACK(setGaPopSize) { options.gaPopSize = parse_long(arg_in); }
DECL_ARG_CALLBACK(setGaIters) { options.gaIters = parse_long(arg_in); }
DECL_ARG_CALLBACK(setGaSurviveRate) { options.gaSurviveRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setGaCrossRate) { options.gaCrossRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setGaMutateRate) { options.gaMutateRate = parse_float(arg_in); }

DECL_ARG_CALLBACK(setLogLevel) { set_log_level((Log_Level) parse_long(arg_in)); }

const argument_bundle argb[] = {
    {"-h", "--help", "Print this help text and exit", false, helpAndExit},
    {"-s", "--settingsFile", "Set settings file (default ./settings.json)", true, setSettingsFile},
    {"-i", "--inputFile", "Set input file", true, setInputFile},
    {"-x", "--xDBFile", "Set xDB file (default ./xDB.json)", true, setXDB},
    {"-o", "--outputFile", "Set output file (default ./output.json)", true, setOutput},
    {"-d", "--chromoLenDev", "Set chromosome length deviation allowance (default 0.20)", true, setChromoLenDev},
    {"-a", "--avgPairDist", "Overwrite default average distance between pairs of CoMs (default 38.0)", true, setAvgPairDist},
    {"-rs", "--randSeed", "Set RNG seed (default 0x1337cafe; setting to 0 uses time as seed)", true, setRandSeed},
    {"-gps", "--gaPopSize", "Set GA population size (default 10000)", true, setGaPopSize},
    {"-git", "--gaIters", "Set GA iterations (default 1000)", true, setGaIters},
    {"-gsr", "--gaSurviveRate", "Set GA survival rate (default 0.1)", true, setGaSurviveRate},
    {"-gcr", "--gaCrossRate", "Set GA surviver cross rate (default 0.60)", true, setGaCrossRate},
    {"-gmr", "--gaMutateRate", "Set GA surviver mutation rate (default 0.3)", true, setGaMutateRate},
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
    JSON j = JSONParser().parse(options.settingsFile);

    if (j["inputFile"] != NULL)
        options.inputFile = j["inputFile"].get<std::string>();

    if (j["xDBFile"] != NULL)
        options.xDBFile = j["xDBFile"].get<std::string>();

    if (j["outputFile"] != NULL)
        options.outputFile = j["inputFile"].get<std::string>();

    // Parse as signed primitives for validation
    if (j["chromoLenDev"] != NULL)
        options.chromoLenDev = j["chromoLenDev"].get<float>();

    if (j["randSeed"] != NULL)
    {
        std::stringstream ss;
        ss << std::hex << j["randSeed"].get<std::string>();
        ss >> options.randSeed;
    }

    // GA params
    if (j["gaPopSize"] != NULL)
        options.gaPopSize = j["gaPopSize"].get<long>();

    if (j["gaIters"] != NULL)
        options.gaIters = j["gaIters"].get<long>();

    if (j["gaSurviveRate"] != NULL)
        options.gaSurviveRate = j["gaSurviveRate"].get<float>();

    if (j["gaCrossRate"] != NULL)
        options.gaCrossRate = j["gaCrossRate"].get<float>();

    if (j["gaMutateRate"] != NULL)
        options.gaMutateRate = j["gaMutateRate"].get<float>();


    if (j["avgPairDist"] != NULL)
        options.avgPairDist = j["avgPairDist"].get<float>();
}

void checkOptions()
{
    // Do basic checks for each option

    // Files
    panic_if(options.xDBFile == "",
             "No xDB file given. Check your settings.json\n");

    panic_if(options.inputFile == "",
             "No input spec file given. Check help using -h\n");

    panic_if(options.settingsFile == "",
             "No settings file file given. Check help using -h\n");

    // Extensions
    if (std::regex_match(
                options.inputFile,
                std::regex("(.*)(\\.csv)$", std::regex::icase)))
    {
        msg("Using CSV input\n");
        options.inputType = OptionPack::InputType::CSV;
    }
    else if (std::regex_match(
                 options.inputFile,
                 std::regex("(.*)(\\.json$)", std::regex::icase)))
    {
        msg("Using JSON input\n");
        options.inputType = OptionPack::InputType::JSON;
    }
    else {
        die("Unrecognized input file type\n");
    }

    // Settings

    panic_if(options.gaPopSize < 0, "Population size cannot be < 0\n");

    panic_if(options.gaIters < 0, "Number of iterations cannot be < 0\n");

    panic_if(options.chromoLenDev < 0.0 ||
             options.chromoLenDev > 1.0,
             "Gene length deviation must be between 0 and 1 inclusive\n");

    // GA params
    panic_if(options.gaSurviveRate < 0.0 ||
             options.gaSurviveRate > 1.0,
             "GA survive rate must be between 0 and 1 inclusive\n");
    panic_if(options.gaCrossRate < 0.0 ||
             options.gaCrossRate > 1.0,
             "GA cross rate must be between 0 and 1 inclusive\n");
    panic_if(options.gaMutateRate < 0.0 ||
             options.gaMutateRate > 1.0,
             "GA mutate rate must be between 0 and 1 inclusive\n");

    panic_if(options.gaCrossRate +
             options.gaMutateRate > 1.0,
             "Sum of GA cross rate + mutate rate must be <= 1\n");

    panic_if(options.avgPairDist < 0, "Average CoM distance must be > 0\n");

}

Points3f parseInput()
{
    using namespace elfin;
    switch (options.inputType)
    {
    case OptionPack::InputType::CSV:
        return CSVParser().parseSpec(options.inputFile);
    case OptionPack::InputType::JSON:
        return JSONParser().parseSpec(options.inputFile);
    default:
        die("Unknown input format\n");
    }
}

} // namespace elfin

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
    using namespace elfin;

    // Get values from settings json first
    parseSettings();

    // Parse user arguments
    parse_args(argc, argv, ARG_BUND_SIZE, argb, argParseFail);

    checkOptions();

    RelaMat relaMat;
    NameIdMap nameIdMap;
    IdNameMap idNameMap;
    RadiiList radiiList;
    JSONParser().parseDB(options.xDBFile, nameIdMap, idNameMap, relaMat, radiiList);

    Gene::setup(&idNameMap);

    Points3f spec = parseInput();

    EvolutionSolver(relaMat,
                    spec,
                    radiiList,
                    options).run();

    wrn("TODO: Output a json containing \"nodes\" solution data, and score\n");

    return 0;
}