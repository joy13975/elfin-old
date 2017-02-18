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
    std::string xDBFile = "";
    std::string inputFile = "";
    enum InputType { Unknown, CSV, JSON };
    InputType inputType = Unknown;
    std::string settingsFile = "./settings.json";
    std::string outputFile = "./output.json";
};

Options options;

DECL_ARG_CALLBACK(setInputFile) { options.inputFile = arg_in; }
DECL_ARG_CALLBACK(setSettingsFile) { options.settingsFile = arg_in; }
DECL_ARG_CALLBACK(helpAndExit); // defined later due to need of bundle size
DECL_ARG_CALLBACK(setOutput) { options.outputFile = arg_in; }
DECL_ARG_CALLBACK(setLogLevel) { set_log_level((Log_Level) parse_long(arg_in)); }

const argument_bundle argb[] = {
    // necessary argument
    {"-i", "--input <file>", "Set input file", true, setInputFile},
    {"-s", "--settings <file>", "Set settings file, default ./settings.json", true, setSettingsFile},
    // optional argument
    {"-h", "--help", "Print this help text and exit", false, helpAndExit},
    {"-lg", "--loglevel <0-4>", "Set log level", true, setLogLevel},
    {"-o", "--output <file>", "Set output file, default ./output.json", true, setOutput}
};
#define ARG_BUND_SIZE (sizeof(argb) / sizeof(argb[0]))

DECL_ARG_CALLBACK(helpAndExit)
{
    raw("elfin: Repeat Modular Protein Designer GA\n");
    raw("Usage: ./elfin [INPUT] [OPTIONS]\n");
    raw("Note: Setting an option will override the same setting in the json file\n");

    const argument_bundle * ptr = argb;
    print_arg_title("INPUT:");
    print_arg_bundles(&ptr, 1);
    print_arg_title("OPTIONS:");
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
        options.inputFile = j["inputFile"].get<std::string>();
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
}

std::vector<float> parseInput()
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

    std::vector<float> specComs = parseInput();

    EvolutionSolver es = EvolutionSolver(rm, specComs);
    es.run();

    return 0;
}