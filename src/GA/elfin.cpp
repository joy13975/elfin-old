#include <string>
#include <regex>
#include <sstream>
#include <csignal>
#include <iostream>
#include <fstream>

#include "data/TypeDefs.hpp"
#include "util.h"
#include "input/SpecParser.hpp"
#include "input/CSVParser.hpp"
#include "input/JSONParser.hpp"
#include "core/EvolutionSolver.hpp"
#include "core/ParallelUtils.hpp"
#include "core/Mathutils.hpp"
#include "core/Kabsch.hpp"

namespace elfin
{

static OptionPack options;

DECL_ARG_CALLBACK(helpAndExit); // defined later due to need of bundle size
DECL_ARG_CALLBACK(setSettingsFile) { options.settingsFile = arg_in; }
DECL_ARG_CALLBACK(setInputFile) { options.inputFile = arg_in; }
DECL_ARG_CALLBACK(setXDB) { options.xDBFile = arg_in; }
DECL_ARG_CALLBACK(setOutputDir) { 
    options.outputDir = arg_in; 
    mkdir_ifn_exists(arg_in);
}

DECL_ARG_CALLBACK(setChromoLenDev) { options.chromoLenDev = parse_float(arg_in); }
DECL_ARG_CALLBACK(setAvgPairDist) { options.avgPairDist = parse_float(arg_in); }
DECL_ARG_CALLBACK(setRandSeed) { options.randSeed = parse_long(arg_in); }

DECL_ARG_CALLBACK(setGaPopSize) { options.gaPopSize = parse_long(arg_in); }
DECL_ARG_CALLBACK(setGaIters) { options.gaIters = parse_long(arg_in); }
DECL_ARG_CALLBACK(setGaSurviveRate) { options.gaSurviveRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setGaCrossRate) { options.gaCrossRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setGaPointMutateRate) { options.gaPointMutateRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setGaLimbMutateRate) { options.gaLimbMutateRate = parse_float(arg_in); }
DECL_ARG_CALLBACK(setScoreStopThreshold) { options.scoreStopThreshold = parse_float(arg_in); }

DECL_ARG_CALLBACK(setLogLevel) { set_log_level((Log_Level) parse_long(arg_in)); }
DECL_ARG_CALLBACK(setRunUnitTests) { options.runUnitTests = true; }

const argument_bundle argb[] = {
    {"-h", "--help", "Print this help text and exit", false, helpAndExit},
    {"-s", "--settingsFile", "Set settings file (default ./settings.json)", true, setSettingsFile},
    {"-i", "--inputFile", "Set input file", true, setInputFile},
    {"-x", "--xDBFile", "Set xDB file (default ./xDB.json)", true, setXDB},
    {"-o", "--outputDir", "Set output directory (default ./out/)", true, setOutputDir},
    {"-d", "--chromoLenDev", "Set chromosome length deviation allowance (default 0.20)", true, setChromoLenDev},
    {"-a", "--avgPairDist", "Overwrite default average distance between pairs of CoMs (default 38.0)", true, setAvgPairDist},
    {"-rs", "--randSeed", "Set RNG seed (default 0x1337cafe; setting to 0 uses time as seed)", true, setRandSeed},
    {"-gps", "--gaPopSize", "Set GA population size (default 10000)", true, setGaPopSize},
    {"-git", "--gaIters", "Set GA iterations (default 1000)", true, setGaIters},
    {"-gsr", "--gaSurviveRate", "Set GA survival rate (default 0.1)", true, setGaSurviveRate},
    {"-gcr", "--gaCrossRate", "Set GA surviver cross rate (default 0.60)", true, setGaCrossRate},
    {"-gmr", "--gaPointMutateRate", "Set GA surviver point mutation rate (default 0.3)", true, setGaPointMutateRate},
    {"-gmr", "--gaLimbMutateRate", "Set GA surviver limb mutation rate (default 0.3)", true, setGaLimbMutateRate},
    {"-stt", "--scoreStopThreshold", "Set GA exit score threshold (default 0.0)", true, setScoreStopThreshold},
    {"-lg", "--logLevel", "Set log level", true, setLogLevel},
    {"-t", "--test", "Run unit tests", false, setRunUnitTests}
};
const size_t ARG_BUND_SIZE = (sizeof(argb) / sizeof(argb[0]));

DECL_ARG_CALLBACK(helpAndExit)
{
    raw("elfin: A Modular Protein Design Tool\n");
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
    panic_if(!file_exists(options.settingsFile.c_str()),
             "Settings file does not exist: \"%s\"\n",
             options.settingsFile.c_str());

    JSON j = JSONParser().parse(options.settingsFile);

    // Perhaps there's a better way to parse these
    // setting-variable names so it never gets out
    // of sync with argb...

    std::ostringstream ss;

    auto jsonToCStr = [&](const JSON & j) {
        ss.str("");

        if (j.is_string())
            ss << j.get<std::string>();
        else
            ss << j;

        return ss.str().c_str();
    };


    if (!j["inputFile"].is_null())
        setInputFile(jsonToCStr(j["inputFile"]));

    if (!j["xDBFile"].is_null())
        setXDB(jsonToCStr(j["xDBFile"]));

    if (!j["outputDir"].is_null())
        setOutputDir(jsonToCStr(j["outputDir"]));

    // Parse as signed primitives for validation
    if (!j["chromoLenDev"].is_null())
        setChromoLenDev(jsonToCStr(j["chromoLenDev"]));

    if (!j["randSeed"].is_null())
        setRandSeed(jsonToCStr(j["randSeed"]));

    // GA params
    if (!j["gaPopSize"].is_null())
        setGaPopSize(jsonToCStr(j["gaPopSize"]));

    if (!j["gaIters"].is_null())
        setGaIters(jsonToCStr(j["gaIters"]));

    if (!j["gaSurviveRate"].is_null())
        setGaSurviveRate(jsonToCStr(j["gaSurviveRate"]));

    if (!j["gaCrossRate"].is_null())
        setGaCrossRate(jsonToCStr(j["gaCrossRate"]));

    if (!j["gaPointMutateRate"].is_null())
        setGaPointMutateRate(jsonToCStr(j["gaPointMutateRate"]));

    if (!j["gaLimbMutateRate"].is_null())
        setGaLimbMutateRate(jsonToCStr(j["gaLimbMutateRate"]));

    if (!j["scoreStopThreshold"].is_null())
        setScoreStopThreshold(jsonToCStr(j["scoreStopThreshold"]));

    if (!j["avgPairDist"].is_null())
        setAvgPairDist(jsonToCStr(j["avgPairDist"]));

}

void checkOptions()
{
    // Do basic checks for each option

    // Files
    panic_if(options.xDBFile == "",
             "No xDB file given. Check your settings.json\n");

    panic_if(!file_exists(options.xDBFile.c_str()),
             "xDB file could not be found\n");

    panic_if(options.inputFile == "",
             "No input spec file given. Check help using -h\n");

    panic_if(!file_exists(options.inputFile.c_str()),
             "Input file could not be found\n");

    panic_if(options.settingsFile == "",
             "No settings file file given. Check help using -h\n");

    panic_if(!file_exists(options.settingsFile.c_str()),
             "Settings file \"%s\" could not be found\n",
             options.settingsFile.c_str());

    panic_if(options.outputDir == "",
             "No output directory given. Check help using -h\n");

    panic_if(!file_exists(options.outputDir.c_str()),
             "Output directory could not be found\n");

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
    panic_if(options.gaPointMutateRate < 0.0 ||
             options.gaPointMutateRate > 1.0,
             "GA point mutate rate must be between 0 and 1 inclusive\n");
    panic_if(options.gaLimbMutateRate < 0.0 ||
             options.gaLimbMutateRate > 1.0,
             "GA limb mutate rate must be between 0 and 1 inclusive\n");

    panic_if(options.gaCrossRate +
             options.gaPointMutateRate +
             options.gaLimbMutateRate > 1.0,
             "Sum of GA cross + point mutate + limb mutate rates must be <= 1\n");

    panic_if(options.avgPairDist < 0, "Average CoM distance must be > 0\n");

}

Points3f parseInput()
{
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

using namespace elfin;

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

EvolutionSolver * es;
bool esStarted = false;

void interruptHandler(int signal)
{
    raw("\n\n");
    wrn("Caught interrupt signal\n");

    // Save latest results
    if (esStarted)
    {
        wrn("Saving latest best solutions and exiting!\n");
        using namespace elfin;

        const Population & p = es->bestSoFar();

        for (int i = 0; i < p.size(); i++)
        {
            std::vector<std::string> nodeNames = p.at(i).getNodeNames();
            JSON nn = nodeNames;
            JSON j;
            j["nodes"] = nn;
            j["score"] = p.at(i).getScore();

            std::ostringstream ss;
            ss << options.outputDir << "/" << &p.at(i) << ".json";
            std::string dump = j.dump();
            const char * data = dump.c_str();
            const size_t len = dump.size();
            write_binary(ss.str().c_str(), data, len);
        }

        delete es;
    }
    else
    {
        wrn("GA did not get to start\n");
    }

    exit(1);
}

int runUnitTests()
{
    msg("Running unit tests...\n");
    int failCount = 0;
    failCount += _testMathUtils();
    failCount += _testKabsch();
    failCount += _testChromosome();
}

int runMetaTests(const Points3f & spec)
{
    msg("Running meta tests...\n");
    int failCount = 0;

    Points3f movedSpec = spec;

    Vector3f rotArr[3] = {
        Vector3f(1.0f, 0.0f, 0.0f),
        Vector3f(0.0f, -0.5177697998f, 0.855519979f),
        Vector3f(0.0f, -0.855519979f, -0.5177697998f)
    };
    Mat3x3 rotAroundX = Mat3x3(rotArr);

    // It seems that Kabsch cannot handle very large
    // translations, but in the scale we're working
    // at it should rarely, if ever, go beyond
    // one thousand Angstroms
    Vector3f tran(-39.0f, 999.3413f, -400.11f);

    for (Point3f &p : movedSpec)
    {
        p = p.dot(rotAroundX);
        p += tran;
    }

    // Test scoring a transformed version of spec
    Crc32 crc;
    const float trxScore = kabschScore(movedSpec, spec, crc);
    if (!float_approximates(trxScore, 0))
    {
        failCount ++;
        wrn("Self score test failed: self score should be 0\n");
    }

    // Test randomiser
    const int N = 10;
    const int randTrials = 50000000;
    const int expectAvg = randTrials / N;
    const float randDevTolerance = 0.05f * expectAvg; //5% deviation

    int randCount[N] = {0};
    for (int i = 0; i < randTrials; i++)
    {
        const int dice = getDice(N);
        if (dice >= N)
        {
            failCount;
            err("Failed to produce correct dice: getDice() produced %d for [0-%d)",
                dice, N);
            break;
        }
        randCount[dice]++;
    }

    for (int i = 0; i < N; i++)
    {
        const float randDev = (float) abs(randCount[i] - expectAvg) / (expectAvg);
        if (randDev > randDevTolerance)
        {
            failCount++;
            err("Too much random deviation: %.3f%% (expecting %d)\n",
                randDev, expectAvg);
        }
    }

    // Test parallel randomiser
#ifndef _NO_OMP
    std::vector<uint> paraRandSeeds = getParaRandSeeds();
    const int nThreads = paraRandSeeds.size();
    const int paraRandN = 8096;
    const long diceLim = 13377331;

    std::vector<uint> rands1(paraRandN);
    #pragma omp parallel for
    for (int i = 0; i < paraRandN; i++)
        rands1.at(i) = getDice(diceLim);

    getParaRandSeeds() = paraRandSeeds;
    std::vector<uint> rands2(paraRandN);
    #pragma omp parallel for
    for (int i = 0; i < paraRandN; i++)
        rands2.at(i) = getDice(diceLim);

    for (int i = 0; i < paraRandN; i++)
    {
        if (rands1.at(i) != rands2.at(i))
        {
            failCount++;
            err("Parallel randomiser failed: %d vs %d\n",
                rands1.at(i), rands2.at(i));
        }
    }
#endif

    return failCount;
}

int main(int argc, const char ** argv)
{
    std::signal(SIGINT, interruptHandler);

    // Default set to warning and above
    set_log_level(LOG_WARN);

    // Parse user arguments first to potentially get a settings file path
    parse_args(argc, argv, ARG_BUND_SIZE, argb, argParseFail);

    // Get values from settings json first
    parseSettings();

    // Parse user arguments a second time to override settings file
    parse_args(argc, argv, ARG_BUND_SIZE, argb, argParseFail);

    checkOptions();

    msg("Using master seed: %d\n", options.randSeed);

    RelaMat relaMat;
    NameIdMap nameIdMap;
    IdNameMap idNameMap;
    RadiiList radiiList;
    JSONParser().parseDB(options.xDBFile, nameIdMap, idNameMap, relaMat, radiiList);

    Gene::setup(&idNameMap);
    setupParaUtils(options.randSeed);

    Points3f spec = parseInput();

    if (options.runUnitTests)
    {
        int failCount = 0;
        failCount += runUnitTests();
        failCount += runMetaTests(spec);

        if (failCount > 0)
        {
            die("Some unit tests failed\n");
        }
        else
        {
            msg("Passed!\n");
        }
    }
    else
    {
        es = new EvolutionSolver(relaMat,
                                 spec,
                                 radiiList,
                                 options);
        esStarted = true;

        es->run();

        const Population & p = es->population();

        // Output best N solutions
        const uint outputN = 3;
        for (int i = 0; i < outputN; i++)
        {
            std::vector<std::string> nodeNames = p.at(i).getNodeNames();
            JSON nn = nodeNames;
            JSON j;
            j["nodes"] = nn;
            j["score"] = p.at(i).getScore();

            std::ostringstream ss;
            ss << options.outputDir << "/" << &p.at(i) << ".json";
            const char * data = j.dump().c_str();
            const size_t len = j.dump().size();
            write_binary(ss.str().c_str(), data, len);
        }

        delete es;
    }

    return 0;
}