#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <tuple>

#include "util.h"
#include "PrimitiveShorthands.hpp"
#include "Geometry.hpp"

#define toCString toString().c_str

namespace elfin
{

class PairRelationship;

// Shorthands
typedef std::map<std::string, int> NameIdMap;
typedef std::map<int, std::string> IdNameMap;

typedef std::vector<PairRelationship *> RelaRow;
typedef std::vector<RelaRow> RelaMat;

typedef std::vector<int> IdRoulette;
typedef std::vector<int> Ids;
typedef std::tuple<int, int> IdPair;
typedef std::vector<IdPair> IdPairs;

template <typename T>
using Matrix = std::vector<std::vector<T>>;

struct Radii
{
	float avgAll;
	float maxCA;
	float maxHeavy;
	Radii(float aa, float mca, float mh) :
		avgAll(aa), maxCA(mca), maxHeavy(mh)
	{}
};
typedef std::vector<Radii> RadiiList;

struct OptionPack
{
	// Input settings
	std::string xDBFile = "./xDB.json";
	std::string inputFile = "";

	enum InputType { Unknown, CSV, JSON };
	InputType inputType = Unknown;
	std::string settingsFile = "./settings.json";
	std::string outputDir = "./out/";

	float chromoLenDev = 0.2;
	// Average CoM distance found by xDBStat.py as of 28/Feb/2017
	float avgPairDist = 38.0f;

	// GA parameters
	uint randSeed = 0x1337cafe;
	long gaPopSize = 10000;
	long gaIters = 1000;
	float gaSurviveRate = 0.1;
	float gaCrossRate = 0.5;
	float gaPointMutateRate = 0.5;
	float gaLimbMutateRate = 0.5;

};

} // namespace elfin

#endif /* include guard */