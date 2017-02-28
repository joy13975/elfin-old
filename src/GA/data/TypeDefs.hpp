#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <vector>
#include <map>
#include <sstream>
#include <string>

#include "util.h"
#include "PrimitiveShorthands.hpp"
#include "Geometry.hpp"

namespace elfin
{

class PairRelationship;

// Shorthands
typedef std::map<std::string, uint> NameIdMap;
typedef std::map<uint, std::string> IdNameMap;
typedef std::vector<PairRelationship *> RelaRow;
typedef std::vector<RelaRow> RelaMat;
typedef std::vector<std::string> Solution;	// A solution is a series of node names

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
	std::string outputFile = "./output.json";

	float chromoLenDev = 0.2;
	// Average CoM distance found by xDBStat.py as of 28/Feb/2017
	float avgPairDist = 38.0f;

	// GA parameters
	uint randSeed = 0x1337cafe;
	long gaPopSize = 10000;
	long gaIters = 1000;
	float gaSurviveRate = 0.1;
	float gaCrossRate = 0.6;
	float gaMutateRate = 0.3;

};

} // namespace elfin

#endif /* include guard */