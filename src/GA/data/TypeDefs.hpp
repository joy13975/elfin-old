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
typedef std::vector<PairRelationship *> RelaRow;
typedef std::vector<RelaRow> RelaMat;
typedef std::vector<std::string> Solution;	// A solution is a series of node names

template <typename T>
using Matrix = std::vector<std::vector<T>>;

struct Radii
{
	const float avgAll;
	const float maxCA;
	const float maxHeavy;
	Radii(float aa, float mca, float mh) :
		avgAll(aa), maxCA(mca), maxHeavy(mh)
	{}
};
typedef std::vector<Radii> RadiiList;

struct Gene
{
	uint nodeId;
	Point3f com;
	Gene(const uint _nodeId,
	     const Point3f _com) :
		nodeId(_nodeId),
		com(_com)
	{}
	Gene(const uint _nodeId,
	     const float x,
	     const float y,
	     const float z) :
		nodeId(_nodeId),
		com(x, y, z)
	{}
};
typedef std::vector<Gene> Genes;

} // namespace elfin

#endif /* include guard */