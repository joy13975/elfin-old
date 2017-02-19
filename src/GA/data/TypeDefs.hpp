#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <vector>
#include <map>
#include <sstream>
#include <string>

#include "util.h"
#include "PrimitiveShorthands.hpp"

namespace elfin
{

class PairRelationship;

// Shorthands
typedef std::vector<PairRelationship *> RelaRow;
typedef std::vector<RelaRow> RelaMat;
typedef std::map<std::string, unsigned int> NameIdMap;
typedef std::vector<float>::const_iterator FloatConstIterator;

typedef std::vector<std::string> Solution;	// A solution is a series of node names

struct Vector3f
{
	float x, y, z;

	Vector3f() :
		x(0), y(0), z(0) {}

	Vector3f(const std::vector<float> & v) :
		Vector3f(v.begin(), v.end())
	{}

	Vector3f(FloatConstIterator begin,
	         FloatConstIterator end)
	{
		if ((end - begin) != 3)
			die("Vector3f() not called with vector range of length 3\n");

		FloatConstIterator itr = begin;
		x = *itr++;
		y = *itr++;
		z = *itr++;
	}

	std::string toString() const
	{
		std::ostringstream ss;
		ss << "v3f[" << x << ", " << y << ", " << z << ']';
		return ss.str();
	}
};

typedef Vector3f Point3f;

struct Mat3x3
{
	Vector3f rows[3];

	Mat3x3(const std::vector<float> & v) :
		Mat3x3(v.begin(), v.end())
	{}

	Mat3x3(FloatConstIterator begin,
	       FloatConstIterator end)
	{
		if ((end - begin) != 9)
			die("Mat3x3() not called with vector range of length 9\n");

		FloatConstIterator itr = begin;
		for (int i = 0; i < 3; i++)
		{
			const FloatConstIterator endItrTmp = itr + 3;
			rows[i] = Vector3f(itr, endItrTmp);
			itr += 3;
		}
	}

	std::string toString() const
	{
		std::ostringstream ss;

		ss << "m3x3[" << std::endl;
		for (int i = 0; i < 3; i++)
			ss << "      row" << (i + 1) << ":" << rows[i].toString() << std::endl;
		ss << " ]";

		return ss.str();
	}
};

typedef std::vector<Point3f> Points3f;

} // namespace elfin

#endif /* include guard */