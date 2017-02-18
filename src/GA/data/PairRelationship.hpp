#ifndef _PAIRRELATIONSHIP_HPP_
#define _PAIRRELATIONSHIP_HPP_

#include <vector>
#include <sstream>

#include "util.h"
#include "TypeDefs.hpp"

namespace elfin
{

// Lazy local typedefs
typedef std::vector<float>::const_iterator FloatConstIterator;

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

class PairRelationship
{
public:
	PairRelationship(
	    const std::vector<float> & comBv,
	    const std::vector<float> & rotv,
	    const std::vector<float> & tranv) :
		comB(Point3f(comBv)),
		rot(Mat3x3(rotv)),
		tran(Vector3f(tranv))
	{};

	virtual ~PairRelationship() {};

	const Point3f comB;
	const Mat3x3 rot;
	const Vector3f tran;

	std::string toString()
	{
		std::ostringstream ss;

		ss << "pr[" << std::endl;
		ss << "    comB:" << comB.toString() << std::endl;
		ss << "    rot:" << rot.toString() << std::endl;
		ss << "    tran:" << tran.toString() << std::endl;
		ss << "]" << std::endl;

		return ss.str();
	}
};

} // namespace elfin

#endif /* include guard */