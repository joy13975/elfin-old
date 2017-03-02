#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

#include <string>
#include <vector>

namespace elfin
{

typedef std::vector<float>::const_iterator FloatConstIterator;
class Mat3x3;

class Vector3f
{
public:
	float x, y, z;

	Vector3f() :
		x(0), y(0), z(0) {}

	Vector3f(float _x, float _y, float _z) :
		x(_x), y(_y), z(_z) {}

	Vector3f(const std::vector<float> & v) :
		Vector3f(v.begin(), v.end())
	{}

	Vector3f(FloatConstIterator begin,
	         FloatConstIterator end);

	std::string toString() const;

	Vector3f operator+(const Vector3f & rhs) const;
	Vector3f operator-(const Vector3f & rhs) const;
	Vector3f operator*(const float f) const;
	Vector3f & operator+=(const Vector3f & rhs);
	Vector3f & operator-=(const Vector3f & rhs);
	float dot(const Vector3f & rhs) const;
	Vector3f dot(const Mat3x3 & rotMat) const;
	float distTo(const Vector3f & rhs) const;

	bool approximates(const Vector3f & ref, double tolerance = 1e-5);
};
typedef Vector3f Point3f;
typedef std::vector<Point3f> Points3f;

class Mat3x3
{
public:
	Vector3f rows[3];

	Mat3x3(Vector3f _rows[3]);

	Mat3x3(const std::vector<float> & v) :
		Mat3x3(v.begin(), v.end())
	{}

	Mat3x3(FloatConstIterator begin,
	       FloatConstIterator end);

	Vector3f dot(const Vector3f & rotMat) const;
	Mat3x3 dot(const Mat3x3 & rotMat) const;
	Mat3x3 transpose() const;

	std::string toString() const;
};


} // namespace elfin

#endif /* include guard */