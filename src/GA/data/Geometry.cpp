#include <sstream>
#include <cmath>

#include "Geometry.hpp"
#include "util.h"

namespace elfin
{


Vector3f::Vector3f(FloatConstIterator begin,
                   FloatConstIterator end)
{
	if ((end - begin) != 3)
		die("Vector3f() not called with vector range of length 3\n");

	FloatConstIterator itr = begin;
	x = *itr++;
	y = *itr++;
	z = *itr++;
}

std::string
Vector3f::toString() const
{
	std::ostringstream ss;
	ss << "v3f[" << x << ", " << y << ", " << z << ']';
	return ss.str();
}

Vector3f
Vector3f::operator+(const Vector3f & rhs) const
{
	return Vector3f(
	           rhs.x + this->x,
	           rhs.y + this->y,
	           rhs.z + this->z);
}

Vector3f
Vector3f::operator-(const Vector3f & rhs) const
{
	return Vector3f(
	           this->x - rhs.x,
	           this->y - rhs.y,
	           this->z - rhs.z);
}

Vector3f
Vector3f::operator*(const float f) const
{
	return Vector3f(
	           f * this->x,
	           f * this->y,
	           f * this->z);
}

Vector3f &
Vector3f::operator+=(const Vector3f& rhs)
{
	this->x += rhs.x;
	this->y += rhs.y;
	this->z += rhs.z;

	return *this;
}

float
Vector3f::dot(const Vector3f & rhs) const
{
	return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z;
}

Vector3f &
Vector3f::rotate(const Mat3x3 & rotMat)
{
	const float tx = this->x * rotMat.rows[0].x +
	                 this->y * rotMat.rows[0].y +
	                 this->z * rotMat.rows[0].z;

	const float ty = this->x * rotMat.rows[1].x +
	                 this->y * rotMat.rows[1].y +
	                 this->z * rotMat.rows[1].z;

	const float tz = this->x * rotMat.rows[2].x +
	                 this->y * rotMat.rows[2].y +
	                 this->z * rotMat.rows[2].z;

	this->x = tx;
	this->y = ty;
	this->z = tz;

	return *this;
}

float
Vector3f::distTo(const Point3f & rhs) const
{
	const float dx = (this->x - rhs.x);
	const float dy = (this->y - rhs.y);
	const float dz = (this->z - rhs.z);
	return sqrt(dx * dx + dy * dy + dz * dz);
}

bool
Vector3f::approximates(const Vector3f & ref, double tolerance)
{
	if (this->x != ref.x ||
	        this->y != ref.y ||
	        this->z != ref.z)
	{
		const float dx = this->x - ref.x;
		const float dy = this->y - ref.y;
		const float dz = this->z - ref.z;
		wrn("Point3f diffs: %.8f, %.8f, %.8f\n",
		    dx, dy, dz);



		if (!float_approximates_err(
		            dx, 0.0, tolerance) ||
		        !float_approximates_err(
		            dy, 0.0, tolerance) ||
		        !float_approximates_err(
		            dz, 0.0, tolerance))
		{
			return false;
		}
	}

	return true;
}

Mat3x3::Mat3x3(Vector3f _rows[3])
{
	rows[0] = _rows[0];
	rows[1] = _rows[1];
	rows[2] = _rows[2];
}
Mat3x3::Mat3x3(FloatConstIterator begin,
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

std::string
Mat3x3::toString() const
{
	std::ostringstream ss;

	ss << "m3x3[" << std::endl;
	for (int i = 0; i < 3; i++)
		ss << "      row" << (i + 1) << ":" << rows[i].toString() << std::endl;
	ss << " ]";

	return ss.str();
}

} // namespace elfin