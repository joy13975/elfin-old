#ifndef _MATHUTILS_HPP_
#define _MATHUTILS_HPP_

#include <cmath>
#include <cstdlib>

#include "../data/TypeDefs.hpp"

// COLLISION_MEASURE is one of {avgAll, maxHeavy, maxCA}
#define COLLISION_MEASURE maxHeavy

namespace elfin
{

inline float
dot(const Vector3f & v1, const Vector3f & v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline void
translate(Point3f & point, const Vector3f & transMat)
{
	point.x += transMat.x;
	point.y += transMat.y;
	point.z += transMat.z;
}

inline void
rotate(Point3f & point, const Mat3x3 & rotMat)
{
	const float tx = point.x * rotMat.rows[0].x +
	                 point.y * rotMat.rows[0].y +
	                 point.z * rotMat.rows[0].z;

	const float ty = point.x * rotMat.rows[1].x +
	                 point.y * rotMat.rows[1].y +
	                 point.z * rotMat.rows[1].z;

	const float tz = point.x * rotMat.rows[2].x +
	                 point.y * rotMat.rows[2].y +
	                 point.z * rotMat.rows[2].z;

	point.x = tx;
	point.y = ty;
	point.z = tz;
}

inline float
norm(const Point3f & p1, const Point3f & p2)
{
	const float dx = (p1.x - p2.x);
	const float dy = (p1.y - p2.y);
	const float dz = (p1.z - p2.z);
	return sqrt(dx * dx + dy * dy + dz * dz);
}

inline bool
collides(const uint newId,
         const Point3f & newCOM,
         const Genes & genes,
         const RadiiList & radiiList)
{
	// Check collision with all nodes up to previous PAIR
	const int lim = genes.size() - 2;
	for (int i = 0; i < lim; i++)
	{
		const float comDist = norm(genes.at(i).com, newCOM);
		const float requiredComDist = radiiList.at(i).COLLISION_MEASURE +
		                              radiiList.at(newId).COLLISION_MEASURE;
		if (comDist < requiredComDist)
			return true;
	}

	return false;
}

inline float
minDistFromLine(const Point3f & point,
                const Points3f & line)
{
	float minDist = INFINITY;

	for (int i = 1; i < line.size(); i++)
	{
		const Vector3f v = Vector3f(line[i].x - line[i - 1].x,
		                            line[i].y - line[i - 1].y,
		                            line[i].z - line[i - 1].z);
		const Vector3f w = Vector3f(point.x - line[i - 1].x,
		                            point.y - line[i - 1].y,
		                            point.z - line[i - 1].z);

		const float c1 = dot(w, v);
		float dist = NAN;
		if (c1 <= 0)
		{
			dist = norm(w, Vector3f(0, 0, 0));
		}
		else
		{
			const float c2 = dot(v, v);
			if (c2 <= c1)
			{
				dist = norm(Vector3f(point.x - line[i].x,
				                     point.y - line[i].y,
				                     point.z - line[i].z),
				            Vector3f(0, 0, 0));
			}
			else
			{
				const float b = c1 / c2;
				const Vector3f pol = Vector3f(line[i - 1].x - b * v.x,
				                              line[i - 1].y - b * v.y,
				                              line[i - 1].z - b * v.z);
				dist = norm(Vector3f(point.x - pol.x,
				                     point.y - pol.y,
				                     point.z - pol.z),
				            Vector3f(0, 0, 0));
			}
		}

		if (dist < minDist)
			minDist = dist;
	}

	return minDist;
}

} // namespace elfin

#endif /* include guard */