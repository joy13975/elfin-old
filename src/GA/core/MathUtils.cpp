#ifdef _TEST_MATHUTILS

#include <cmath>

#include "MathUtils.hpp"
#include "util.h"
#include "../data/TypeDefs.hpp"

using namespace elfin;

int failCount = 0;

int main(int argc, const char ** argv)
{
	msg("Testing MathUtils\n");

	// Translate 1
	Point3f a(1.0f, 2.0f, 3.0f);
	Vector3f t(9.0f, 9.0f, 9.0f);
	msg("Translate point %s using %s = ",
	    a.toString().c_str(),
	    t.toString().c_str());

	translate(a, t);
	raw("%s\n", a.toString().c_str());
	if (!a.approximates(Point3f(10.0f, 11.0f, 12.0f)))
		failCount++;

	// Translate 2
	t = Vector3f(-3.0f, 100.0f, 493.1337f);
	msg("Translate point %s using %s = ",
	    a.toString().c_str(),
	    t.toString().c_str());

	translate(a, t);
	raw("%s\n", a.toString().c_str());
	if (!a.approximates(Point3f(7.0f, 111.0f, 505.1337f)))
		failCount++;

	// Rotate 1
	Vector3f rotRows1[3] = {
		Vector3f(1.0f, 0.0f, 0.0f),
		Vector3f(0.0f, 1.0f, 0.0f),
		Vector3f(0.0f, 0.0f, 1.0f)
	};
	Mat3x3 r = Mat3x3(rotRows1);
	msg("Rotate point %s using %s = ",
	    a.toString().c_str(),
	    r.toString().c_str());

	rotate(a, r);
	raw("%s\n", a.toString().c_str());
	if (!a.approximates(Point3f(7.0f, 111.0f, 505.1337f)))
		failCount++;

	// Rotate 2
	Vector3f rotRows2[3] = {
		Vector3f(0.4f, 0.5f, 0.0f),
		Vector3f(0.5f, 1.0f, 0.0f),
		Vector3f(0.0f, 0.0f, 1.0f)
	};
	r = Mat3x3(rotRows2);
	msg("Rotate point %s using %s = ",
	    a.toString().c_str(),
	    r.toString().c_str());

	rotate(a, r);
	raw("%s\n", a.toString().c_str());
	if (!a.approximates(Point3f(58.3f, 114.5f, 505.1337f)))
		failCount++;

	// Rotate + translate
	t = Vector3f(-9.32f, 1.001f, -0.1337f);
	Vector3f rotRows3[3] = {
		Vector3f(0.4f, 0.1f, 0.3f),
		Vector3f(0.5f, 0.1f, 0.53f),
		Vector3f(0.9f, 0.0f, 0.01f)
	};
	r = Mat3x3(rotRows3);
	msg("Rotate point %s using %s\nand then translate with %s = \n",
	    a.toString().c_str(),
	    r.toString().c_str(),
	    t.toString().c_str());

	rotate(a, r);
	translate(a, t);
	raw("%s\n", a.toString().c_str());
	if (!a.approximates(Point3f(176.99011f, 309.321861f, 57.387637f)))
		failCount++;

	if (failCount > 0)
	{
		err("Some tests failed - total %d\n", failCount);
	}
	else
	{
		msg("All tests passed!\n");
	}
	return 0;
}

#endif //ifdef _TEST_MATHUTILS
