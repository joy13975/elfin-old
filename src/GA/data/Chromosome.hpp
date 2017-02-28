#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_

#include <cmath>
#include <string>

#include "TypeDefs.hpp"
#include "Gene.hpp"

namespace elfin
{

class Chromosome
{
public:
	Chromosome();
	// Chromosome(const Chromosome & rhs);
	// Chromosome & operator=(const Chromosome & rhs);
	virtual ~Chromosome() {};

	void score(const Points3f & ref);
	bool operator>(const Chromosome & rhs) const;
	bool operator<(const Chromosome & rhs) const;
	float getScore() const;
	std::string toString() const;
	void randomise();
	void mutate();

	static void setup(const uint minLen,
	                  const uint maxLen,
	                  const RelaMat & relaMat,
	                  const RadiiList & radiiList);
	static uint calcExpectedLength(const Points3f & lenRef,
	                               const float avgPairDist);
private:
	Genes myGenes;
	static bool setupDone;
	static uint myMinLen;
	static uint myMaxLen;
	static const RelaMat * myRelaMat;
	static const RadiiList * myRadiiList;

	float myScore = NAN;

	bool pointMutate();
	void limbMutate();
	void syntehsise();
};

} // namespace elfin

#endif /* include guard */