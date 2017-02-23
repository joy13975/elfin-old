#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_

#include <cmath>

#include "TypeDefs.hpp"

namespace elfin
{

class Chromosome
{
public:

	Chromosome();
	Chromosome(const Genes & geneRef);
	virtual ~Chromosome() {};

	void score(const Points3f & ref);

	static uint calcExpectedLength(const Points3f & lenRef,
	                               const float avgPairDist);

	static Genes genRandomGenes(
	    const uint minLen,
	    const uint maxLen,
	    const RelaMat & relaMat,
	    const RadiiList & radiiList);
private:
	Genes myGenes;
	float myScore = NAN;

	void superimpose(const Points3f & ref);
};

} // namespace elfin

#endif /* include guard */