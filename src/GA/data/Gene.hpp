#ifndef _GENE_HPP_
#define _GENE_HPP_

#include <cmath>

#include "TypeDefs.hpp"

namespace elfin
{

typedef std::vector<uint> Chromosome; 		// A chromosome is a series of node indexes

class Gene
{
public:

	Gene();
	Gene(const Chromosome & chromoRef);
	virtual ~Gene() {};

	void score(const Points3f & ref);

	static uint calcExpectedLength(const Points3f & lenRef);

	static Chromosome randomChromosome(
	    const uint expLen,
	    const float lenDev,
	    const RelaMat & relaMat);
private:
	Chromosome myChromo;
	float myScore = NAN;

	float calcRMSD(const Points3f & ref, Points3f & toMove);
	void superimpose(const Points3f & ref, Points3f & toMove);
	Points3f growShape();
};

} // namespace elfin

#endif /* include guard */