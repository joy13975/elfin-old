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
	void autoMutate();
	void randomise();
	bool pointMutate();
	bool limbMutate();
	void setGenes(const Genes & genes);
	Genes & getGenes();

	static void setup(const uint minLen,
	                  const uint maxLen,
	                  const RelaMat & relaMat,
	                  const RadiiList & radiiList);
	static uint calcExpectedLength(const Points3f & lenRef,
	                               const float avgPairDist);
	static bool synthesiseReverse(Genes & genes);
	static bool synthesise(Genes & genes);
private:
	Genes myGenes;
	static bool setupDone;
	static uint myMinLen;
	static uint myMaxLen;
	static const RelaMat * myRelaMat;
	static const RadiiList * myRadiiList;
	static std::vector<std::tuple<uint, uint>> myNeighbourCounts;
	static IdRoulette myGlobalRoulette;

	float myScore = NAN;

	Genes genRandomGenesReverse(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());
	Genes genRandomGenes(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());
};

} // namespace elfin

#endif /* include guard */