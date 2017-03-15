#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_

#include <cmath>
#include <string>

#include "TypeDefs.hpp"
#include "Gene.hpp"
#include "../core/Checksum.hpp"

namespace elfin
{

class Chromosome
{
public:
	Chromosome();
	Chromosome(const Chromosome & rhs);
	Chromosome(const Genes & genes);
	// Chromosome & operator=(const Chromosome & rhs);
	virtual ~Chromosome() {};

	void score(const Points3f & ref);
	bool operator>(const Chromosome & rhs) const;
	bool operator<(const Chromosome & rhs) const;

	// Getter & setters
	float getScore() const;
	Genes & genes();
	const Genes & genes() const;
	Crc32 checksum() const;
	std::vector<std::string> getNodeNames() const;

	std::string toString() const;
	bool cross(const Chromosome & father,
	           const Chromosome & mother,
	           const IdPairs & crossingIds);
	void inheritMutate(const Chromosome & parent);
	void autoMutate();
	void randomise();
	bool pointMutate();
	bool limbMutate();

	IdPairs findCompatibleCrossings(const Chromosome & other) const;

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
	float myScore = NAN;
	Crc32 myChecksum = 0;

	static bool setupDone;
	static uint myMinLen;
	static uint myMaxLen;
	static const RelaMat * myRelaMat;
	static const RadiiList * myRadiiList;
	static IdPairs myNeighbourCounts;
	static IdRoulette myGlobalRoulette;

	Genes genRandomGenesReverse(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());
	Genes genRandomGenes(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());
};

int _testChromosome();
} // namespace elfin

#endif /* include guard */