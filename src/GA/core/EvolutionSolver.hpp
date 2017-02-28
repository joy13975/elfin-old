#ifndef _EVOLUTIONSOLVER_HPP_
#define _EVOLUTIONSOLVER_HPP_

#include "../data/TypeDefs.hpp"
#include "../data/Chromosome.hpp"

namespace elfin
{

typedef std::vector<Chromosome> Population;

class EvolutionSolver
{
public:
	EvolutionSolver(const RelaMat & relaMat,
	                const Points3f & spec,
	                const RadiiList & radiiList,
	                const OptionPack & options);
	virtual ~EvolutionSolver() {};

	void run();
private:
	const RelaMat & myRelaMat;
	const Points3f & mySpec;
	const RadiiList & myRadiiList;
	const OptionPack & myOptions;

	uint myExpectedTargetLen;
	ulong mySruviverCutoff;
	ulong myCrossCutoff;
	ulong myMutateCutoff;

	double myStartTimeInUs = 0;
	Population myPopulation;

	void pruneColliders();
	void evolvePopulation();
	void rankPopulation();
	void scorePopulation();
	void initPopulation();

	void printStartMsg();
	void printEndMsg();
	void startTimer();
	void printTiming();
};

} // namespace elfin

#endif /* include guard */