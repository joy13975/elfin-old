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

	const Population * population() const;
	const Population & bestSoFar() const;

	void run();
private:
	const RelaMat & myRelaMat;
	const Points3f & mySpec;
	const RadiiList & myRadiiList;
	const OptionPack & myOptions;

	uint myExpectedTargetLen;
	ulong myNonSurviverCount;
	ulong mySurviverCutoff;
	ulong myCrossCutoff;
	ulong myPointMutateCutoff;
	ulong myLimbMutateCutoff;

	double myStartTimeInUs = 0;
	Population myPopulationBuffers[2]; // double buffer
	const Population *myCurrPop;
	Population *myBuffPop;
	Population myBestSoFar; // Currently used for emergency output

	void initPopulation();
	void evolvePopulation();
	void scorePopulation();
	void rankPopulation();
	void selectParents();
	void swapPopBuffers();

	void printStartMsg();
	void printEndMsg();
	void startTimer();
	void printTiming();
};

} // namespace elfin

#endif /* include guard */