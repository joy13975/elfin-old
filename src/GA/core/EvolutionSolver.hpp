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
	                const float chromoLenDev,
	                const float avgPairDist);
	virtual ~EvolutionSolver() {};

	void run(const ulong popSize, const ulong nIters);
private:
	const RelaMat & myRelaMat;
	const Points3f & mySpec;
	const RadiiList & myRadiiList;
	const float myChromoLenDev;
	const uint myExpectedTargetLen;

	double myStartTimeInUs = 0;
	Population myPopulation;

	void evolvePopulation();
	void rankPopulation();
	void scorePopulation();
	void initPopulation(const ulong popSize);

	void printStartMsg(const ulong popSize, const ulong nIters);
	void printEndMsg();
	void startTimer();
	void printTiming();
};

} // namespace elfin

#endif /* include guard */