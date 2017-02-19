#ifndef _EVOLUTIONSOLVER_HPP_
#define _EVOLUTIONSOLVER_HPP_

#include "../data/TypeDefs.hpp"
#include "../data/Gene.hpp"

namespace elfin
{

typedef std::vector<Gene> Population;		// A population is a pool of genes

class EvolutionSolver
{
public:
	EvolutionSolver(const RelaMat & rm, const Points3f & spec);
	virtual ~EvolutionSolver() {};

	void run(const ulong popSize, const ulong nIters);

private:
	const RelaMat myRelaMat;
	const Points3f mySpec;
	const float myGeneLenDev = 0.2f;

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