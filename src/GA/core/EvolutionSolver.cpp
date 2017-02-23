#include <cmath>
#include <sstream>

#include "EvolutionSolver.hpp"
#include "util.h"

namespace elfin
{

// Constructors

EvolutionSolver::EvolutionSolver(const RelaMat & relaMat,
                                 const Points3f & spec,
                                 const RadiiList & radiiList,
                                 const float chromoLenDev,
                                 const float avgPairDist) :
	myRelaMat(relaMat),
	mySpec(spec),
	myRadiiList(radiiList),
	myChromoLenDev(chromoLenDev),
	myExpectedTargetLen(Chromosome::calcExpectedLength(spec, avgPairDist))
{
}

// Public methods

void
EvolutionSolver::run(const ulong popSize,
                     const ulong nIters)
{
	this->printStartMsg(popSize, nIters);

	this->startTimer();

	initPopulation(popSize);

	for (int i = 0; i < nIters; i++)
	{
		scorePopulation();

		rankPopulation();

		evolvePopulation();
	}

	this->printEndMsg();
}

// Private methods

void
EvolutionSolver::evolvePopulation()
{
	// Make new generation by discarding unfit
	// genes and mutating/reproducing fit genes
	// to fill up discarded slots
	wrn("TODO: population evolution\n");
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)
	wrn("TODO: population ranking\n");
}

void
EvolutionSolver::scorePopulation()
{
	for (auto & chromo : myPopulation)
		chromo.score(mySpec);
}

void
EvolutionSolver::initPopulation(const ulong popSize)
{
	const uint minLen = myExpectedTargetLen * (1 - myChromoLenDev);
	const uint maxLen = myExpectedTargetLen * (1 + myChromoLenDev);

	myPopulation = Population();
	myPopulation.reserve(popSize);
	for (int i = 0; i < popSize; i++)
		myPopulation.emplace_back(
		    Chromosome::genRandomGenes(
		        minLen,
		        maxLen,
		        myRelaMat,
		        myRadiiList)
		);
}

void
EvolutionSolver::printStartMsg(const ulong popSize,
                               const ulong nIters)
{
	for (auto & p : mySpec)
		dbg("Spec Point: %s\n", p.toString().c_str());

	msg("Expecting length: %u\n", myExpectedTargetLen);
	msg("Using deviation allowance: %.1f%%\n", myChromoLenDev * 100);

	// Want auto significant figure detection with streams
	std::ostringstream psStr;
	psStr << (float) (popSize / 1000.0f) << "k";
	std::ostringstream niStr;
	niStr << (float) (nIters / 1000.0f) << "k";

	msg("EvolutionSolver starting with popSize %s, %s iterations\n",
	    psStr.str().c_str(), niStr.str().c_str());
}

void
EvolutionSolver::printEndMsg()
{
	msg("EvolutionSolver finished: ");
	this->printTiming();
}

void
EvolutionSolver::startTimer()
{
	myStartTimeInUs = get_timestamp_us();
}

void
EvolutionSolver::printTiming()
{
	const double timeElapsedInUs = get_timestamp_us() - myStartTimeInUs;
	const uint minutes = std::floor(timeElapsedInUs / 1e6 / 60.0f);
	const uint seconds = std::floor(fmod(timeElapsedInUs / 1e6, 60.0f));
	const uint milliseconds = std::floor(fmod(timeElapsedInUs / 1e3, 1000.0f));
	raw("%um %us %ums\n",
	    minutes, seconds, milliseconds);
}

} // namespace elfin
