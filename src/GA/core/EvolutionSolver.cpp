#include <cmath>
#include <sstream>

#include "EvolutionSolver.hpp"
#include "util.h"

namespace elfin
{

// Constructors

EvolutionSolver::EvolutionSolver(const RelaMat & _relaMat, const Points3f & _spec) :
	myRelaMat(_relaMat), mySpec(_spec)
{
}

// Public methods

void EvolutionSolver::run(const ulong popSize, const ulong nIters)
{
	this->printStartMsg(popSize, nIters);

	for (auto & p : mySpec)
		dbg("Spec Point: %s\n", p.toString().c_str());

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

void EvolutionSolver::evolvePopulation()
{
	// Make new generation by discarding unfit
	// genes and mutating/reproducing fit genes
	// to fill up discarded slots
	wrn("TODO: population evolution\n");
}

void EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)
	wrn("TODO: population ranking\n");
}

void EvolutionSolver::scorePopulation()
{
	for (auto & gene : myPopulation)
		gene.score(mySpec);
}

void EvolutionSolver::initPopulation(const ulong popSize)
{
	const uint expLen = Gene::calcExpectedLength(mySpec);

	myPopulation = Population();
	for (int i = 0; i < popSize; i++)
		myPopulation.emplace_back(
		    Gene::randomChromosome(
		        expLen,
		        myGeneLenDev,
		        myRelaMat)
		);
}

void EvolutionSolver::printStartMsg(const ulong popSize, const ulong nIters)
{
	// Want auto significant figure detection with streams
	std::ostringstream psStr;
	psStr << (float) (popSize / 1000.0f) << "k";
	std::ostringstream niStr;
	niStr << (float) (nIters / 1000.0f) << "k";

	msg("EvolutionSolver running with popSize %s, %s iterations\n",
	    psStr.str().c_str(), niStr.str().c_str());
}

void EvolutionSolver::printEndMsg()
{
	msg("EvolutionSolver finished: ");
	this->printTiming();
}

void EvolutionSolver::startTimer()
{
	myStartTimeInUs = get_timestamp_us();
}

void EvolutionSolver::printTiming()
{
	const double timeElapsedInUs = get_timestamp_us() - myStartTimeInUs;
	const uint minutes = std::floor(timeElapsedInUs / 1e6 / 60.0f);
	const uint seconds = std::floor(fmod(timeElapsedInUs / 1e6, 60.0f));
	const uint milliseconds = std::floor(fmod(timeElapsedInUs / 1e3, 1000.0f));
	raw("%um %us %ums\n",
	    minutes, seconds, milliseconds);
}

} // namespace elfin
