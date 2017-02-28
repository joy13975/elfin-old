#include <cmath>
#include <sstream>
#include <algorithm>

#include "EvolutionSolver.hpp"
#include "util.h"


namespace elfin
{

// Constructors

EvolutionSolver::EvolutionSolver(const RelaMat & relaMat,
                                 const Points3f & spec,
                                 const RadiiList & radiiList,
                                 const OptionPack & options) :
	myRelaMat(relaMat),
	mySpec(spec),
	myRadiiList(radiiList),
	myOptions(options)
{
	double timeSeed = get_timestamp_us();
	std::srand(options.randSeed == 0 ? timeSeed : options.randSeed);

	mySruviverCutoff = std::round(options.gaSurviveRate * options.gaPopSize);

	const ulong nNonSurvivers = (options.gaPopSize - mySruviverCutoff);
	myCrossCutoff = mySruviverCutoff + std::round(options.gaCrossRate * nNonSurvivers);
	myMutateCutoff = std::min(
	                     (ulong) std::round(myCrossCutoff + options.gaMutateRate * nNonSurvivers),
	                     nNonSurvivers);

	myExpectedTargetLen = Chromosome::calcExpectedLength(spec, options.avgPairDist);
	const float minTargetLen = myExpectedTargetLen * (1 - myOptions.chromoLenDev);
	const float maxTargetLen = myExpectedTargetLen * (1 + myOptions.chromoLenDev);

	Chromosome::setup(minTargetLen, maxTargetLen, myRelaMat, myRadiiList);
}

// Public methods

void
EvolutionSolver::run()
{
	this->printStartMsg();

	this->startTimer();

	initPopulation();

	const int genDispDigits = std::ceil(std::log(myOptions.gaIters) / std::log(10));
	char * genMsgFmt;
	asprintf(&genMsgFmt,
	         "Generation #%%%dd: fittest=%%.2f\n", genDispDigits);
	for (int i = 0; i < myOptions.gaIters; i++)
	{
		scorePopulation();

		rankPopulation();

		evolvePopulation();

		pruneColliders();

		msg(genMsgFmt, i, myPopulation.front().getScore());
	}

	this->printEndMsg();
}

// Private methods

void
EvolutionSolver::pruneColliders()
{
	// Instead of enforcing collision-free mutation,
	// it is simpler to mutate without checking and
	// just assign a bad fitness to shapes that end
	// up in self collision
	wrn("TODO: invalidate chromosomes that have collision\n");
}

void
EvolutionSolver::evolvePopulation()
{
	// Make new generation by discarding unfit
	// genes and mutating/reproducing fit genes
	// to fill up discarded slots

	wrn("TODO: cross population\n");

	for (int i = myCrossCutoff; i < myMutateCutoff; i++)
		myPopulation.at(i).mutate();

	for (int i = myMutateCutoff; i < myOptions.gaPopSize; i++)
		myPopulation.at(i).randomise();
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)

	std::sort(myPopulation.begin(),
	          myPopulation.end());

	// Print score
	for (int i = 0; i < myPopulation.size(); i++)
		prf("myPopulation[%d] Score: %.2f\n", i, myPopulation.at(i).getScore());
}

void
EvolutionSolver::scorePopulation()
{
	for (auto & chromo : myPopulation)
		chromo.score(mySpec);
}

void
EvolutionSolver::initPopulation()
{
	myPopulation = Population(myOptions.gaPopSize);

	for (int i = 0; i < myOptions.gaPopSize; i++)
		myPopulation.at(i).randomise();
}

void
EvolutionSolver::printStartMsg()
{
	for (auto & p : mySpec)
		dbg("Spec Point: %s\n", p.toString().c_str());

	msg("Expecting length: %u\n", myExpectedTargetLen);
	msg("Using deviation allowance: %.1f%%\n", myOptions.chromoLenDev * 100);

	// Want auto significant figure detection with streams
	std::ostringstream psStr;
	if (myOptions.gaPopSize > 1000)
		psStr << (float) (myOptions.gaPopSize / 1000.0f) << "k";
	else
		psStr << myOptions.gaPopSize;

	std::ostringstream niStr;
	if (myOptions.gaIters > 1000)
		niStr << (float) (myOptions.gaIters / 1000.0f) << "k";
	else
		niStr << myOptions.gaIters;


	msg("EvolutionSolver starting with following settings:\n"
	    "Population size:      %s\n"
	    "Iterations:           %s\n"
	    "Survive cutoff:       %u\n"
	    "Cross cutoff:         %u\n"
	    "Mutate cutoff:        %u\n"
	    "New species:          %u\n",
	    psStr.str().c_str(),
	    niStr.str().c_str(),
	    mySruviverCutoff,
	    myCrossCutoff,
	    myMutateCutoff,
	    myOptions.gaPopSize - myMutateCutoff);
}

void
EvolutionSolver::printEndMsg()
{
	msg("EvolutionSolver finished: ");
	this->printTiming();

	// Print best N solutions
	const uint N = 3;

	for (int i = 0; i < N; i++)
	{
		const auto & p = myPopulation.at(i);
		msg("Solution #%d score %.2f: \n%s\n",
		    p.getScore(),
		    i,
		    p.toString().c_str());
	}
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
