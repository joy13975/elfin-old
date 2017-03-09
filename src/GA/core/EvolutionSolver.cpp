#include <cmath>
#include <sstream>
#include <algorithm>
#include <omp.h>

#ifdef _NO_OMP
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

#include "EvolutionSolver.hpp"
#include "util.h"
#include "Parallel.hpp"

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

	myNonSurviverCount = (options.gaPopSize - mySruviverCutoff);
	myCrossCutoff = mySruviverCutoff + std::round(options.gaCrossRate * myNonSurviverCount);
	myPointMutateCutoff = myCrossCutoff + std::round(options.gaPointMutateRate * myNonSurviverCount);
	myLimbMutateCutoff = std::min(
	                         (ulong) (myPointMutateCutoff + std::round(options.gaLimbMutateRate * myNonSurviverCount)),
	                         (ulong) options.gaPopSize);

	myExpectedTargetLen = Chromosome::calcExpectedLength(spec, options.avgPairDist);
	const float minTargetLen = myExpectedTargetLen * (1 - myOptions.chromoLenDev);
	const float maxTargetLen = myExpectedTargetLen * (1 + myOptions.chromoLenDev);

	Chromosome::setup(minTargetLen, maxTargetLen, myRelaMat, myRadiiList);
}

const Population &
EvolutionSolver::getPopulation() const
{
	return myPopulation;
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
	         "Generation #%%%dd: best=%%.2f   worst=%%.2f   time taken=%%dms\n", genDispDigits);
	for (int i = 0; i < myOptions.gaIters; i++)
	{
		const double genStartTime = get_timestamp_us();
		scorePopulation();

		rankPopulation();

		evolvePopulation();

		msg(genMsgFmt, i,
		    myPopulation.front().getScore(),
		    myPopulation.back().getScore(),
		    (long) std::round((get_timestamp_us() - genStartTime) / 1e3));
	}

	// Must do final score & rank
	scorePopulation();

	rankPopulation();

	this->printEndMsg();
}

// Private methods

void
EvolutionSolver::evolvePopulation()
{

#ifdef _DO_TIMING // For crossing pair section
	const double startTime1 = get_timestamp_us();
#endif

	// Make new generation by discarding unfit
	// genes and mutating/reproducing fit genes
	// to fill up discarded slots

	// First compute possible crossing parents
	msg("Computing cross-able parents: 0%% Done");

	using ParentIdVector = std::vector<std::tuple<IdPair, IdPairs>>;
	ParentIdVector parentIds;

	const uint nPossibleCrossings = (mySruviverCutoff * (mySruviverCutoff + 1)) / 2;
	const uint gaCrossBlock = nPossibleCrossings / 10;

	// Credit to Stack Overflow response by Z boson at
	// http://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
	#pragma omp declare reduction (merge : ParentIdVector : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

	#pragma omp parallel for reduction(merge: parentIds)
	for (int i = 0; i < mySruviverCutoff; i++)
	{
		for (int j = i + 1; j < mySruviverCutoff; j++)
		{
			// If two chromosomes have at least one identical
			// gene and the length sum is legal
			const Chromosome & father =  myPopulation.at(i);
			const Chromosome & mother =  myPopulation.at(j);

			IdPairs crossingIds = father.findCompatibleCrossings(mother);
			if (crossingIds.size() > 0)
				parentIds.push_back(
				    std::make_tuple(
				        std::make_tuple(i, j),
				        crossingIds));

			const uint cid = nPossibleCrossings -
			                 ((mySruviverCutoff - i) * (mySruviverCutoff - i + 1)) / 2
			                 + (j - i);

			if (cid % gaCrossBlock == 0)
			{
				ERASE_LINE();
				msg("Computing cross-able parents: %.2f%% Done",
				    (float) cid / nPossibleCrossings);
			}
		}
	}
	ERASE_LINE();
	msg("Computing cross-able parents: 100%% Done\n");

#ifdef _DO_TIMING
	msg("Section time: %dms\n", (long) ((get_timestamp_us() - startTime1) / 1e3));
#endif

#ifdef _DO_TIMING // For reproduction section
	const double startTime2 = get_timestamp_us();
#endif

	// Probabilistic generation evolution
	msg("Evolution: %.2f%% Done", (float) 0.0f);

	ulong crossCount = 0, pmCount = 0, lmCount = 0, randCount = 0;
	const uint gaPopBlock = myOptions.gaPopSize / 10;

	OMP_PAR_FOR
	for (int i = mySruviverCutoff; i < myOptions.gaPopSize; i++)
	{
		const ulong dice = (ulong) mySruviverCutoff +
		                   std::round(
		                       ((float) std::rand() / RAND_MAX) * myNonSurviverCount);

		// Replicate a high ranking parent
		myPopulation.at(i) = Chromosome(
		                         myPopulation.at(
		                             ((float) std::rand() / RAND_MAX) * mySruviverCutoff
		                         )
		                     );

		if (dice < myCrossCutoff)
		{
			// Pick random parent pair and cross
			if (parentIds.size() > 0)
			{
				const uint randId = std::rand() % parentIds.size();

				const IdPair & parentTuple = std::get<0>(parentIds.at(randId));
				const IdPairs & crossingIds = std::get<1>(parentIds.at(randId));

				uint fatherId, motherId;
				std::tie(fatherId, motherId) = parentTuple;
				const Chromosome & father =  myPopulation.at(fatherId);
				const Chromosome & mother =  myPopulation.at(motherId);

				if (!myPopulation.at(i).cross(father, mother, crossingIds))
					myPopulation.at(i).inheritMutate((std::rand() % 2) ?
					                                 father : mother);
			}
			else
			{
				myPopulation.at(i).inheritMutate(myPopulation.at((std::rand() % myCrossCutoff)));
			}
			crossCount++;
		}
		else if (dice < myPointMutateCutoff)
		{
			if (!myPopulation.at(i).pointMutate())
				myPopulation.at(i).randomise();
			pmCount++;
		}
		else if (dice < myLimbMutateCutoff)
		{
			if (!myPopulation.at(i).limbMutate())
				myPopulation.at(i).randomise();
			lmCount++;
		}
		else
		{
			myPopulation.at(i).randomise();
			randCount++;
		}

		if (i % gaPopBlock == 0)
		{
			ERASE_LINE();
			msg("Evolution: %.2f%% Done", (float) i / myOptions.gaPopSize);
		}
	}

	ERASE_LINE();
	msg("Evolution: 100%% Done\n");

	dbg("Mutation rates: cross %.2f, pm %.2f, lm %.2f, rand %.2f, mySruviverCutoff: %d\n",
	    (float) crossCount / myNonSurviverCount,
	    (float) pmCount / myNonSurviverCount,
	    (float) lmCount / myNonSurviverCount,
	    (float) randCount / myNonSurviverCount,
	    mySruviverCutoff);

#ifdef _DO_TIMING
	msg("Section time: %dms\n", (long) ((get_timestamp_us() - startTime2) / 1e3));
#endif
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)

	std::sort(myPopulation.begin(),
	          myPopulation.end());
}

void
EvolutionSolver::scorePopulation()
{
	msg("Scoring: 0%% Done");
	const uint scoreBlock = myPopulation.size() / 10;

	OMP_PAR_FOR
	for (int i = 0; i < myPopulation.size(); i++)
	{
		myPopulation.at(i).score(mySpec);

		if (i % scoreBlock == 0)
		{
			ERASE_LINE();
			msg("Scoring: %.2f%% Done",
			    (float) i / myPopulation.size());
		}
	}
	ERASE_LINE();
	msg("Scoring: 100%% Done\n");
}

void
EvolutionSolver::initPopulation()
{
	myPopulation = Population(myOptions.gaPopSize);

	const uint block = myOptions.gaPopSize / 10;

	msg("Initialising population: %.2f%% Done", 0.0f);

	OMP_PAR_FOR
	for (int i = 0; i < myOptions.gaPopSize; i++)
	{
		myPopulation.at(i).randomise();
		if (i % block == 0)
		{
			ERASE_LINE();
			msg("Initialising population: %.2f%% Done", (float) i / myOptions.gaPopSize);
		}
	}

	ERASE_LINE();
	msg("Initialising population: 100%% done\n");
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
	    "Population size:            %s\n"
	    "Iterations:                 %s\n"
	    "Survive cutoff:             %u\n"
	    "Cross cutoff:               %u\n"
	    "Point Mutate cutoff:        %u\n"
	    "Limb Mutate cutoff:         %u\n"
	    "New species:                %u\n",
	    psStr.str().c_str(),
	    niStr.str().c_str(),
	    mySruviverCutoff,
	    myCrossCutoff,
	    myPointMutateCutoff,
	    myLimbMutateCutoff,
	    myOptions.gaPopSize - myLimbMutateCutoff);

	OMP_PAR
	{
		if (omp_get_thread_num() == 0)
			msg("Running with %d threads\n", omp_get_num_threads());
	}
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
