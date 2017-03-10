#include <cmath>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <stdlib.h>
#include <unordered_map>

#ifdef _NO_OMP
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

#include "EvolutionSolver.hpp"
#include "util.h"
#include "ParallelUtils.hpp"

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

		evolvePopulation();

		scorePopulation();

		rankPopulation();

		selectParents();

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
EvolutionSolver::selectParents()
{
	TIMING_START(startTimeSelectParents);
	{
		// Ensure variety within survivors using hashmap
		// and crc as key
		using CrcMap = std::unordered_map<Crc32, Chromosome>;
		CrcMap crcMap;
		ulong uniqueCount = 0;

		// We don't want parallelism here because
		// the loop must priotise low indexes
		for (int i = 0; i < myPopulation.size(); i++)
		{
			const Crc32 crc = myPopulation.at(i).checksum();
			if (crcMap.find(crc) == crcMap.end())
			{
				// This individual is a new one - record
				crcMap[crc] = Chromosome(myPopulation.at(i));
				uniqueCount++;

				if (uniqueCount >= mySruviverCutoff)
					break;
			}
		}

		// Insert map-value-indexed individual back into population
		ulong popIndex = 0;
		for (CrcMap::iterator it = crcMap.begin(); it != crcMap.end(); ++it)
		{
			myPopulation.at(popIndex++) = Chromosome(it->second);
		}

		// Sort survivors
		std::sort(myPopulation.begin(),
		          myPopulation.begin() + uniqueCount);
	}
	TIMING_END("selecting", startTimeSelectParents);
}

void
EvolutionSolver::evolvePopulation()
{
	using ParentIdVector = std::vector<std::tuple<IdPair, IdPairs>>;
	ParentIdVector parentIds;

	TIMING_START(startTimeCrossingPairs);
	{
		// Make new generation by discarding unfit
		// genes and mutating/reproducing fit genes
		// to fill up discarded slots

		// First compute possible crossing parents
		msg("Computing cross-able parents: 0%% Done");

		const ulong nPossibleCrossings = (mySruviverCutoff * (mySruviverCutoff + 1)) / 2;
		const ulong gaCrossBlock = nPossibleCrossings / 10;

		// Credit to Stack Overflow response by Z boson at
		// http://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
		#pragma omp declare reduction 										\
		(																	\
		        mergeParentIds : ParentIdVector : 							\
		        omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())	\
		)

		#pragma omp parallel for reduction(mergeParentIds: parentIds)
		for (int i = 0; i < mySruviverCutoff; i++)
		{
			for (int j = i; j < mySruviverCutoff; j++)
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

				const ulong cid = nPossibleCrossings -
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
	}
	TIMING_END("crossing-pairs", startTimeCrossingPairs);


	TIMING_START(startTimeMutation);
	{
		// Probabilistic generation evolution
		msg("Evolution: %.2f%% Done", (float) 0.0f);

		ulong crossCount = 0, pmCount = 0, lmCount = 0, randCount = 0;
		const ulong gaPopBlock = myOptions.gaPopSize / 10;

		OMP_PAR_FOR
		for (int i = mySruviverCutoff; i < myOptions.gaPopSize; i++)
		{
			const ulong dice = mySruviverCutoff +
			                   getDice(myNonSurviverCount);

			if (dice < myCrossCutoff)
			{
				// Pick random parent pair and cross
				if (parentIds.size() > 0)
				{
					const ulong randId = getDice(parentIds.size());

					// Can't use std::tie() because we want references not values
					const IdPair & parentTuple = std::get<0>(parentIds.at(randId));
					const IdPairs & crossingIds = std::get<1>(parentIds.at(randId));

					const Chromosome & father = myPopulation.at(std::get<0>(parentTuple));
					const Chromosome & mother = myPopulation.at(std::get<1>(parentTuple));

					if (!myPopulation.at(i).cross(father, mother, crossingIds))
						myPopulation.at(i).inheritMutate(getDice(2) == 0 ?
						                                 father : mother);
				}
				else
				{
					// Pick a random parent to inherit from and then mutate
					myPopulation.at(i).inheritMutate(
					    myPopulation.at(getDice(myCrossCutoff))
					);
				}
				crossCount++;
			}
			else
			{
				// Replicate a high ranking parent
				myPopulation.at(i) = Chromosome(
				                         myPopulation.at(getDice(mySruviverCutoff))
				                     );

				if (dice < myPointMutateCutoff)
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
					// Individuals not covered by specified mutation
					// rates undergo random destructive mutation
					myPopulation.at(i).randomise();
					randCount++;
				}
			}

			if (i % gaPopBlock == 0)
			{
				ERASE_LINE();
				msg("Evolution: %.2f%% Done", (float) i / myOptions.gaPopSize);
			}
		}

		ERASE_LINE();
		msg("Evolution: 100%% Done\n");

		// Keep some actual counts to make sure the RNG is working
		// correctly
		dbg("Mutation rates: cross %.2f, pm %.2f, lm %.2f, rand %.2f, survivalCount: %d\n",
		    (float) crossCount / myNonSurviverCount,
		    (float) pmCount / myNonSurviverCount,
		    (float) lmCount / myNonSurviverCount,
		    (float) randCount / myNonSurviverCount,
		    mySruviverCutoff);
	}
	TIMING_END("mutation", startTimeMutation);
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)
	TIMING_START(startTimeRanking);
	{
		std::sort(myPopulation.begin(),
		          myPopulation.end());
	}
	TIMING_END("ranking", startTimeRanking);
}

void
EvolutionSolver::scorePopulation()
{
	TIMING_START(startTimeScoring);
	{
		msg("Scoring: 0%% Done");
		const ulong scoreBlock = myPopulation.size() / 10;

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
	TIMING_END("scoring", startTimeScoring);
}

void
EvolutionSolver::initPopulation()
{
	TIMING_START(startTimeInit);
	{

		myPopulation = Population(myOptions.gaPopSize);

		const ulong block = myOptions.gaPopSize / 10;

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
	TIMING_END("init", startTimeInit);
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
	const ulong N = 3;

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
	const ulong minutes = std::floor(timeElapsedInUs / 1e6 / 60.0f);
	const ulong seconds = std::floor(fmod(timeElapsedInUs / 1e6, 60.0f));
	const ulong milliseconds = std::floor(fmod(timeElapsedInUs / 1e3, 1000.0f));
	raw("%um %us %ums\n",
	    minutes, seconds, milliseconds);
}

} // namespace elfin
