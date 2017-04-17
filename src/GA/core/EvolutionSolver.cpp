#include <cmath>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <stdlib.h>
#include <unordered_map>
#include <limits>

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
	mySurviverCutoff = std::round(options.gaSurviveRate * options.gaPopSize);

	myNonSurviverCount = (options.gaPopSize - mySurviverCutoff);
	myCrossCutoff = mySurviverCutoff + std::round(options.gaCrossRate * myNonSurviverCount);
	myPointMutateCutoff = myCrossCutoff + std::round(options.gaPointMutateRate * myNonSurviverCount);
	myLimbMutateCutoff = std::min(
	                         (ulong) (myPointMutateCutoff + std::round(options.gaLimbMutateRate * myNonSurviverCount)),
	                         (ulong) options.gaPopSize);

	myExpectedTargetLen = Chromosome::calcExpectedLength(spec, options.avgPairDist);
	const float minTargetLen = myExpectedTargetLen * (1 - myOptions.chromoLenDev);
	const float maxTargetLen = myExpectedTargetLen * (1 + myOptions.chromoLenDev);

	Chromosome::setup(minTargetLen, maxTargetLen, myRelaMat, myRadiiList);
}

const Population *
EvolutionSolver::population() const
{
	return myCurrPop;
}

const Population &
EvolutionSolver::bestSoFar() const
{
	return myBestSoFar;
}

// Public methods

void
EvolutionSolver::run()
{
	this->printStartMsg();

	this->startTimer();

	initPopulation();

	const int nBestSoFar = 3;
	myBestSoFar.resize(nBestSoFar);

	float lastGenBestScore = std::numeric_limits<double>::infinity();
	int stagnantCount = 0;

	const int genDispDigits = std::ceil(std::log(myOptions.gaIters) / std::log(10));
	char * genMsgFmt;
	asprintf(&genMsgFmt,
	         "Generation #%%%dd: best=%%.2f (%%.2f/module), worst=%%.2f, time taken=%%.0fms\n", genDispDigits);
	char * avgTimeMsgFmt;
	asprintf(&avgTimeMsgFmt,
		"Avg Times: Evolve=%%.0f,Score=%%.0f,Rank=%%.0f,Select=%%.0f,Gen=%%.0f\n");
	for (int i = 0; i < myOptions.gaIters; i++)
	{
		const double genStartTime = get_timestamp_us();

		{
			evolvePopulation();

			scorePopulation();

			rankPopulation();

			selectParents();

			swapPopBuffers();
		}

		const float genBestScore = myCurrPop->front().getScore();
		const ulong genBestChromoLen = myCurrPop->front().genes().size();
		const float genWorstScore = myCurrPop->back().getScore();
		const double genTime = ((get_timestamp_us() - genStartTime) / 1e3);
		msg(genMsgFmt, i,
		    genBestScore,
		    genBestScore / genBestChromoLen,
		    genWorstScore,
		    genTime);
		msg(avgTimeMsgFmt,
		    (float) myTotEvolveTime / (i+1),
                    (float) myTotScoreTime / (i+1),
                    (float) myTotRankTime / (i+1),
                    (float) myTotSelectTime / (i+1),
                    (float) myTotGenTime / (i+1));

		myTotGenTime += genTime;

		// Can stop loop if best score is low enough
		if (genBestScore < myOptions.scoreStopThreshold)
		{
			msg("Score stop threshold %.2f reached\n", myOptions.scoreStopThreshold);
			break;
		}
		else
		{
			for (int i = 0; i < nBestSoFar; i++)
				myBestSoFar.at(i) = myCurrPop->at(i);

			if (float_approximates(genBestScore, lastGenBestScore))
			{
				stagnantCount++;
			}
			else
			{
				stagnantCount = 0;
			}

			lastGenBestScore = genBestScore;

			if (stagnantCount >= myOptions.maxStagnantGens)
			{
				wrn("Solver stopped because max stagnancy is reached (%d)\n", myOptions.maxStagnantGens);
				break;
			}
			else
			{
				msg("Current stagnancy: %d, max: %d\n", stagnantCount, myOptions.maxStagnantGens);
			}
		}
	}

	this->printEndMsg();
}

// Private methods
#ifdef _VTUNE
#include <ittnotify.h>
#endif

void
EvolutionSolver::evolvePopulation()
{
#ifdef _VTUNE
	__itt_resume();  // start VTune, again use 2 underscores
#endif

	CrossingVector possibleCrossings;

	TIMING_START(startTimeEvolving);
	{
		// Probabilistic evolution
		msg("Evolution: %.2f%% Done", (float) 0.0f);

		ulong crossCount = 0, pmCount = 0, lmCount = 0, randCount = 0;
		const ulong gaPopBlock = myOptions.gaPopSize / 10;
		ulong crossFailCount = 0;

		//OMP_PAR_FOR
		Chromosome * myBuffPopData = myBuffPop->data();
		size_t myBuffPopSize = myBuffPop->size();
		const Chromosome * myCurrPopData = myCurrPop->data();
		size_t myCurrPopSize = myCurrPop->size();
		Chromosome & (Chromosome::*assign)(Chromosome const&) = &Chromosome::operator=;

		#pragma omp target teams distribute parallel for simd schedule(runtime) map(myBuffPopData[0:myBuffPopSize], myCurrPopData[0:myCurrPopSize])
		for (int i = 0; i < mySurviverCutoff; i++)
			(myBuffPopData[i].*assign)(myCurrPopData[i]);

		OMP_PAR_FOR
		for (int i = mySurviverCutoff; i < myOptions.gaPopSize; i++)
		{
			Chromosome & chromoToEvolve = myBuffPop->at(i);
			const ulong evolutionDice = mySurviverCutoff +
			                            getDice(myNonSurviverCount);

			if (evolutionDice < myCrossCutoff)
			{
				long motherId, fatherId;
				if (getDice(2))
				{
					motherId = getDice(mySurviverCutoff);
					fatherId = getDice(myOptions.gaPopSize);
				}
				else
				{
					motherId = getDice(myOptions.gaPopSize);
					fatherId = getDice(mySurviverCutoff);
				}

				const Chromosome & mother = myCurrPop->at(motherId);
				const Chromosome & father = myCurrPop->at(fatherId);

				if (!mother.cross(father, chromoToEvolve))
				{
					// Pick a random parent to inherit from and then mutate
					chromoToEvolve = mother.mutateChild();
					crossFailCount++;
				}
				crossCount++;
			}
			else
			{
				// Replicate a high ranking parent
				const ulong parentId = getDice(mySurviverCutoff);
				chromoToEvolve = myCurrPop->at(parentId);

				if (evolutionDice < myPointMutateCutoff)
				{
					if (!chromoToEvolve.pointMutate())
						chromoToEvolve.randomise();
					pmCount++;
				}
				else if (evolutionDice < myLimbMutateCutoff)
				{
					if (!chromoToEvolve.limbMutate())
						chromoToEvolve.randomise();
					lmCount++;
				}
				else
				{
					// Individuals not covered by specified mutation
					// rates undergo random destructive mutation
					chromoToEvolve.randomise();
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
		dbg("Mutation rates: cross %.2f (fail=%d), pm %.2f, lm %.2f, rand %.2f, survivalCount: %d\n",
		    (float) crossCount / myNonSurviverCount,
		    crossFailCount,
		    (float) pmCount / myNonSurviverCount,
		    (float) lmCount / myNonSurviverCount,
		    (float) randCount / myNonSurviverCount,
		    mySurviverCutoff);
	}
	myTotEvolveTime += TIMING_END("evolving", startTimeEvolving);

#ifdef _VTUNE
	__itt_pause(); // stop VTune
#endif
}

void
EvolutionSolver::scorePopulation()
{
	TIMING_START(startTimeScoring);
	{
		msg("Scoring: 0%% Done");
		const ulong scoreBlock = myOptions.gaPopSize / 10;

		OMP_PAR_FOR
		for (int i = 0; i < myOptions.gaPopSize; i++)
		{
			myBuffPop->at(i).score(mySpec);
			if (i % scoreBlock == 0)
			{
				ERASE_LINE();
				msg("Scoring: %.2f%% Done",
				    (float) i / myOptions.gaPopSize);
			}
		}
		ERASE_LINE();
		msg("Scoring: 100%% Done\n");
	}
	myTotScoreTime += TIMING_END("scoring", startTimeScoring);
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)
	TIMING_START(startTimeRanking);
	{
		std::sort(myBuffPop->begin(),
		          myBuffPop->end());
	}
	myTotRankTime += TIMING_END("ranking", startTimeRanking);
}

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
		for (int i = 0; i < myBuffPop->size(); i++)
		{
			const Crc32 crc = myBuffPop->at(i).checksum();
			if (crcMap.find(crc) == crcMap.end())
			{
				// This individual is a new one - record
				crcMap[crc] = myBuffPop->at(i);
				uniqueCount++;

				if (uniqueCount >= mySurviverCutoff)
					break;
			}
		}

		// Insert map-value-indexed individual back into population
		ulong popIndex = 0;
		for (CrcMap::iterator it = crcMap.begin(); it != crcMap.end(); ++it)
			myBuffPop->at(popIndex++) = it->second;

		// Sort survivors
		std::sort(myBuffPop->begin(),
		          myBuffPop->begin() + uniqueCount);
	}
	myTotSelectTime += TIMING_END("selecting", startTimeSelectParents);
}

void
EvolutionSolver::swapPopBuffers()
{
	const Population * tmp = myCurrPop;
	myCurrPop = myBuffPop;
	myBuffPop = const_cast<Population *>(tmp);
}

void
EvolutionSolver::initPopulation()
{
	TIMING_START(startTimeInit);
	{
		myPopulationBuffers[0] = Population(myOptions.gaPopSize);
		myPopulationBuffers[1] = Population(myOptions.gaPopSize);
		myCurrPop = &(myPopulationBuffers[0]);
		myBuffPop = &(myPopulationBuffers[1]);

		const ulong block = myOptions.gaPopSize / 10;

		msg("Initialising population: %.2f%% Done", 0.0f);

		OMP_PAR_FOR
		for (int i = 0; i < myOptions.gaPopSize; i++)
		{
			myBuffPop->at(i).randomise();
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

	// We filled buffer first (because myCurrPop shouldn't be modified)
	swapPopBuffers();
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
	    mySurviverCutoff,
	    myCrossCutoff,
	    myPointMutateCutoff,
	    myLimbMutateCutoff,
	    myOptions.gaPopSize - myLimbMutateCutoff);

	#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			msg("Running with %d threads\n", omp_get_max_threads());
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
		const auto & p = myCurrPop->at(i);
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
