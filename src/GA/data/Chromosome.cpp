#include "Chromosome.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include "../core/MathUTils.hpp"
#include "../core/Kabsch.hpp"
#include "../data/PairRelationship.hpp"

namespace elfin
{

bool Chromosome::setupDone = false;
uint Chromosome::myMinLen = 0;
uint Chromosome::myMaxLen = 0;
const RelaMat * Chromosome::myRelaMat = NULL;
const RadiiList * Chromosome::myRadiiList = NULL;

// Constructors
Chromosome::Chromosome()
{
	panic_if(!setupDone, "Chromosome::setup() not called!\n");
}

void
Chromosome::setup(const uint minLen,
                  const uint maxLen,
                  const RelaMat & relaMat,
                  const RadiiList & radiiList)
{
	if (setupDone)
		die("Chromosome::setup() called second time!\n");

	myMinLen = minLen;
	myMaxLen = maxLen;
	myRelaMat = &relaMat;
	myRadiiList = &radiiList;

	setupDone = true;
}

// Public methods

void
Chromosome::score(const Points3f & ref)
{
	myScore = kabschScore(myGenes, ref);
	// Use Kabsch to move genes
	// superimpose(ref);

	// Use sum of distances - not different effect from RMSD
	// myScore = 0.0f;

	// for (auto & g : myGenes)
	// 	myScore += minDistFromLine(g.com(), ref);
}

bool
Chromosome::operator>(const Chromosome & rhs) const
{
	return myScore > rhs.getScore();
}

bool
Chromosome::operator<(const Chromosome & rhs) const
{
	return myScore < rhs.getScore();
}

float
Chromosome::getScore() const
{
	return myScore;
}

void
Chromosome::mutate() 
{
	// Try point mutate first, if not possible then 
	// do limb mutate

	if(!pointMutate())
		limbMutate();
}

std::string
Chromosome::toString() const
{
	std::stringstream ss;
	ss << "Chromosome " << this << ":\n";

	const int N = myGenes.size();
	for (int i = 0; i < N; i++)
	{
		ss << "Node #" << i << " / " << N << ": "
		   << myGenes.at(i).toString() << std::endl;
	}

	return ss.str();
}

/*
 * Calculate expected length as total point
 * displacements over avg pair module distance
 *
 * Note: another possible heuristic is to insert
 * 	an extra point count only between pairs of
 * 	points that are too far apart, i.e. 2x avg dist
 */
uint
Chromosome::calcExpectedLength(const Points3f & lenRef,
                               const float avgPairDist)
{
	float sumDist = 0.0f;

	for (std::vector<Point3f>::const_iterator i = lenRef.begin() + 1; // !!
	        i != lenRef.end();
	        ++i)
		sumDist += (i - 1)->distTo(*i);

	return (uint) round(sumDist / avgPairDist);
}

/*
 * Build a random list of genes that is valid
 * i.e. non colliding and conforms to relaMat
 *
 * The final length is random - accept as long
 * as it is within [minLen, maxLen]
 *
 * Note: for performance, might need to strip
 * collision checking from this function
 */
void
Chromosome::randomise()
{
	const size_t dim = myRelaMat->size();

	// A roulette wheel represents the probability of
	// each node being picked as the next node, based
	// on the number of neighbours they have.
	//
	// This is so as to not pick "dead-end" nodes
	// often, which can result in very "boring" shapes
	// e.g. formed by repetition of just one node
	//
	// Once a dead-end node is picked, further nodes
	// are simply repeated and we don't want this to
	// happen often.
	std::vector<uint> rouletteWheel;

	do
	{
		myGenes.clear();
		rouletteWheel.clear();
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (myRelaMat->at(i).at(j))
					rouletteWheel.push_back(i);

		// Pick random starting node
		const uint firstNodeId = rouletteWheel.at(
		                             std::rand() % rouletteWheel.size());
		myGenes.emplace_back(firstNodeId, 0, 0, 0);

		while (myGenes.size() < myMaxLen)
		{
			rouletteWheel.clear();
			const Gene currGene = myGenes.back();

			// Compute whether each neighbour is colliding
			const RelaRow & rr = myRelaMat->at(currGene.nodeId());
			for (int i = 0; i < dim; i++)
			{
				const PairRelationship * prPtr = rr.at(i);

				// Create roulette based on number of neighbours
				// of the current neighbour being considered
				if (prPtr &&
				        !collides(i,
				                  prPtr->comB,
				                  myGenes,
				                  *myRadiiList))
					for (int j = 0; j < dim; j++)
						if (myRelaMat->at(i).at(j))
						{
							raw_at(LOG_PROOF, "%d, ", i);
							rouletteWheel.push_back(i);
						}
			}

			if (rouletteWheel.size() == 0)
				break;

			// Pick a random valid neighbour
			const uint nextNodeId = rouletteWheel.at(std::rand() %
			                        rouletteWheel.size());

			const PairRelationship * nextNodePR = rr.at(nextNodeId);

			// Grow shape
			for (auto & g : myGenes)
			{
				g.com().rotate(nextNodePR->rot);
				g.com() += nextNodePR->tran;
			}

			myGenes.emplace_back(nextNodeId, 0, 0, 0);
		}
	}
	while (myGenes.size() < myMinLen);
}

// Private methods

bool
Chromosome::pointMutate()
{
	const size_t dim = myRelaMat->size();

	// Find nodes that can be swapped
	// First uint is index from myGenes
	// Second uint is nodeId to swap to
	std::vector<std::tuple<uint, uint>> swappableIds;

	for (int i = 1; i < myGenes.size() - 1; i++)
	{
		// For all neighbours of previous node
		// find those that has nodes[i+1] as
		// one of their RHS neighbours
		for (int j = 0; j < dim; j++)
			if (j != myGenes.at(i).nodeId() && // Make sure it's not the original one
			        myRelaMat->at(myGenes.at(i - 1).nodeId()).at(j) != NULL &&
			        myRelaMat->at(j).at(myGenes.at(i + 1).nodeId()) != NULL)
				swappableIds.push_back(std::make_tuple(i, j));
	}

	// Pick a random one, or report impossible
	if (swappableIds.size() > 0)
	{
		uint swapIndex, swapToNodeId;
		std::tie(swapIndex, swapToNodeId) = swappableIds.at(std::rand() % swappableIds.size());
		myGenes.at(swapIndex).nodeId() = swapToNodeId;
		syntehsise();
		return true;
	}
	else
	{
		return false;
	}
}

void
Chromosome::limbMutate()
{
	// Pick a random direction and node to start

	// Re-generate that whole "limb"
	// wrn("TODO: limbMutate()\n");

}

void
Chromosome::syntehsise()
{
	if (myGenes.size() == 0)
		return;

	auto & firstGene = myGenes.front();
	firstGene.com().x = (firstGene.com().y = (firstGene.com().z = 0));

	for (int i = 1; i < myGenes.size(); i++)
	{
		const auto & prevGene = myGenes.at(i - 1);
		auto & currGene = myGenes.at(i);

		// Current node is always the tip at origin
		currGene.com().x = (currGene.com().y = (currGene.com().z = 0));

		// Grow shape
		const PairRelationship * nextNodePR =
		    myRelaMat->at(prevGene.nodeId()).at(currGene.nodeId());
		for (int j = 0; j < i; j++)
		{
			auto & g = myGenes.at(j);
			g.com().rotate(nextNodePR->rot);
			g.com() += nextNodePR->tran;
		}
	}
}

} // namespace elfin