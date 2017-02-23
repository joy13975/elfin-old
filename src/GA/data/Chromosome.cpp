#include "Chromosome.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "../core/MathUTils.hpp"
#include "../data/PairRelationship.hpp"

namespace elfin
{

// Constructors
Chromosome::Chromosome()
{}

Chromosome::Chromosome(const Genes & geneRef) :
	myGenes(geneRef)
{}

// Public methods

void Chromosome::score(const Points3f & ref)
{
	// Use Kabsch to move genes
	superimpose(ref);

	// Use sum of distances - not different effect from RMSD
	myScore = 0.0f;

	for (auto & g : myGenes)
		myScore += minDistFromLine(g.com, ref);
}

/*
 * Calculate expected length as total point
 * displacements over avg pair module distance
 *
 * Note: another possible heuristic is to insert
 * 	an extra point count only between pairs of
 * 	points that are too far apart, i.e. 2x avg dist
 */
uint Chromosome::calcExpectedLength(const Points3f & lenRef,
                                    const float avgPairDist)
{
	float sumDist = 0.0f;

	for (std::vector<Point3f>::const_iterator i = lenRef.begin() + 1; // !!
	        i != lenRef.end();
	        ++i)
		sumDist += norm(*(i - 1), *i);

	return (uint) round(sumDist / avgPairDist);
}

/*
 * Build a random list of genes that is valid
 * i.e. non colliding and conforms to relaMat
 *
 * The final length is random - accept as long
 * as it is within [minLen, maxLen]
 */
Genes Chromosome::genRandomGenes(
    const uint minLen,
    const uint maxLen,
    const RelaMat & relaMat,
    const RadiiList & radiiList)
{
	const size_t dim = relaMat.size();
	double seed = get_timestamp_us();
	std::srand(seed);

	Genes genes;
	std::vector<uint> rouletteWheel;

	do
	{
		genes.clear();
		rouletteWheel.clear();
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (relaMat.at(i).at(j))
					rouletteWheel.push_back(i);
		dbg("Random generation try\n");

		// Pick random starting node
		const uint firstNodeId = rouletteWheel.at(
		                             std::rand() % rouletteWheel.size());
		genes.emplace_back(firstNodeId, 0, 0, 0);

		while (genes.size() < maxLen)
		{
			rouletteWheel.clear();
			const Gene currGene = genes.back();

			// Compute whether each neighbour is colliding
			const RelaRow & rr = relaMat.at(currGene.nodeId);
			prf("\nWheel for node=%d: \n", currGene.nodeId);
			for (int i = 0; i < dim; i++)
			{
				const PairRelationship * prPtr = rr.at(i);

				// Create roulette based on number of neighbours
				// of the current neighbour being considered
				if (prPtr &&
				        !collides(i,
				                  prPtr->comB,
				                  genes,
				                  radiiList))
					for (int j = 0; j < dim; j++)
						if (relaMat.at(i).at(j))
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
			prf("Pick node: %d\n", nextNodeId);

			const PairRelationship * nextNodePR = rr.at(nextNodeId);

			// Grow shape
			for (auto & g : genes)
			{
				rotate(g.com, nextNodePR->rot);
				translate(g.com, nextNodePR->tran);
			}

			genes.emplace_back(nextNodeId, 0, 0, 0);
		}
	}
	while (genes.size() < minLen);

	return genes;
}

// Private meth

void Chromosome::superimpose(const Points3f & ref)
{
	wrn("TODO: superimpose obj to reference using Kasbch\n");
}

} // namespace elfin