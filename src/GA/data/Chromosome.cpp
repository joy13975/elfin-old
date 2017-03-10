#include "Chromosome.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include "../core/MathUTils.hpp"
#include "../core/Kabsch.hpp"
#include "../data/PairRelationship.hpp"
#include "../input/JSONParser.hpp"
#include "../core/ParallelUtils.hpp"

namespace elfin
{

// Static variables

bool Chromosome::setupDone = false;
uint Chromosome::myMinLen = 0;
uint Chromosome::myMaxLen = 0;
const RelaMat * Chromosome::myRelaMat = NULL;
const RadiiList * Chromosome::myRadiiList = NULL;
IdPairs Chromosome::myNeighbourCounts;
IdRoulette Chromosome::myGlobalRoulette;

// Constructors
Chromosome::Chromosome()
{
	panic_if(!setupDone, "Chromosome::setup() not called!\n");
}

Chromosome::Chromosome(const Chromosome & rhs)
{
	myGenes = rhs.myGenes;
	myScore = rhs.myScore;
}


Chromosome::Chromosome(const Genes & genes)
{
	myGenes = genes;
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

	// Compute neighbour counts
	const uint dim = myRelaMat->size();
	myNeighbourCounts.resize(dim, std::make_tuple(0, 0));
	for (int i = 0; i < dim; i++)
	{
		uint lhs = 0, rhs = 0;
		for (int j = 0; j < dim; j++)
		{
			if (myRelaMat->at(j).at(i) != NULL)
				lhs++;
			if (myRelaMat->at(i).at(j) != NULL)
				rhs++;
		}

		myNeighbourCounts.at(i) = std::make_tuple(lhs, rhs);
	}

	// Compute global roulette as neighbour count
	myGlobalRoulette.clear();
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < std::get<1>(myNeighbourCounts.at(i)); j++)
			myGlobalRoulette.push_back(i);

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

Genes &
Chromosome::genes()
{
	return myGenes;
}

const Genes &
Chromosome::getGenes() const
{
	return myGenes;
}

std::vector<std::string>
Chromosome::getNodeNames() const
{
	std::vector<std::string> out;

	for (const auto & gene : myGenes)
		out.push_back(Gene::inm->at(gene.nodeId()));

	return out;
}


bool
Chromosome::cross(const Chromosome & father,
                  const Chromosome & mother,
                  const IdPairs & crossingIds)
{
	// cross can fail - if resultant genes collide during synth
	const uint maxTries = 10;
	for (int i = 0; i < maxTries; i++)
	{
		// Pick random crossing point
		const IdPair & crossPoint = crossingIds.at(getDice(crossingIds.size()));
		const uint fatherGeneId = std::get<0>(crossPoint);
		const uint motherGeneId = std::get<1>(crossPoint);

		const Genes & fatherG = father.getGenes();
		const Genes & motherG = mother.getGenes();

		Genes newGenes;
		newGenes.insert(newGenes.end(), fatherG.begin(), fatherG.begin() + fatherGeneId);
		newGenes.insert(newGenes.end(), motherG.begin() + motherGeneId, motherG.end());

		// dbg("Crossing at father[%d] and mother[%d]\n", fatherGeneId, motherGeneId);
		// dbg("Father: \n%s\nMother: %s\n", father.toCString(), mother.toCString());
		// dbg("New Genes: \n%s\n", genesToString(newGenes).c_str());

		if (synthesise(newGenes))
		{
			myGenes = newGenes;
			return true;
		}
	}

	return false;
}

void
Chromosome::inheritMutate(const Chromosome & parent)
{
	*this = Chromosome(parent);
	autoMutate();
}

void
Chromosome::autoMutate()
{
	// Try point mutate first, if not possible then
	// do limb mutate

	if (!pointMutate())
	{
		prf("\npointMutate failed\n");
		if (!limbMutate())
		{
			prf("\nlimbMutate failed\n");
			randomise();
		}
	}
}

std::string
Chromosome::toString() const
{
	std::stringstream ss;

	ss << "Chromosome " << this << ":\n";
	ss << genesToString(myGenes);

	return ss.str();
}

IdPairs
Chromosome::findCompatibleCrossings(const Chromosome & other) const
{
	IdPairs crossingIds;

	const Genes & otherG = other.getGenes();
	const uint otGeneLen = otherG.size();
	const uint myGeneLen = myGenes.size();

	// In below comments gene1 = this, gene2 = other
	for (int i = 0; i < myGeneLen; i++)
	{
		// Using i as gene1 left limb cutoff
		{
			const uint leftLimbLen = i + 1; // This includes the cutoff point
			const uint maxJ = std::min(
			                      std::max(
			                          otGeneLen - (myMinLen - leftLimbLen) - 1, // -1 to account for the duplicate cross point node
			                          (uint) 0),
			                      (uint) otGeneLen - 1);
			const uint minJ = std::min(
			                      std::max(
			                          otGeneLen - (myMaxLen - leftLimbLen) - 1, // -1 to account for the duplicate cross point node
			                          (uint) 0),
			                      (uint) otGeneLen - 1);
			for (int j = minJ; j < maxJ; j++)
			{
				if (myGenes.at(i).nodeId() == otherG.at(j).nodeId())
				{
					// Test 0-i=(left limb), j-end=(right limb)
					const uint childLen = leftLimbLen + (otGeneLen - j - 1);

#ifdef _TEST_CHROMO
					if (childLen < myMinLen || childLen > myMaxLen)
					{
						die("Fatal: length invalid (i=leftLimb) childLen=%d, myMinLen=%d, myMaxLen=%d\n",
						    childLen, myMinLen, myMaxLen);
					}
#endif

					crossingIds.push_back(std::make_tuple(i, j));
				}
			}
		}

		// Using i as gene1 right limb cutoff
		{
			const uint rightLimbLen = myGeneLen - i; // This includes the cutoff point
			const uint minJ = std::min(
			                      std::max(
			                          (myMinLen - rightLimbLen),
			                          (uint) 0),
			                      (uint) otGeneLen - 1);
			const uint maxJ = std::min(
			                      std::max(
			                          (myMaxLen - rightLimbLen),
			                          (uint) 0),
			                      (uint) otGeneLen - 1);
			for (int j = minJ; j < maxJ; j++)
			{
				if (myGenes.at(i).nodeId() == otherG.at(j).nodeId())
				{
					// Test 0-j=(left limb), i-myGenLen=(right limb)
					const uint childLen = j + rightLimbLen;

#ifdef _TEST_CHROMO
					if (childLen < myMinLen || childLen > myMaxLen)
					{
						die("Fatal: length invalid (i=rightLimb) childLen=%d, myMinLen=%d, myMaxLen=%d\n",
						    childLen, myMinLen, myMaxLen);
					}
#endif

					crossingIds.push_back(std::make_tuple(i, j));
				}
			}
		}
	}

	return crossingIds;
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

void
Chromosome::randomise()
{
	do
	{
		myGenes = genRandomGenes(); // Auto synth'd
	}
	while (myGenes.size() < myMinLen);
}

// Private methods

enum PointMutateMode {
	SwapMode,
	InsertMode,
	DeleteMode,
	EnumSize
};

const PointMutateMode pmModeArr[] =
{
	PointMutateMode::SwapMode,
	PointMutateMode::InsertMode,
	PointMutateMode::DeleteMode
};

bool
Chromosome::pointMutate()
{
	// There are 3 ways for point mutate:
	// 1. Swap with another node
	// 2. Insert a node
	// 3. Delete the node
	// As of now it uses equal probability.
	// Could be opened up as a setting.
	const size_t dim = myRelaMat->size();
	const size_t myGeneSize = myGenes.size();
	std::vector<PointMutateMode> modes(pmModeArr,
	                                   pmModeArr + sizeof(pmModeArr) / sizeof(pmModeArr[0]));

	while (modes.size() > 0)
	{
		// Draw a random mode without replacement
		const int modeIndex = getDice(modes.size());
		const PointMutateMode pmMode = modes.at(modeIndex);
		modes.erase(modes.begin() + modeIndex);

		// Try to perform the pointMutate in the chosen mode
		switch (pmMode)
		{
		case PointMutateMode::SwapMode: // Swap mode
		{
			// First int is index from myGenes
			// Second int is nodeId to swap to
			IdPairs swappableIds;
			for (int i = 0; i < myGeneSize; i++)
			{
				// For all neighbours of previous node
				// find those that has nodes[i+1] as
				// one of their RHS neighbours
				for (int j = 0; j < dim; j++)
				{
					// Make sure it's not the original one
					if (j != myGenes.at(i).nodeId())
					{
						// Check whether i can be exchanged for j
						if (
						    (i < myGeneSize) // i can be myGeneSize, which is used for insertion check
						    &&
						    (i == 0 || // Pass if i is the left end
						     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(j) != NULL)
						    &&
						    (i == myGeneSize - 1 || // Pass if i is the right end
						     myRelaMat->at(j).at(myGenes.at(i + 1).nodeId()) != NULL)
						)
						{
							// Make sure resultant shape won't collide with itself
							Genes testGenes(myGenes);
							testGenes.at(i).nodeId() = j;

							// dbg("checking swap at %d/%d of %s\n",
							//     i, myGeneSize, toString().c_str());
							if (synthesise(testGenes))
								swappableIds.push_back(std::make_tuple(i, j));
						}
					}
				}
			}

			// Pick a random one, or report impossible
			if (swappableIds.size() > 0)
			{
				uint swapIndex, swapToNodeId;
				std::tie(swapIndex, swapToNodeId) = swappableIds.at(getDice(swappableIds.size()));
				myGenes.at(swapIndex).nodeId() = swapToNodeId;

				synthesise(myGenes); // This is guaranteed to succeed
				return true;
			}
		}

		case PointMutateMode::InsertMode: // Insert mode
		{
			IdPairs insertableIds;
			if (myGeneSize < myMaxLen)
			{
				for (int i = 0; i < myGeneSize; i++)
				{
					for (int j = 0; j < dim; j++)
					{
						// Check whether j can be inserted before i
						if (
						    (i == 0 || // Pass if inserting at the left end
						     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(j) != NULL)
						    &&
						    (i == myGeneSize || // Pass if appending at the right end
						     myRelaMat->at(j).at(myGenes.at(i).nodeId()) != NULL)
						)
						{
							// Make sure resultant shape won't collide with itself
							Genes testGenes(myGenes);
							testGenes.insert(testGenes.begin() + i, //This is insertion before i
							                 Gene(j));

							// dbg("checking insertion at %d/%d of %s\n",
							//     i, myGeneSize, toString().c_str());
							if (synthesise(testGenes))
								insertableIds.push_back(std::make_tuple(i, j));
						}
					}
				}

				// Pick a random one, or report impossible
				if (insertableIds.size() > 0)
				{
					uint insertAt, insertNodeId;
					std::tie(insertAt, insertNodeId) = insertableIds.at(getDice(insertableIds.size()));
					myGenes.insert(myGenes.begin() + insertAt, //This is insertion before i
					               Gene(insertNodeId));

					synthesise(myGenes); // This is guaranteed to succeed
					return true;
				}
			}
		}

		case PointMutateMode::DeleteMode: // Delete mode
		{
			Ids deletableIds;
			if (myGeneSize > myMinLen)
			{
				for (int i = 0; i < myGeneSize; i++)
				{
					// Check whether i can be deleted
					if (
					    (i < myGeneSize) // i can be myGeneSize, which is used for insertion check
					    &&
					    (i == 0 || i == myGeneSize - 1 || // Pass if i is at either end
					     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(myGenes.at(i + 1).nodeId()) != NULL)
					)
					{
						// Make sure resultant shape won't collide with itself
						Genes testGenes(myGenes);
						testGenes.erase(testGenes.begin() + i);

						// dbg("checking deletion at %d/%d of %s\n",
						//     i, myGeneSize, toString().c_str());
						if (synthesise(testGenes)) deletableIds.push_back(i);
					}
				}

				// Pick a random one, or report impossible
				if (deletableIds.size() > 0)
				{
					myGenes.erase(myGenes.begin() + deletableIds.at(getDice(deletableIds.size())));

					synthesise(myGenes); // This is guaranteed to succeed
					return true;
				}
			}
		}

		default:
		{
			// Fell through all cases without mutating
			// Do nothing unless pmMode is strange
			panic_if(pmMode < 0 || pmMode >= PointMutateMode::EnumSize,
			         "Invalid pmMode in Chromosome::pointMutate()\n");
		}
		}
	}

	return false;
}

bool
Chromosome::limbMutate()
{
	const size_t N = myGenes.size();

	// Pick a node that can host an alternative limb
	uint severId = -1;
	bool mutateLeftLimb = false;

	// Try max 10 times
	const int maxTries = 10;
	for (int i = 0; i < maxTries; i++)
	{
		const uint geneId = getDice(N - 1) + 1;
		const uint nodeId = myGenes.at(geneId).nodeId();

		uint lhs, rhs;
		std::tie(lhs, rhs) = myNeighbourCounts.at(nodeId);

		if (lhs == 1 && rhs == 1)
			continue;

		if (lhs == 1)
			mutateLeftLimb = false;
		else if (rhs == 1)
			mutateLeftLimb = true;
		else
			mutateLeftLimb = getDice(2);

		severId = geneId;
		break;
	}

	if (severId == -1)
		return false;

	// Server the limb
	if (mutateLeftLimb)
		myGenes.erase(myGenes.begin(), myGenes.begin() + severId);
	else
		myGenes.erase(myGenes.begin() + severId + 1, myGenes.end());

	const uint severedLen = N - myGenes.size();

	// Re-generate that whole "limb"
	Genes newGenes;
	for (int i = 0; i < maxTries; i++)
	{
		newGenes = mutateLeftLimb ?
		           genRandomGenesReverse(myMaxLen, myGenes) :
		           genRandomGenes(myMaxLen, myGenes);

		if (newGenes.size() >= myMinLen)
			break;
	}

	if (newGenes.size() < myMinLen)
		return false;

	myGenes = newGenes;

	return true;
}

Genes
Chromosome::genRandomGenesReverse(
    const uint genMaxLen,
    Genes genes)
{
	const size_t dim = myRelaMat->size();

	if (genes.size() == 0)
	{
		// Pick random starting node
		const uint firstNodeId = myGlobalRoulette.at(getDice(myGlobalRoulette.size()));
		genes.emplace_back(firstNodeId, 0, 0, 0);
	}
	else
	{
		synthesiseReverse(genes);
	}

	// Reverse order so growth tip is at back
	std::reverse(genes.begin(), genes.end());

	while (genes.size() < genMaxLen)
	{
		std::vector<uint> rouletteWheel;
		const Gene & currGene = genes.back();

		// Compute whether each neighbour is colliding
		for (int i = 0; i < dim; i++)
		{
			if (i == currGene.nodeId())
				continue;

			const PairRelationship * prPtr = myRelaMat->at(i).at(currGene.nodeId());
			if (prPtr == NULL)
				continue;

			const Point3f checkpoint = prPtr->tran;

			// Create roulette based on number of neighbours
			// of the current neighbour being considered
			if (!collides(i,
			              checkpoint,
			              genes.begin(),
			              genes.end() - 2,
			              *myRadiiList))
				for (int j = 0; j < std::get<0>(myNeighbourCounts.at(i)); j++)
					rouletteWheel.push_back(i);
		}

		if (rouletteWheel.size() == 0)
			break;

		// Pick a random valid neighbour
		const uint nextNodeId = rouletteWheel.at(getDice(rouletteWheel.size()));

		const PairRelationship * nextNodePR = myRelaMat->at(nextNodeId).at(currGene.nodeId());

		// Grow shape
		for (auto & g : genes)
		{
			g.com() -= nextNodePR->tran;
			g.com() = g.com().dot(nextNodePR->rotInv);
		}

		genes.emplace_back(nextNodeId, 0, 0, 0);
	}

	// Reverse the reverse!
	std::reverse(genes.begin(), genes.end());

	return genes;
}

/*
 * Generate a sequence of non-colliding nodes, with an
 * optional starting genes object and a maximum allowed
 * length. Minimum length should be checked outside this
 * function as the randomness might fail to produce
 * a sufficiently long chain.
 */
Genes
Chromosome::genRandomGenes(
    const uint genMaxLen,
    Genes genes)
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

	if (genes.size() == 0)
	{
		// Pick random starting node
		const uint firstNodeId = myGlobalRoulette.at(getDice(myGlobalRoulette.size()));
		genes.emplace_back(firstNodeId, 0, 0, 0);
	}
	else
	{
		synthesise(genes);
	}

	while (genes.size() < genMaxLen)
	{
		std::vector<uint> rouletteWheel;
		const Gene & currGene = genes.back();

		// Compute whether each neighbour is colliding
		const std::vector<PairRelationship *> rr = myRelaMat->at(currGene.nodeId());
		for (int i = 0; i < dim; i++)
		{
			if (i == currGene.nodeId())
				continue;
			const PairRelationship * prPtr = rr.at(i);

			// Create roulette based on number of neighbours
			// of the current neighbour being considered
			if (prPtr && !collides(i,
			                       prPtr->comB,
			                       genes.begin(),
			                       genes.end() - 2,
			                       *myRadiiList))
				for (int j = 0; j < std::get<1>(myNeighbourCounts.at(i)); j++)
					rouletteWheel.push_back(i);
		}

		if (rouletteWheel.size() == 0)
			break;

		// Pick a random valid neighbour
		const uint nextNodeId = rouletteWheel.at(getDice(rouletteWheel.size()));

		const PairRelationship * nextNodePR = rr.at(nextNodeId);

		// Grow shape
		for (auto & g : genes)
		{
			g.com() = g.com().dot(nextNodePR->rot);
			g.com() += nextNodePR->tran;
		}

		genes.emplace_back(nextNodeId, 0, 0, 0);
	}

	return genes;
}

bool
Chromosome::synthesiseReverse(Genes & genes)
{
	const uint N = genes.size();
	if (N == 0)
		return true;

	// Zero all coords so they don't interfere with
	// collision check before being synth'd
	for (auto & g : genes)
		g.com().x = (g.com().y = (g.com().z = 0));

	for (int i = N - 1; i > 0; i--)
	{
		const auto & lhsGene = genes.at(i - 1);
		const auto & rhsGene = genes.at(i);

		// Check collision
		const PairRelationship * newNodePr =
		    myRelaMat->at(lhsGene.nodeId()).at(rhsGene.nodeId());

		if (newNodePr == NULL)
		{
			// Fatal failure; diagnose!
			err("Synthesise(): impossible pair! %d(%s) <-x-> %d(%s)\n",
			    lhsGene.nodeId(),
			    Gene::inm->at(lhsGene.nodeId()).c_str(),
			    rhsGene.nodeId(),
			    Gene::inm->at(rhsGene.nodeId()).c_str());

			err("Erroneous genes:\n%s\n", genesToString(genes).c_str());


			die("Fatal error in synthesiseReverse(): should never use impossible pair\n");
		}

		const Point3f checkpoint = newNodePr->tran;

		if (collides(lhsGene.nodeId(),
		             checkpoint,
		             genes.begin() + i + 2,
		             genes.end(),
		             *myRadiiList))
			return false;

		// Grow shape
		for (int j = N - 1; j > i - 1; j--)
		{
			auto & g = genes.at(j);
			g.com() -= newNodePr->tran;
			g.com() = g.com().dot(newNodePr->rotInv);
		}
	}

	return true;
}


bool
Chromosome::synthesise(Genes & genes)
{
	if (genes.size() == 0)
		return true;

	// Zero all coords so they don't interfere with
	// collision check before being synth'd
	for (auto & g : genes)
		g.com().x = (g.com().y = (g.com().z = 0));

	for (int i = 1; i < genes.size(); i++)
	{
		const auto & lhsGene = genes.at(i - 1);
		const auto & rhsGene = genes.at(i);

		// Check collision
		const PairRelationship * newNodePr =
		    myRelaMat->at(lhsGene.nodeId()).at(rhsGene.nodeId());

		if (newNodePr == NULL)
		{
			// Fatal failure; diagnose!
			err("Synthesise(): impossible pair! %d(%s) <-x-> %d(%s)\n",
			    lhsGene.nodeId(),
			    Gene::inm->at(lhsGene.nodeId()).c_str(),
			    rhsGene.nodeId(),
			    Gene::inm->at(rhsGene.nodeId()).c_str());

			err("Erroneous genes:\n%s\n", genesToString(genes).c_str());

			die("Fatal error in synthesise(): should never use impossible pair\n");
		}

		if (collides(rhsGene.nodeId(),
		             newNodePr->comB,
		             genes.begin(),
		             genes.begin() + i - 2,
		             *myRadiiList))
			return false;

		// Grow shape
		for (int j = 0; j < i; j++)
		{
			auto & g = genes.at(j);
			g.com() = g.com().dot(newNodePr->rot);
			g.com() += newNodePr->tran;
		}
	}

	return true;
}

} // namespace elfin

#ifdef _TEST_CHROMO

int main(int argc, const char ** argv)
{
	using namespace elfin;

	// Load necessary data to setup Gene
	RelaMat relaMat;
	NameIdMap nameIdMap;
	IdNameMap idNameMap;
	RadiiList radiiList;
	JSONParser().parseDB("../../res/xDB.json", nameIdMap, idNameMap, relaMat, radiiList);

	Gene::setup(&idNameMap);
	Chromosome::setup(0, 100, relaMat, radiiList);

	std::string l10Test1NameArr[] = {
		"D54",
		"D54_j1_D79",
		"D79_j2_D14",
		"D14_j1_D79",
		"D79_j1_D54",
		"D54_j1_D79",
		"D79_j1_D54",
		"D54",
		"D54",
		"D54_j1_D79",
		"D79",
		"D79",
		"D79_j1_D54",
		"D54_j1_D79",
		"D79",
		"D79",
		"D79_j1_D54",
		"D54_j1_D79",
		"D79",
		"D79_j1_D54",
		"D54"
	};

	const Point3f l10Solution1Arr[]
	{
		Point3f(-5.328838348, -102.4714661, -141.4473877),
		Point3f(-3.610388756, -131.6208344, -156.3851471),
		Point3f(-0.08456516266, -162.524353, -180.9295959),
		Point3f(-43.3194313, -181.9328766, -195.2945099),
		Point3f(-83.6757431, -184.3944855, -196.3811035),
		Point3f(-98.81578827, -214.2256165, -167.1357422),
		Point3f(-113.0508194, -239.6369324, -140.295639),
		Point3f(-98.62628937, -233.7519531, -114.0887299),
		Point3f(-87.21094513, -224.6260529, -101.7331085),
		Point3f(-65.63604736, -209.1479645, -82.47885895),
		Point3f(-50.25890732, -195.6275482, -55.12817383),
		Point3f(-29.9008255, -189.4302826, -47.3944397),
		Point3f(-9.871655464, -172.0654907, -57.61371613),
		Point3f(-24.33072662, -130.3964539, -63.00215149),
		Point3f(-30.3760376, -98.14616394, -72.52536011),
		Point3f(-46.71737671, -82.48023987, -72.98136902),
		Point3f(-63.45676422, -71.50189972, -52.82184601),
		Point3f(-40.53733063, -51.61549759, -20.36203384),
		Point3f(-28.71621323, -38.23616409, 8.768226624),
		Point3f(-4.185769081, -25.74235153, 15.7895956),
		Point3f(0, 0, 0)
	};

	const Points3f l10Solution1(l10Solution1Arr,
	                            l10Solution1Arr +
	                            sizeof(l10Solution1Arr) / sizeof(l10Solution1Arr[0]));

	// Parse name into ID, then gene vector
	Genes genes;
	for (const auto & name : l10Test1NameArr)
		genes.emplace_back(nameIdMap.at(name), 0, 0, 0);

	Chromosome chromo(genes);

	// Test synthesis
	uint failCount = 0;

	if (!chromo.synthesise(chromo.getGenes()))
	{
		failCount++;
		wrn("Could not synthesise known spec!\n");
	}
	else
	{
		// Check that exact coordinates match
		for (int i = 0; i < chromo.getGenes().size(); i++)
		{
			if (!chromo.getGenes().at(i).com().approximates(l10Solution1.at(i)))
				failCount++;
		}
	}

	msg("%s\n", chromo.toString().c_str());

	// Test verdict
	if (failCount == 0)
		msg("Passed!\n");
	else
		err("Failed! failCount=%d\n", failCount);

	return 0;
}

#endif //_TEST_CHROMO