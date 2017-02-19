#include "Gene.hpp"

#include <algorithm>

namespace elfin
{

// Constructors
Gene::Gene()
{}

// Public methods

void Gene::score(const Points3f & ref)
{
	wrn("TODO: score function\n");

	// Grow actual shape using chromosome
	Points3f shape = growShape();

	// Use Kabsch to move
	superimpose(ref, shape);

	// RMSD
	myScore = calcRMSD(ref, shape);
}

uint Gene::calcExpectedLength(const Points3f & lenRef)
{
	// Calculate expected length as total point
	// displacements over avg pair module distance

	wrn("TODO: calculate expected gene length\n");

	return 0;
}

Gene::Gene(const Chromosome & chromoRef) :
	myChromo(chromoRef)
{}

Chromosome Gene::randomChromosome(
    const uint expLen,
    const float lenDev,
    const RelaMat & relaMat)
{
	wrn("TODO: random gene generation, and seeding\n");

	Chromosome chromo;

	return chromo;
}

// Private methods

float Gene::calcRMSD(const Points3f & ref, Points3f & toMove)
{
	wrn("TODO: calculate RMSD using line nearest normal method\n");

	return NAN;
}

void Gene::superimpose(const Points3f & ref, Points3f & toMove)
{
	wrn("TODO: superimpose obj to reference using Kasbch\n");
}

Points3f Gene::growShape()
{
	Points3f shape;

	wrn("TODO: grow shape using chromosome\n");

	return shape;
}

} // namespace elfin