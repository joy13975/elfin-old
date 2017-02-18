#include "EvolutionSolver.hpp"
#include "util.h"

namespace elfin
{

EvolutionSolver::EvolutionSolver(const RelaMat & _relaMat, const Spec & _spec) :
	relaMat(_relaMat), spec(_spec)
{

}

void EvolutionSolver::run()
{
	msg("elfin EvolutionSolver.run()\n");



	msg("elfin EvolutionSolver finished\n");
}

} // namespace elfin
