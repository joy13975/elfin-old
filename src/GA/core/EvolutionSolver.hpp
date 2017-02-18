#ifndef _EVOLUTIONSOLVER_HPP_
#define _EVOLUTIONSOLVER_HPP_

#include "../data/TypeDefs.hpp"

namespace elfin
{

class EvolutionSolver
{
public:
	EvolutionSolver(const RelaMat & rm, const Spec & spec);
	virtual ~EvolutionSolver() {};

	void run();

private:
	const RelaMat relaMat;
	const Spec spec;
};

} // namespace elfin

#endif /* include guard */