#ifndef _PARALLELUTILS_HPP_
#define _PARALLELUTILS_HPP_

#include <vector>
#include <cmath>

#include "../data/PrimitiveShorthands.hpp"

// Macro shorthands
#define OMP_PAR _Pragma("omp parallel")

#ifdef _NO_OMP
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

// Use dynamic by default because for elfin,
// workloads can vary by a large margin
// as there are many different branches leading
// to different amount of work (though untested)
#define OMP_FOR _Pragma("omp for")
#define OMP_PAR_FOR _Pragma("omp parallel for")

#ifdef _DO_TIMING

#include "util.h"
#define TIMING_START(varName) \
	const double varName = get_timestamp_us();
#define TIMING_END(sectionName, varName) \
	msg("Section (%s) time: %dms\n", sectionName, (long) ((get_timestamp_us() - varName) / 1e3));

#else //ifdef _DO_TIMING

#define TIMING_START(varName)
#define TIMING_END(sectionName, varName)

#endif //ifdef _DO_TIMING


namespace elfin
{

void setupParaUtils(uint globalSeed);

std::vector<uint> & getParaRandSeeds();

inline ulong getDice(ulong ceiling)
{
	return (ulong) std::floor(
	           (
	               (float)
	               (ceiling - 1) *
	               rand_r(
	                   &(getParaRandSeeds().at(omp_get_thread_num()))
	               )
	               / RAND_MAX
	           )
	       );
}
} // namespace elfin

#endif /* include guard */