#ifndef _PARALLELUTILS_HPP_
#define _PARALLELUTILS_HPP_

#include <vector>
#include <cmath>
#include <omp.h>

#include "../data/PrimitiveShorthands.hpp"

#ifdef _NO_OMP
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_num_devices() { return 0; }
#endif

#define OMP_PAR_FOR _Pragma("omp target teams distribute parallel for simd schedule(runtime)")

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
