#ifndef _PARALLEL_HPP_
#define _PARALLEL_HPP_

#define OMP_PAR _Pragma("omp parallel")
#define OMP_FOR _Pragma("omp for schedule(dynamic)")
#define OMP_PAR_FOR _Pragma("omp parallel for schedule(dynamic)")

namespace elfin
{

template<typename BlockFuncType, typename CBFuncType>
void timedBlockSingleThread(BlockFuncType codeBlock, CBFuncType cb)
{
	const double start_time = get_timestamp_us();
	codeBlock();
	const double end_time = get_timestamp_us();
	cb(end_time - start_time);
}

} // namespace elfin

#endif /* include guard */