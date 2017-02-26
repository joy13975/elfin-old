#ifndef _KABSCH_HPP_
#define _KABSCH_HPP_

#include <vector>

#include "../data/TypeDefs.hpp"

namespace elfin
{
bool Kabsch(
    Points3f const & mobile,
    Points3f const & ref,
    Matrix<double> & rot,
    Vector3f & tran,
    double & rms);

bool RosettaKabsch(
    std::vector < std::vector < double > > const & x,
    std::vector < std::vector < double > > const & y,
    int const n,
    int const mode,
    double *rms,
    std::vector < double > & t,
    std::vector < std::vector < double > > & u );
} // namespace elfin

#endif /* include guard */