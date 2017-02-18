#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <vector>
#include <map>

namespace elfin
{

class PairRelationship;

typedef std::vector<PairRelationship *> RelaRow;
typedef std::vector<RelaRow> RelaMat;
typedef std::map<std::string, unsigned int> NameIdMap;

typedef std::vector<float> Spec;

} // namespace elfin

#endif /* include guard */