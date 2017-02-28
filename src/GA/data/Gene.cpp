#include "Gene.hpp"

namespace elfin
{

bool Gene::setupDone = false;
const IdNameMap * Gene::inm = NULL;

uint &
Gene::nodeId()
{
	return myNodeId;
}

const uint &
Gene::nodeId() const
{
	return myNodeId;
}

Point3f &
Gene::com()
{
	return myCom;
}

const Point3f &
Gene::com() const
{
	return myCom;
}


Gene::Gene(const uint _nodeId,
           const Point3f _com) :
	myNodeId(_nodeId),
	myCom(_com)
{
	panic_if(!setupDone,
	         "Gene::setup() must be callsed first!\n");
}

Gene::Gene(const uint _nodeId,
           const float x,
           const float y,
           const float z) :
	myNodeId(_nodeId),
	myCom(x, y, z)
{
	panic_if(!setupDone,
	         "Gene::setup() must be callsed first!\n");
}

std::string
Gene::toString() const
{
	std::stringstream ss;
	ss << "ID: " << myNodeId << " ";
	ss << "Name: " << inm->at(myNodeId) << " ";
	ss << "CoM: " << myCom.toString();

	return ss.str();
}

void
Gene::setup(const IdNameMap * _inm)
{
	inm = _inm;
	setupDone = true;
}
} // namespace elfin