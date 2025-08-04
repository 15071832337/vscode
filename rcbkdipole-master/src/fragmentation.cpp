/*
 * Virtual class to hide different fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include <sstream>
#include "config.hpp"
#include "fragmentation.hpp"


using namespace Amplitude;


FragmentationFunction::FragmentationFunction()
{
    order = NLO;

}

std::string FragmentationFunction::GetString()
{
    return "not specified";
}

void FragmentationFunction::SetOrder(Order o)
{
    order = o;
}

Order FragmentationFunction::GetOrder()
{
    return order;
}

void FragmentationFunction::Test()
{
	cerr << "FragmentationFunction::Test() is not implemented" << endl;
}


std::string ParticleStr(Hadron h)
{
	switch(h)
	{
		case P:
			return "proton";
		case PI0:
			return "pi0";
		case H:
			return "charged hadron";
		case HM:
			return "h-";
		case HP:
			return "h+";
		case PIP:
			return "pi+";
		case PIM:
			return "pi-";
		default:
			std::stringstream ss;
			ss << "unkown particle type " << h;
			return ss.str();
	}
	
	return "ERROR!";
	
}
