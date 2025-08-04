#ifndef _KKP_H
#define _KKP_H

/*
 * AmplitudeLib PDF
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "config.hpp"
#include "fragmentation.hpp"

/**
 * KKP fragmentation function.
 * Uses fragmentation_kkp.f downloaded from
 * http://www.desy.de/~poetter/kkp.html
 */
class KKP : public FragmentationFunction
{
    public:
        double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double x, double qs);
        KKP();
        std::string GetString();
        
        void Test();
    private:


};


#endif
