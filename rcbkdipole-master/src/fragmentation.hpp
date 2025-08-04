#ifndef _FRAGMENTATION_H
#define _FRAGMENTATION_H

/*
 * Fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "config.hpp"

/**
 * Virtual fragmentation function class.
 *
 * All fragmentation functions are inherited from this and implement
 * Evaluate() function
 */

class FragmentationFunction
{
    public:
        // D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
        /**
         * Evaluate fragmentation function
         *
         * Returns D(z) for parton->hadron fragmentation (note: not z*D(z))
         * @param p Parton
         * @param h Hadron
         * @param z Longitudinal momentun fraction
         * @param qs Scale [GeV]
         */
        virtual double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double z, double qs)=0;
        FragmentationFunction();
        virtual std::string GetString();
        Amplitude::Order GetOrder();

        /**
         * Select between LO and NLO distributions.
         *
         * Check that the FF you use actually supports both NLO and LO!
         */
        virtual void SetOrder(Amplitude::Order order_);
        
        virtual void Test();
    protected:
        Amplitude::Order order; //! Select NLO or LO fragmentation function

};


std::string ParticleStr(Amplitude::Hadron h);

#endif
