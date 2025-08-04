/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

/** @file */

#ifndef _CONFIG_HPP
#define _CONFIG_HPP 


#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#include <string>
using std::string;
#include <cmath>

typedef unsigned int uint;


// Constants in Amplitude namespace, avoid overlapping with other programs
// using this class

namespace Amplitude
{

    const double ALPHA_e = 1.0/137.035999679; 
    //const double e = sqrt(4.0*M_PI*ALPHA_e);

    const double FMGEV = 5.068;



    // Other constants

    const double eps=0.000001;

    // Inline functions

    #define LINEINFO __FILE__ << ":" << __LINE__

    inline double SQR2(const double x) { return x*x; }

    enum Parton
    {
        UVAL,   ///< valence u
        DVAL,   ///< valence d
        USEA,   ///< sea u
        DSEA,   ///< sea d
        U,      ///< u (valence + sea)
        D,      ///< d (valence + sea)
        S,      ///< strange
        C,      ///< charm
        B,      ///< bottom
        G,      ///< gluon
        UBAR,   ///< anti-u
        DBAR,   ///< anti-d
        SBAR,   ///< anti-s
		CBAR,   ///< anti-c
		BBAR,   ///< anti-b
        LIGHT	///< light quarks = u+d+s
    };

    /**
     * Selection between NLO and LO distributions
     */
    enum Order
    {
        LO,
        NLO
    };

    enum RunningAlphas
    {
        RUNNING,
        FIXED
    };

    /**
     * Representation at which the dipole amplitude is evaluated
     */
    enum Representation
    {
        FUNDAMENTAL,
        ADJOINT
    };

    enum Hadron
    {
        PI,   ///< pi+, pi- (charged pion)
        PIP,  ///< pi+
        PIM,  ///< pi-
        PI0,  ///< pi0 (neutral pion)
        K,    ///< K+, K- (charged kaon)
        KP,   ///< K+
        KM,   ///< K-
        K0,   ///< K0, bar K0 (neutral kaon)
        P,    ///< p, bar p (proton, antiproton)
        PP,   ///< p+
        PM,   ///< p-
        NE,   ///< neutron, antineutron
        H,    ///< sum of charged hadrons
        HP,   ///< postiive hadron
        HM    ///< negative hadron
    };

    /**
     * Fourier transfer method used
     */
    enum FT_Method
    {
            GSL,            ///< Direct calculation using GSL integration method
            ACC_SERIES      ///< Series method implemented in fourier/fourier.c
    };

    /**
     * Virtual photon polarizatoin states
     */
    enum Polarization
    {
        L,  ///<< Longitudinal polarization
        T   ///<< Transverse polarization
    };

}



#endif
