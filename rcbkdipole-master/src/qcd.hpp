#ifndef _QCD_HPP
#define _QCD_HPP

#include "config.hpp"

/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
 */

/**
 * QCD functions and constants
 *
 * Strong coupling constant, Nc, Nf etc. parameters
 */

 class QCD
 {
     public:
         QCD();
         
         /**
          * Strong coupling constant
          * @param musqr scale [GeV^2]
          */
         double Alphas(double musqr);
         /**
          * Number of quark flavors
          */
         int Nf();
         /**
          * Number of colors
          */
         int Nc();
 
         /**
          * Fundametal Casimir
          */
         double Cf();
 
         /**
          * Set running/fixed coupling
          */
         void SetRunningCoupling(Amplitude::RunningAlphas rc_);
     private:
         int nc;
         int nf;
         double lqcd;    // lambda_QCD
         double maxalpha;    // freeze coupling to maxalpha in infrared
         Amplitude::RunningAlphas rc;
         
 
 };


#endif

