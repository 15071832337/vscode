/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "qcd.hpp"
#include "/Users/suyoupeng/Downloads/rcbkdipole-master/src/config.hpp"

using namespace Amplitude;


QCD::QCD()
{
    // Initialize default values
    nc=3;
    nf=3;
    rc=RUNNING;
    lqcd = 0.241;
    maxalpha=0.7;
}

/*
 * Strong coupling
 */
double QCD::Alphas(double Qsqr)
{
	if (rc == FIXED)
        return 0.2;
	
	if (Qsqr <= lqcd*lqcd)
		return maxalpha;

    double alpha = 12.0*M_PI/( (11.0*Nc()-2.0*Nf())*log(Qsqr/(lqcd*lqcd)) );
    if (alpha > maxalpha)
        return maxalpha;

    return alpha;
}

double QCD::Cf()
{
    return (Nc()*Nc()-1.0)/(2.0*Nc());
}

int QCD::Nc()
{
    return nc;
}

int QCD::Nf()
{
    return nf;
}

void QCD::SetRunningCoupling(Amplitude::RunningAlphas rc_)
{
    rc=rc_;
}

