/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi> 2011-2014
 */

#include "pdf.hpp"

using namespace std;
using namespace Amplitude;

PDF::PDF()
{
	initialized=false;
	A=1;
}

void PDF::Initialize(int param)
{
}

std::string PDF::GetString()
{
    return "not specified";
}

PDF::~PDF()
{
}

void PDF::PlotPdf(double Q)
{
    cout << "# PDF at Q=" << Q << "GeV" << endl;
    cout <<"# x     up    down   s  gluon" << endl;
    for (double x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << xq(x,Q,U) << " " <<  xq(x,Q,D) << " "
        << xq(x,Q,S) << " " << xq(x,Q,G) << endl;
    }

}

/*
 * Double parton distribution function
 * Following e.g. ref 1009.6123
 * Returns  x_1x_2f(x_1,x_2)  with kinematical constraint x_1+x_2 < 1
 */

double PDF::Dpdf(double x1, double x2, double scale, Parton p1, Parton p2)
{
	if (x1 + x2 >= 1.0) return 0;
	return 0.5 * (
				xq(x1, scale , p1)/x1
				 * xq( x2/(1.0-x1), scale, p2 ) / ( x2/(1.0-x1) ) 
			   + xq(x2, scale , p2)/x2
				* xq( x1/(1.0-x2), scale, p1 ) / ( x1/(1.0-x2) )
			) * x1*x2;
}	

double PDF::MaxX()
{
    cerr << "PDF::MaxX() is not implemented! " << LINEINFO << endl;
    return 1;
}
double PDF::MaxQ()
{
    cerr << "PDF::MaxQ() is not implemented! " << LINEINFO << endl;
    return 1e9;
}

double PDF::MinX()
{
    cerr << "PDF::MinX() is not implemented! " << LINEINFO << endl;
    return 0;
}
double PDF::MinQ()
{
    cerr << "PDF::MinQ() is not implemented! " << LINEINFO << endl;
    return 0;
}

void PDF::SetOrder(Order o)
{
	cerr << "PDF::SetOrder is not implemented! " << LINEINFO << endl;
}

void PDF::Test()
{
	cerr << "PDF::Test is not implemented! " << LINEINFO << endl;
}

int PDF::SetA(int A_)
{
	cerr << "PDF::SetA() is not implemented, most likely the current PDF (" << GetString() <<") only supports proton (A=1)" << endl;
	return -1;
}
