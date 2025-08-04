/*
 * KKP Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include "kkp.hpp"
#include "config.hpp"
#include "fragmentation.hpp"

extern "C"
{
    // parton labels in output array:
    //  0    1    2    3    4    5    6    7    8     9    10
    //  g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
    
    void kkp_(int& h, int& set, double& x, double &qs, double output[]);
}

using namespace Amplitude;

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
double KKP::Evaluate(Parton p, Hadron h, double x, double qs)
{
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }
    int hadron=0; 
    double coeff=1.0;
    double partons[11];
    partons[0]=0; partons[1]=1; partons[2]=2; partons[3]=3; partons[10]=10;
    if (h==PI) { coeff=2.0; hadron=1; }
    else if (h==K) hadron=2;
    else if (h==K0) hadron=3;
    else if (h==P) hadron=4;
    else if (h==PI0) hadron=5;
    else if (h==NE) hadron=6;
    else if (h==H) hadron=7;
    else
    {
        cerr <<"Hadron type " << h << " (" << ParticleStr(h) << ") is not supported by KKP fragmentation! "
            << LINEINFO << endl;
        return 0;
    }

	int set=0;
	if (order == NLO)
		set=1; 
    kkp_(hadron, set, x, qs, partons);
    double result=0;
    if (p==G) result = partons[0];
    else if (p==U) result = partons[1];
    else if (p==UBAR) result  =partons[2];
    else if (p==D) result = partons[3];
    else if (p==DBAR) result = partons[4];
    else if (p==S) result = partons[5];
    else if (p==SBAR) result = partons[6];
    else if (p==C) result = partons[7];
    else if (p==CBAR) result = partons[8];
    else if (p==B) result = partons[9];
    else if (p==BBAR) result = partons[10];
    else
        cerr << "Parton " << p << " is not supported! " << LINEINFO << endl;
    
    return result*coeff;

}



std::string KKP::GetString()
{
    return "KKP";
}

KKP::KKP()
{

}

/*
 * Test
 */
void KKP::Test()
{
	cout <<"#------- Testing KKP FF" << endl;
	cout << "# NLO: " << endl;
	Order orig_order = order;
	
	SetOrder(NLO);
	double result, cor;
	result=0.1*Evaluate(U, PI0, 0.1, std::sqrt(20)); cor=0.82866;
	cout <<"z*D_{u->pi0}(Q^2=20GeV^2, x=0.1) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
		
	result=0.2*Evaluate(D, PI, 0.2, std::sqrt(10)); cor=0.85339;
	cout <<"z*D_{d->pi^+ + pi^-}(Q^2=10GeV^2, x=0.2) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
		
	result=0.5*Evaluate(G, PI0, 0.5, std::sqrt(100)); cor=0.082828;
	cout <<"z*D_{g->pi0}(Q^2=100GeV^2, x=0.5) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	
	SetOrder(orig_order);
	
}
 
