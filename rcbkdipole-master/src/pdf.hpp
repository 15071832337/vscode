#ifndef _PDF_H
#define _PDF_H

/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "config.hpp"
#include <string>

/**
 * Parton distribution function
 *
 * Virtual class from which all PDF classes are inherited.
 */
class PDF
{
    public:
		PDF();
        ~PDF();
        /**
         * Evaluate x*pdf
         *
         * Return x*f(x,Q^2) for the given parton. Every PDF class
         * must implement this.
         * @param x Bjorken-x
         * @param q Scale [GeV]
         * @param p Parton type
         */
        virtual double xq(double x, double q, Amplitude::Parton p)=0;    // return x*q(x,q), q in GeV
        // x_1*x_2*f(x_1,x_2) dpfd with kinematical constraint x_1+x_2 < 1
        double Dpdf(double x1, double x2, double q, Amplitude::Parton p1, Amplitude::Parton p2);	
        virtual void Initialize(int param=-1);
        virtual std::string GetString();
        virtual double MinX();
        virtual double MinQ();   // GeV
        virtual double MaxX();
        virtual double MaxQ();
        virtual void SetOrder(Amplitude::Order o);
        virtual int SetA(int A_);	// set A for nuclear pdf, returns -1 if does not support other than proton
        
        virtual void Test();

        void PlotPdf(double Q);
	protected:
		bool initialized;
		Amplitude::Order order;
		int A;
};


#endif
