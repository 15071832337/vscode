//Estimate dihadron correlation cross section numeric multiDim integral
//in ROOT 
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <ctime>
#include <sstream>
#include <unistd.h>

#include "/Users/suyoupeng/Downloads/rcbkdipole-master/src/amplitudelib.hpp"

#include "TF1.h"
#include "TF2.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"

#include "Math/SpecFuncMathCore.h" 
#include "Math/SpecFuncMathMore.h" 

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"

//#include "dssWrapper.C"

using namespace std;

const double PI=3.1415926;

string datafile="/Users/suyoupeng/Downloads/rcbkdipole-master/data/proton/mve.dat";
// Read data
AmplitudeLib N(datafile);

double Scorre(double xg, double r, int rep){

  double Nval=N.N(r,xg);

  if(rep>0) return 1-Nval;
  else if(rep<0) return 1-2*Nval+Nval*Nval;
  else
    cout<<"sth wrong!"<<endl;

  return 0;
}

double functor(double *x, double *par){

  double xg = par[0];
  double kt = par[1];
  double rep = par[2];
  double r = x[0];

  return Scorre(xg, r, rep)*TMath::BesselJ0(kt*r)*r*2*PI;
}



void main()
{

  double xbj=0.01;
  N.InitializeInterpolation(xbj);
  cout << "r [1/GeV]   N(r,x=" << xbj << ")" << endl;
  const double Ns =  1.0-std::exp(-0.5);
  cout << " Q_s^2(x=0.01) = " << N.SaturationScale(0.01, Ns) << " GeV^2" <<endl;
  cout << "Q_s^2(x=0.0001) = " << N.SaturationScale(0.0001, Ns) << " GeV^2" << endl;

  TF1 *fun = new TF1("fun", functor, 1e-1, 50, 3);
  fun->SetParameters(0.01, 2.0, -1);
  cout<<fun->Integral(1e-6, 100)<<endl;
  

  TCanvas *canvas = new TCanvas("canvas", "Function Plot", 800, 600);
  fun->SetLineColor(kBlue); // 设置函数颜色为蓝色
  fun->SetLineWidth(4);     // 设置函数线条宽度
  fun->Draw(); 
  
}
