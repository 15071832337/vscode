//Estimate dihadron correlation cross section numeric multiDim integral in ROOT 
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

string datafile="/Users/suyoupeng/Downloads/rcbkdipole-master/data/proton/mvgamma.dat";
// Read data
AmplitudeLib N(datafile);

double Scorre(double xg, double r, int rep){

  double Nval=N.N(r,xg);

  if(rep>0) return 1-Nval; //fundamental rep
  else if(rep<0) return 1-2*Nval+Nval*Nval; //adjoint rep
  else
    cout<<"sth wrong!"<<endl;

  return 0;
}

// The function is used for performing Fourier transforms.
double functor(double *x, double *par){
  /*
  par[0] xg : Bjorken x
  par[1] kt : transverse momentum
  par[2] rep:different rep, 
    1 fundamental rep
    -1 adjoint rep
  */
  double xg = par[0];
  double kt = par[1];
  double rep = par[2];
  double r = x[0];

  return Scorre(xg, r, rep)*TMath::BesselJ0(kt*r)*r*2*TMath::Pi();
}

// This function  is used to calculate the S(k) after performing the Fourier transform.
int main()
{
  string outputFileName = "S_k_in_fund_mv_1e-2.txt";
  ofstream outputFile(outputFileName);

  double xbj=0.01;
  //double xbj=0.00001;

  N.InitializeInterpolation(xbj);
  cout << "r [1/GeV]   N(r,x=" << xbj << ")" << endl;
  const double Ns =  1.0-std::exp(-0.5);
  cout << " Q_s^2(x=0.01) = " << N.SaturationScale(0.01, Ns) << " GeV^2" <<endl;
  cout << "Q_s^2(x=0.0001) = " << N.SaturationScale(0.0001, Ns) << " GeV^2" << endl;

  TF1 *fun = new TF1("fun", functor, 1e-1, 50, 3);
  fun->SetParameters(0.01, 2.0, 1);
  
  fun->Draw();

  TGraph* graph = new TGraph();
  graph->SetNameTitle("S(k)", "; k [GeV]; S(k)");
   
  double minr = N.MinR(); double maxr=N.MaxR();

  int i = 1;
  for(double k=0.01; k<maxr; k=k+0.01){
    fun->SetParameters(0.01, k, 1);
    double S_k = fun -> Integral(1e-6, 100);
    graph->SetPoint(i, k, S_k);
    outputFile << std::scientific << std::setprecision(9) << k << " " << S_k  << endl;
    i++;
    cout<<k<<endl;
  }

  graph->GetXaxis()->SetRangeUser(1e-1, 15); 
  graph->GetYaxis()->SetRangeUser(1e-4, 1e+2);
  graph->Draw("APL");

  outputFile.close();
  return 0;
}

// This function is used to plot the distribution of S(k) with respect to  k .
void draw(){
  //fund rep = 1
  string datafile1 = "S_k_in_fund_mv_1e-2.txt"; 
  string datafile2 = "S_k_in_fund_mve_1e-2.txt";
  string datafile3 = "S_k_in_fund_mvgamma_1e-2.txt"; 

  //adjoint rep = -1
  string datafile4 = "S_k_in_adj_mv_1e-2.txt"; 
  string datafile5 = "S_k_in_adj_mve_1e-2.txt";
  string datafile6 = "S_k_in_adj_mvgamma_1e-2.txt"; 

  TGraph* graph1 = new TGraph();
  graph1->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph2 = new TGraph();
  graph2->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph3 = new TGraph();
  graph3->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");


  TGraph* graph4 = new TGraph();
  graph4->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph5 = new TGraph();
  graph5->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph6 = new TGraph();
  graph6->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");


  ifstream inputFile1(datafile1);
  ifstream inputFile2(datafile2);
  ifstream inputFile3(datafile3);

  ifstream inputFile4(datafile4);
  ifstream inputFile5(datafile5);
  ifstream inputFile6(datafile6);

  double k, Sp;
  int i = 1;
  while (inputFile1 >> k >> Sp)
  {
      graph1->SetPoint(i, k, Sp);
      i++;
  }
  inputFile1.close();


  while (inputFile2 >> k >> Sp)
  {
      graph2->SetPoint(i, k, Sp);
      i++;
  }
  inputFile2.close();

  while (inputFile3 >> k >> Sp)
  {
      graph3->SetPoint(i, k, Sp);
      i++;
  }
  inputFile3.close();

  while (inputFile4 >> k >> Sp)
  {
      graph4->SetPoint(i, k, Sp);
      i++;
  }
  inputFile4.close();


  while (inputFile5 >> k >> Sp)
  {
      graph5->SetPoint(i, k, Sp);
      i++;
  }
  inputFile5.close();



  while (inputFile6 >> k >> Sp)
  {
      graph6->SetPoint(i, k, Sp);
      i++;
  }
  inputFile6.close();

  graph1->SetLineColor(kRed); 
  graph1->SetLineStyle(2);      
  graph1->SetLineWidth(6);     


  graph2->SetLineColor(kBlack); 
  graph2->SetLineStyle(2);      
  graph2->SetLineWidth(6);     

  graph3->SetLineColor(kBlue); 
  graph3->SetLineStyle(1);      
  graph3->SetLineWidth(6);  
  
  graph4->SetLineColor(kRed); 
  graph4->SetLineStyle(2);      
  graph4->SetLineWidth(3);     

  graph5->SetLineColor(kBlack); 
  graph5->SetLineStyle(2);      
  graph5->SetLineWidth(3);     

  graph6->SetLineColor(kBlue); 
  graph6->SetLineStyle(1);      
  graph6->SetLineWidth(3);    

  graph1->GetXaxis()->SetRangeUser(1e-1, 15); 
  graph1->GetYaxis()->SetRangeUser(1e-4, 20+1e+2);  

  TCanvas* canvas = new TCanvas("canvas", "Dipole amplitude", 800, 600);
  graph1->Draw("AL");
  graph2->Draw("SAME");
  graph3->Draw("SAME");
  graph4->Draw("SAME");
  graph5->Draw("SAME");
  graph6->Draw("SAME");

  TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
  legend->AddEntry(graph1, "Model:MV,x = 1e-2,fundamental rep", "l");   
  legend->AddEntry(graph2, "Model:MV^{e},x = 1e-2,fundamental rep", "l");   
  legend->AddEntry(graph3, "Model:MV^{#gamma},x = 1e-2,fundamental rep", "l");   
  legend->AddEntry(graph4, "Model:MV,x = 1e-2,adjoint rep", "l");   
  legend->AddEntry(graph5, "Model:MV^{e},x = 1e-2,adjoint rep", "l");   
  legend->AddEntry(graph6, "Model:MV^{#gamma},x = 1e-2,adjoint rep", "l");  
  legend->SetBorderSize(0);                           
  legend->SetFillColor(0);                             
  legend->Draw();  
  
}


// This is the function of α_s with respect to Q^{2}, which has been extracted by fitting the data of α_s as it varies with Q^{2} 
// The fitting parameters are already provided within it. 
// beta0 = 9.2791 and lambdaQCD = 0.218344
// Q: momentum transfer
double alphas(Double_t Q) {
  // p0 = 9.2791
  // p1 = 0.218344
  double mu = Q;
  double beta0 = 9.2791;
  double lambdaQCD = 0.218344;
  return 2 * TMath::Pi() / (beta0 * TMath::Log(mu / lambdaQCD));
}

double Alphas(double Q){
  double lqcd = 0.241;
  double alpha = 12.0*TMath::Pi()/( (11.0*3-2.0*3)*log(Q/(lqcd*lqcd)) );
  return alpha;
}


// This function is used to calculate Equation 9 from the literature.
// k : transverse momentum
// x : Bjorken x
double Phip_k(double k,double x){

  double xbj=x;
  //double xbj=0.00001;
  N.InitializeInterpolation(xbj);
  TF1 *fun = new TF1("fun", functor, 1e-1, 50, 3);
  //fun->SetParameters(0.01, k, 1);
  //double S_k = fun -> Integral(1e-6, 100); 

  //sigma0/2  mb：
  // MV 18.81 
  // MV gamma 16.45
  // MV e 16.36
  //alpha_s 0.7
  //Cf = 4/3 fundamental rep 
  double sigma0Over2 = 18.81; //mb
  double CF = 4.0 / 3.0;
  double alpha_s = 0.1181;
  
  fun->SetParameters(xbj, k, -1);
  double S_k = fun -> Integral(1e-6, 100);
  return (CF*sigma0Over2/(8*pow(TMath::Pi(),3)*alpha_s))*pow(k,4) * S_k ;
  //return (CF*sigma0Over2/(8*pow(TMath::Pi(),3)))*pow(k,4) * S_k ;
}

//This function is used to calculate the integrand in Equation 13 from the literature.
//par[0] x_g : Bjorken x
// qt : transverse momentum
double PhipOverqTSquare(double *x,double *par){
  double qt = x[0];
  double x_g = par[0];
  return 2 * Phip_k(qt,x_g) / (qt);
}

// This function represents Equation 13.
// kt : transverse momentum
// x : Bjorken x
double xg(double kt,double x){
  TF1 *fun = new TF1("fun", PhipOverqTSquare, 1, 1e+3,1);
  fun->SetParameters(x);
  fun->Draw();
  double xg_kt = fun -> Integral(1e-8, kt);
  return xg_kt;
}

// This function is used to calculate the gluon distribution function depicted in Figure 4 of the article.
void draw_xg_kt(){

  string outputFileName = "xg1e-2_MV.txt";
  ofstream outputFile(outputFileName);


  TGraph* graph = new TGraph();
  graph->SetNameTitle("Gluon distribution function at x = 10^{-2}", "; k [GeV]; S(k)");

  int k =1;
  for(double i = 0;i<= 20;i=i+0.5){
    outputFile << std::scientific << std::setprecision(9) << i*i << " " << xg(i*i,0.01) << endl;
    //outputFile << std::scientific << std::setprecision(9) << i*i << " " << xg(i*i,0.01) /alphas(i*i) << endl;
    //graph->SetPoint(k, i*i,xg(i*i));
    k++;
    cout<<"-------------"<<i<<"-------------"<<endl;
  }
  graph->GetXaxis()->SetRangeUser(0, 1e+3); 
  graph->GetYaxis()->SetRangeUser(1e+0, 1e+2);

  TH1D *h = new TH1D("h","Dipole amplitude",1,1,1000);
  h->GetYaxis()->SetRangeUser(1,1e+2);
  h->SetStats(0);
  h->Draw();
  graph->Draw("SAME");

  outputFile.close();
}

// This represents the relationship between Bjorken x , transverse momentum pt, and rapidity y.
// sqrt_s : collision energy
// k : transverse momentum
double x_g(double k){
  double sqrt_s = 200;
  double y = 3.3;
  return (k/sqrt_s) * TMath::Exp(-y);
}


// This function is used to calculate the invariant yields for \( pp \) collisions as depicted in Figure 6.
void dNoverdyd2kt(){ //we use σinel = 42 mb at RHIC [39], σinel = 44.4 mb at Tevatron [40] and σinel = 70 mb at LHC energies [41]).
  string outputFileName = "MVgamma_dNOverdydpT_y=3_3.txt";
  ofstream outputFile(outputFileName);
  double sigma0Over2 = 16.45;
  double sigma_inel = 42;
  double A = pow(2*TMath::Pi(),2);
  
  TF1 *fun = new TF1("fun", functor, 1e-1, 50, 3);
  
  for(double i = 1;i<= 10;i=i+0.1){
    //outputFile << std::scientific << std::setprecision(9) << i*i << " " << xg(i*i) << endl;
    double x = x_g(i);
    fun->SetParameters(x, i, -1);

    double S_k = fun -> Integral(1e-6, 100);
    //double result = (sigma0Over2/sigma_inel) * 1/A * (xg(i*i,x) / alphas(i*i)) * S_k;
    double result = (sigma0Over2/sigma_inel) * 1/A * (xg(i*i,x)) * S_k;
    outputFile << std::scientific << std::setprecision(9) << i << " " << result << endl;
  
    cout<<"-------------"<<i<<"-------------"<<endl;
    cout<<"-------------"<<result<<"-------------"<<endl;

  }

  outputFile.close();
}

void drawdNOverdyd2kt(){
  string datafile1 = "dNOverdydpT_y=3_8.txt"; 
  string datafile2 = "dNOverdydpT_y=3_3.txt"; //MVgamma_dNOverdydpT_y=3_3
 

  string datafile3 = "MVe_dNOverdydpT_y=3_3.txt"; 
  string datafile4 = "MVe_dNOverdydpT_y=3_8.txt"; 

  string datafile5 = "MV_dNOverdydpT_y=3_3.txt"; 
  string datafile6 = "MV_dNOverdydpT_y=3_8.txt"; 

  string datafile7 = "data_MVgamma_dNOverdydpT_y=3_3.txt"; 
  
  ifstream inputFile1(datafile1);
  ifstream inputFile2(datafile2);

  ifstream inputFile3(datafile3);
  ifstream inputFile4(datafile4);

  ifstream inputFile5(datafile5);
  ifstream inputFile6(datafile6);

  ifstream inputFile7(datafile7);
  
  TGraph* graph1 = new TGraph();
  graph1->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph2 = new TGraph();
  graph2->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph3 = new TGraph();
  graph3->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph4 = new TGraph();
  graph4->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph5 = new TGraph();
  graph3->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph6 = new TGraph();
  graph4->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");

  TGraph* graph7 = new TGraph();
  graph7->SetNameTitle("dN/dyd^{2}k_{T}", "dN/dyd^{2}k_{T}; p_{T} [GeV]; dN/dyd^{2}k_{T}");


  double k, Sp;
  int i = 1;
  while (inputFile1 >> k >> Sp)
  {
      graph1->SetPoint(i, k, Sp);
      i++;
  }

  while (inputFile2 >> k >> Sp)
  {
      graph2->SetPoint(i, k, Sp);
      i++;
  }
  while (inputFile3 >> k >> Sp)
  {
      graph3->SetPoint(i, k, Sp);
      i++;
  }
  while (inputFile4 >> k >> Sp)
  {
      graph4->SetPoint(i, k, Sp);
      i++;
  }
  while (inputFile5 >> k >> Sp)
  {
      graph5->SetPoint(i, k, Sp);
      i++;
  }
  while (inputFile6 >> k >> Sp)
  {
      graph6->SetPoint(i, k, Sp);
      i++;
  }

  while (inputFile7 >> k >> Sp)
  {
      graph7->SetPoint(i, k, Sp);
      i++;
  }


  TH1D *h = new TH1D("h","dN/dyd^{2}p_{T}; p_{T} [GeV/c]; dN/dyd^{2}p_{T} [1/GeV^{2}]",1,1,10);
  h->GetYaxis()->SetRangeUser(1e-8,1);
  h->SetStats(0);
  h->Draw();

  
  graph1->SetLineColorAlpha(kYellow,1); 
  graph1->SetLineStyle(2);      
  graph1->SetLineWidth(6);   
  graph2->SetLineColorAlpha(kYellow,1); 
  graph2->SetLineStyle(1);      
  graph2->SetLineWidth(6);   

  graph3->SetLineColorAlpha(kBlue,1); 
  graph3->SetLineStyle(1);      
  graph3->SetLineWidth(6);   
  graph4->SetLineColorAlpha(kBlue,1); 
  graph4->SetLineStyle(2);      
  graph4->SetLineWidth(6);  

  graph5->SetLineColorAlpha(kBlack,1); 
  graph5->SetLineStyle(1);      
  graph5->SetLineWidth(6);   
  graph6->SetLineColorAlpha(kBlack,1); 
  graph6->SetLineStyle(2);      
  graph6->SetLineWidth(6);  

  graph7->SetLineColorAlpha(kRed,1); 
  graph7->SetLineStyle(2);      
  graph7->SetLineWidth(6); 



  graph1->Draw("SAME");
  graph2->Draw("SAME");

  graph3->Draw("SAME");
  graph4->Draw("SAME");

  graph5->Draw("SAME");
  graph6->Draw("SAME");

  graph7->Draw("SAME");

  TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
  legend->AddEntry(graph1, "MV^{#gamma},y = 3.8 ", "lp"); 
  legend->AddEntry(graph2, "MV^{#gamma},y = 3.3 ", "lp");
  legend->AddEntry(graph4, "MV^{e},y = 3.8 ", "lp");
  legend->AddEntry(graph3, "MV^{e},y = 3.3 ", "lp"); 
  legend->AddEntry(graph6, "MV,y = 3.8 ", "lp");
  legend->AddEntry(graph5, "MV,y = 3.3 ", "lp");
  legend->AddEntry(graph7, "Fig6:MV^{#gamma},y = 3.3 ", "lp");


  legend->Draw("SAME");

}


Double_t alphas(Double_t *x, Double_t *params) {
  // p0 = 9.2791
  // p1 = 0.218344
  Double_t mu = x[0];
  Double_t beta0 = params[0];
  Double_t lambdaQCD = params[1];
  return 2 * TMath::Pi() / (beta0 * TMath::Log(mu / lambdaQCD));
}

//This function is used to fit the alpha_s  data to extract the parameters beta and Lambda_{QCD}
void fitAlphas() {

  string datafile1 = "alpha_s.txt"; 
  ifstream inputFile1(datafile1);
  TGraph* graph1 = new TGraph();
  graph1->SetNameTitle("#alpha_{s}(Q)", "#alpha_{s}(Q); Q [GeV]; #alpha_{s}(Q)");

  TGraph* graph2 = new TGraph();
  graph2->SetNameTitle("#alpha_{s}(Q)", "#alpha_{s}(Q); Q [GeV]; #alpha_{s}(Q)");

  TH1D *h = new TH1D("#alpha_{s}(Q)", "#alpha_{s}(Q); Q [GeV]; #alpha_{s}(Q)",1,1,1000);
  h->GetYaxis()->SetRangeUser(0.05,0.4);
  h->SetStats(0);
  h->Draw();

  double q, alpha_s;
  int i = 1;
  while (inputFile1 >> q >> alpha_s)
  {
      graph1->SetPoint(i, q, alpha_s);
      i++;
  }
  for(int k = 1;k<1000;k++){
    graph2->SetPoint(k, k, Alphas(k));
  }

  graph2->SetLineColorAlpha(kBlack,1); 
  graph2->SetLineStyle(1);      
  graph2->SetLineWidth(6);  


  graph1->Draw("SAME");

  TF1 *func = new TF1("alphas", alphas, 1, 1000, 2); 
  func->SetParameters(1, 0.2); 

  func->SetNpx(100);
  func->SetLineWidth(2);
  func->SetLineColor(kRed);

  
  double x1 = 0;
  double x2 = 1000;
  double y = 0.1181; 

  TLine *line = new TLine(x1, y, x2, y);
  line->SetLineColor(kBlue);  
  line->SetLineStyle(1);  
  line->SetLineWidth(2);  

  graph1->Fit("alphas", "R");
  func->Draw("same");
  graph2->Draw("same");
  line->Draw("same");

  TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
  legend->AddEntry(graph1, "data ", "l"); 
  legend->AddEntry(func, "#alpha_{s}(Q^{2})", "l"); 
  legend->AddEntry(line, "#alpha_{s} = 0.1181", "l");   
  legend->AddEntry(graph2, "#alpha_{s} from (7)", "l"); 
  legend->Draw();
 
}



// This function is used to plot the dependence of the calculated invariant yields on transverse momentum.
void draw2(){
  string datafile1 = "xg1e-5_MV.txt"; 
  string datafile4 = "xg1e-5_MV_corrected.txt"; 
  string datafile2 = "xg1e-5_MVgamma.txt"; 
  string datafile3 = "data_MVgamma_x=1e-4.txt"; 
  string datafile5 = "xg1e-5_MVgamma_correted.txt"; 
  string datafile6 = "data_MV_x=1e-4.txt"; 
  string datafile7 = "data_MV_x=1e-2.txt"; 
  string datafile8 = "xg1e-2_MV.txt"; 
  string datafile9 = "xg1e-2_MV_corrected.txt"; 

  string datafile10 = "data_MVgamma_x=1e-2.txt"; 
  string datafile11 = "xg1e-2_MVgamma.txt"; 

  string datafile12 = "xg1e-2_MVgamma_corrected.txt"; 


  ifstream inputFile1(datafile1);
  ifstream inputFile2(datafile2);
  ifstream inputFile3(datafile3);
  ifstream inputFile4(datafile4);
  ifstream inputFile5(datafile5);
  ifstream inputFile6(datafile6);
  ifstream inputFile7(datafile7);
  ifstream inputFile8(datafile8);
  ifstream inputFile9(datafile9);

  ifstream inputFile10(datafile10);
  ifstream inputFile11(datafile11);

  ifstream inputFile12(datafile12);
  

  TGraph* graph1 = new TGraph();
  graph1->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph2 = new TGraph();
  graph2->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph3 = new TGraph();
  graph3->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph4 = new TGraph();
  graph4->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph5 = new TGraph();
  graph5->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph6 = new TGraph();
  graph6->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph7 = new TGraph();
  graph7->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph8 = new TGraph();
  graph8->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph9 = new TGraph();
  graph9->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph10 = new TGraph();
  graph10->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  TGraph* graph11 = new TGraph();
  graph11->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");

  TGraph* graph12 = new TGraph();
  graph12->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  
  double k, Sp;
  int i = 1;
  while (inputFile1 >> k >> Sp)
  {
      graph1->SetPoint(i, k, Sp*0.258);
      i++;
  }
  while (inputFile2 >> k >> Sp)
  {
      graph2->SetPoint(i, k, Sp*0.258);
      i++;
  }
  while (inputFile3 >> k >> Sp)
  {
      graph3->SetPoint(i, k, Sp);
      i++;
  }
  while (inputFile4 >> k >> Sp)
  {
      graph4->SetPoint(i, k, Sp*0.258);
      i++;
  }
  while (inputFile5 >> k >> Sp)
  {
      graph5->SetPoint(i, k, Sp*0.258);
      i++;
  }
  
  while (inputFile6 >> k >> Sp)
  {
      graph6->SetPoint(i, k, Sp);
      i++;
  }

  while (inputFile7 >> k >> Sp)
  {
      graph7->SetPoint(i, k, Sp);
      i++;
  }

  while (inputFile8 >> k >> Sp)
  {
      graph8->SetPoint(i, k, Sp*0.258);
      i++;
  }
  while (inputFile9 >> k >> Sp)
  {
      graph9->SetPoint(i, k, Sp*0.258);
      i++;
  }

  while (inputFile10 >> k >> Sp)
  {
      graph10->SetPoint(i, k, Sp);
      i++;
  }

  while (inputFile11 >> k >> Sp)
  {
      graph11->SetPoint(i, k, Sp*0.258);
      i++;
  }

  while (inputFile12 >> k >> Sp)
  {
      graph12->SetPoint(i, k, Sp*0.258);
      i++;
  }



  graph1->SetLineColorAlpha(kRed,0.35); 
  graph1->SetLineStyle(2);      
  graph1->SetLineWidth(6);   

  graph4->SetLineColorAlpha(kRed,0.75); 
  graph4->SetLineStyle(2);      
  graph4->SetLineWidth(6);  

  graph6->SetLineColorAlpha(kRed,1); 
  graph6->SetLineStyle(1);      
  graph6->SetLineWidth(6);

  graph2->SetLineColorAlpha(kBlack,0.35); 
  graph2->SetLineStyle(2);      
  graph2->SetLineWidth(6);  

  graph3->SetLineColorAlpha(kBlack,1); 
  graph3->SetLineStyle(1);      
  graph3->SetLineWidth(6);   

  graph5->SetLineColorAlpha(kBlack,0.75); 
  graph5->SetLineStyle(2);      
  graph5->SetLineWidth(6);


  graph7->SetLineColorAlpha(kBlue, 1); 
  graph7->SetLineStyle(1);      
  graph7->SetLineWidth(6);

  graph8->SetLineColorAlpha(kBlue, 0.35); 
  graph8->SetLineStyle(2);      
  graph8->SetLineWidth(6);

  graph9->SetLineColorAlpha(kBlue, 0.75); 
  graph9->SetLineStyle(2);      
  graph9->SetLineWidth(6);

  graph10->SetLineColorAlpha(kGreen, 1); 
  graph10->SetLineStyle(1);      
  graph10->SetLineWidth(6);

  graph11->SetLineColorAlpha(kGreen, 0.35); 
  graph11->SetLineStyle(2);      
  graph11->SetLineWidth(6);

  graph12->SetLineColorAlpha(kGreen, 0.75); 
  graph12->SetLineStyle(2);      
  graph12->SetLineWidth(6);

  

  TH1D *h = new TH1D("h","Gluon distribution function; Q^{2} [GeV^{2}]; xg(x,Q^{2})",1,1,1000);
  h->GetYaxis()->SetRangeUser(1,1e+2);
  h->SetStats(0);
  h->Draw();

  
  graph1->Draw("SAME");
  graph2->Draw("SAME");
  graph3->Draw("SAME");
  graph4->Draw("SAME");
  graph5->Draw("SAME");
  graph6->Draw("SAME");
  graph7->Draw("SAME");
  graph8->Draw("SAME");
  graph9->Draw("SAME");
  graph10->Draw("SAME");
  graph11->Draw("SAME");
  graph12->Draw("SAME");
  

  TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
  legend->AddEntry(graph1, "x 0.075 MV,x=1e-5", "l");   
  legend->AddEntry(graph4, "x 0.16 corrected:MV,x=1e-5", "l"); 
  legend->AddEntry(graph6, "data(from fig 4):MV,x=1e-5", "l"); 

  legend->AddEntry(graph2, "x 0.075 MV^{#gamma},x=1e-5", "l");   
  legend->AddEntry(graph5, "x 0.16 corrected:MV^{#gamma},x=1e-5 ", "l"); 
  legend->AddEntry(graph3, "data(from fig 4):MV^{#gamma},x=1e-5", "l"); 

  legend->AddEntry(graph8, "x 0.15 MV,x=1e-2", "l");  
  legend->AddEntry(graph9, "x 0.2158 corrected:MV,x=1e-2", "l");  
  legend->AddEntry(graph7, "data(from fig 4):MV,x=1e-2", "l"); 
  legend->AddEntry(graph10, "data(from fig 4):MV^{#gamma},x=1e-2", "l"); 
  legend->AddEntry(graph11, "x 0.21 MV^{#gamma},x=1e-2", "l"); 
  legend->AddEntry(graph12, "x 0.32 corrected:MV^{#gamma},x=1e-2", "l"); 

  legend->SetBorderSize(0);                           
  legend->SetFillColor(0);                             
  legend->Draw();  

}

void test(){
  string datafile9 = "xg1e-2_MV_corrected.txt"; 
  ifstream inputFile9(datafile9);
  TGraph* graph9 = new TGraph();
  graph9->SetNameTitle("#tilde{S}(k)", "Two dimensional Fourier transform of S(r); k [GeV]; #tilde{S}(k)");
  double k, Sp;
  int i = 1;
  while (inputFile9 >> k >> Sp)
  {
      graph9->SetPoint(i, k, Sp*0.0258);
      i++;
  }
  TH1D *h = new TH1D("h","Gluon distribution function; Q^{2} [GeV^{2}]; xg(x,Q^{2})",1,1,1000);
  h->GetYaxis()->SetRangeUser(1,1e+2);
  h->SetStats(0);
  //h->Draw();
  graph9->Draw();


}