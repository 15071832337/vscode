#include <gsl/gsl_errno.h>
#include <string>
#include <cassert>

#include "amplitudelib.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <TSystem.h>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_min.h>
#include <sstream>
#include "/Users/suyoupeng/Downloads/rcbkdipole-master/src/single_inclusive.hpp"


using namespace std;

string datafile="/Users/suyoupeng/Downloads/rcbkdipole-master/data/proton/mvgamma.dat";
string datafile_probe="/Users/suyoupeng/Downloads/rcbkdipole-master/data/proton/mvgamma.dat";

SingleInclusive::SingleInclusive(AmplitudeLib* N_)
{
    N=N_; 
}
// Read data
AmplitudeLib N(datafile);
AmplitudeLib N2(datafile_probe);
SingleInclusive xs(&N);


double dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N,AmplitudeLib* N2, double scale );


double Inthelperf_xg(double qsqr, void* p);

struct Inthelper_xg
{
	double xbj,q; AmplitudeLib* N;
};

double Scorre(double xg, double r, int rep){

    double Nval=N.N(r,xg);
  
    if(rep>0) return 1-Nval; //fundamental rep
    else if(rep<0) return 1-2*Nval+Nval*Nval; //adjoint rep
    else
      cout<<"sth wrong!"<<endl;
  
    return 0;
}

double functor(double *x, double *par){
    double xg = par[0];
    double kt = par[1];
    double rep = par[2];
    double r = x[0];
    return Scorre(xg, r, rep)*TMath::BesselJ0(kt*r)*r*2*TMath::Pi();
}

double S_k(double q,double xbj,int rep){
    TF1 *fun = new TF1("fun", functor, 1e-1, 100, 3);
    fun->SetParameters(xbj, q, rep);
    double S_k = fun -> Integral(1e-6, 100);
    return S_k;
}

double Nc(){
    return 3;
}
double Cf()
{
    return (Nc()*Nc()-1.0)/(2.0*Nc());
}
double Nf(){
    return 3;
}
double Alphas(double Qsqr)
{
    double lqcd = 0.241;
    double maxalpha = 0.7;
	
	if (Qsqr <= lqcd*lqcd)
		return maxalpha;

    double alpha = 12.0*M_PI/( (11.0*Nc()-2.0*Nf())*log(Qsqr/(lqcd*lqcd)) );
    if (alpha > maxalpha)
        return maxalpha;
    return alpha;
}
double SQR(double x){
    return x*x;
}
//====================Dipole_UGD================================//
double Dipole_UGD(double q, double xbj, double scale_, double S_T)
{
	double scale;
	if (scale_<0) scale=q*q; else scale=scale_;
	double alphas = Alphas(scale);
	return Cf() / (8.0 * M_PI*M_PI*M_PI) * S_T/alphas * std::pow(q,4) * S_k(q, xbj, -1);
}

double Inthelperf_xg(double qsqr, void* p)
{
	Inthelper_xg* par = (Inthelper_xg*) p;
	return 1.0/qsqr * Dipole_UGD(std::sqrt(qsqr), par->xbj, N.SQR(par->q),1);
    //return 2.0/(qsqr) * Dipole_UGD(qsqr, par->xbj, N.SQR(par->q),18.81);
}


double xg(double x, double q)
{
	Inthelper_xg par; 
	par.xbj=x; par.q=q;
	gsl_function fun; fun.function=Inthelperf_xg;
	fun.params=&par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(100);
	int status = gsl_integration_qag(&fun, 0, N.SQR(q), 0, 0.01,100, GSL_INTEG_GAUSS51, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		cerr << "UGD integral failed at "<< ", x=" << x
			<< ", k_T=" << q << " res " << result << " relerr " << std::abs(abserr/result)  << endl;
    }
	return result;
}

/*
 * k_T factorization result for the dN/(d^2kt dy)
 * \alpha_s/S_T 2/C_f 1/k_T^2 \int d^2 q UGD(q)/q^2 UQD(k-q)/|k-q|^2
 * UGD is AmplitudeLib::UGD()
 * 
 * If N2!=NULL, it is used to compute \phi_2
 */
struct Inthelper_ktfact{double x1, x2, pt, qt, scale; AmplitudeLib *N1, *N2; };
double Inthelperf_ktfact_q(double q, void* p);
double Inthelperf_ktfact_phi(double phi, void* p);
const int INTPOINTS_KTFACT = 4;	

double Inthelperf_ktfact_phi(double phi, void* p)
{
	Inthelper_ktfact* par = (Inthelper_ktfact*)p;
	double kt_m_qt = std::sqrt( SQR(par->pt) + SQR(par->qt) - 2.0*par->pt*par->qt*std::cos(phi));
	// In this limit we would have ugd(0)/0 -> 0
	if (kt_m_qt < 1e-5)
		return 0;	

	double ugd1,ugd2;
	if (par->N2==NULL)
	{
		par->N1->InitializeInterpolation(par->x1);

        //Dipole_UGD(std::sqrt(qsqr), par->xbj, N.SQR(par->q),1);

		ugd1 = Dipole_UGD(par->qt, par->x1, par->scale,1);
		par->N1->InitializeInterpolation(par->x2);
		ugd2 = Dipole_UGD(kt_m_qt, par->x2, par->scale,1);
	} 
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				ugd1 = Dipole_UGD(par->qt, par->x1, par->scale,1);
			}
			#pragma omp section
			{
				ugd2 = Dipole_UGD(kt_m_qt, par->x2, par->scale,1);
			}
		}
		
	}
	if (ugd1 < 0 or ugd2 < 0) return 0;
	return par->qt * ugd1/SQR(par->qt) * ugd2/SQR(kt_m_qt);
}

double Inthelperf_ktfact_q(double q, void *p)
{
	Inthelper_ktfact* par = (Inthelper_ktfact*)p;
	par->qt=q;
	
	gsl_function fun; fun.function=Inthelperf_ktfact_phi;
	fun.params=par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1);
	int status = gsl_integration_qag(&fun, 0, M_PI, 0, 0.05,
		1, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
    result *= 2.0;	// as we integrate [0,pi], not [0, 2pi]   
   /* if (status)
    {
		cerr << "kt-factorization phi integral failed at " << LINEINFO <<", q=" << q
			<< ", k_T=" << par->pt << " relerr " << std::abs(abserr/result) << endl;
    }*/
	return result;
}


double dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N,AmplitudeLib* N2, double scale )
{
	double x1 = pt*std::exp(-y)/sqrts;
	double x2 = pt*std::exp(y)/sqrts;
	double y1 = std::log(N->X0()/x1);
	double y2;
	if (N2==NULL)
		y2 = std::log(N->X0()/x2);
	else
		y2 = std::log(N2->X0()/x2);
	
	if (y1<0 or y2<0)
	{
		cerr << "Evolution variables y1=" << y1 <<", y2=" << y2 <<" is not possible to compute, too large x. pt=" << pt << ", y=" << y << ", sqrts=" << sqrts << endl;
		return 0;
		if (y1<0) y1=0; if (y2<0)y2=0;
		
	}
	
	Inthelper_ktfact par; par.pt=pt; par.N1=N;
    par.x1=x1; par.x2=x2;
	par.N2=N2;
    if (scale<0)  
        par.scale=pt*pt;
    else
        par.scale=scale;
	gsl_function fun; fun.params=&par;
	fun.function=Inthelperf_ktfact_q;
	
	if (N2!=NULL) // Initialize interpolators only once
	{
		N->InitializeInterpolation(x1);
		N2->InitializeInterpolation(x2);
	}
	
	double maxq = std::max(3.5*pt, 35.0);
    if (maxq>80)
        maxq=80;

    
    double result, abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.001, maxq, 0, 0.05,		
		INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		//#pragma omp critical
		//cerr << "kt-factorization q integral failed at " << LINEINFO <<", pt=" << pt <<", x1=" << x1 <<", x2=" << x2 
		//	<< ", result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    result *= 2.0/(Cf()*SQR(pt));
    double alphas=Alphas(scale);
    result *= alphas;
    
    return result;
}

void testf_phi(){
    Inthelper_ktfact params;
    params.pt = 1.0;     // 示例值
    params.qt = 0.5;     // 示例值
    params.x1 = 0.01;    // 示例值
    params.x2 = 0.01;    // 示例值
    params.scale = 1.0; // 示例值
    params.N1 = &N;  // 初始化你的 UGD 对象
    params.N2 = &N2; // 如果不使用并行计算

    double phi = 0.5; // 示例角度
    double result = Inthelperf_ktfact_phi(phi, &params); // 调用函数
    std::cout << "Result at phi=" << phi << ": " << result << std::endl;

    string outputFileName1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/f_phi.txt";
    ofstream outputFile1(outputFileName1);

    for(double phi = -TMath::Pi();phi<=TMath::Pi();phi=phi+0.1){
        double result = Inthelperf_ktfact_phi(phi, &params);
        outputFile1 << std::scientific << std::setprecision(9) << phi << " " << result  << endl;
    }
    outputFile1.close();
}


void testf_qqt(){
    Inthelper_ktfact params;
    params.pt = 1.0;      // 示例值
    params.qt = 0.5;       // 会被 Inthelperf_ktfact_q 修改
    params.x1 = 0.01;      // 示例值
    params.x2 = 0.01;      // 示例值
    params.scale = 1.0;   // 示例值
    params.N1 = &N;    // 初始化你的 UGD 对象
    params.N2 = &N2;   // 如果不使用并行计算

    //double q = 0.5;  // 要计算的 qt 值
    //double result = Inthelperf_ktfact_q(q, &params);
    //std::cout << "Integral result at q=" << q << ": " << result << std::endl;

    string outputFileName1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/f_qqt.txt";
    ofstream outputFile1(outputFileName1);
    for (double pt=0.1; pt<=6; pt+=0.1){
        double result = Inthelperf_ktfact_q(pt, &params);
        outputFile1 << std::scientific << std::setprecision(9) << pt << " " << result  << endl;
    }
    outputFile1.close();
}


double testf_dNdydpt(double pt){
    double y = 0;
    double sqrts = 7000;
    Inthelper_ktfact params;
    params.pt = pt;      // 示例值
    //params.qt = 0.5;       // 会被 Inthelperf_ktfact_q 修改
    params.x1 = 0.01;      // 示例值
    params.x2 = 0.01;      // 示例值
    params.scale = 1.0;   // 示例值
    params.N1 = &N;    // 初始化你的 UGD 对象
    params.N2 = &N2;   // 如果不使用并行计算

    double x1 = params.pt*std::exp(-y)/sqrts;
	double x2 = params.pt*std::exp(y)/sqrts;

    
	gsl_function fun; fun.params=&params;
	fun.function=Inthelperf_ktfact_q;
    fun.function=Inthelperf_ktfact_q;

    N.InitializeInterpolation(x1);
	N2.InitializeInterpolation(x2);

    double maxq = std::max(3.5*params.pt, 35.0);

    double result, abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.001, maxq, 0, 0.05,		
		INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);

    result *= 2.0/(Cf()*SQR(pt));
    double alphas=Alphas(params.scale);
    result *= alphas;
    return result;

}


void testf_dNdydpt_answer(){

    string outputFileName1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/f_dNdydpt.txt";
    ofstream outputFile1(outputFileName1);
    for (double pt=0.1; pt<=6; pt+=0.1){
        double result = testf_dNdydpt(pt);
        outputFile1 << std::scientific << std::setprecision(9) << pt << " " << result  << endl;
    }
    outputFile1.close();

}

void test(){

    double sqrts = 7000;
    double minpt = 1;
    double maxpt = 6;
    double ptstep = 0.1;
    double y = 0 ;

    //Target

    //Probe

    cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
    cout << "# Using k_T factorization (gluon production); sigma02 = 1, S_T=1" << endl;
    cout << "# Gluon production (parton level)" << endl;
    cout << "# Probe: " << N2.GetString() << endl << "# Target: " << N.GetString() << endl;
    cout << "# NOTICE: results must be multiplied by (\\sigma_0/2)^2/S_T" << endl;
    cout << "# p_T   dN/(d^2 p_T dy) partonlevel   " << endl;


    string outputFileName1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/MVgamma_dNdy.txt";
    ofstream outputFile1(outputFileName1);

    for (double pt=minpt; pt<=maxpt; pt+=ptstep)
    {
        double partonresult = dHadronMultiplicity_dyd2pt_ktfact_parton(y, pt, sqrts, &N,&N2,-1);
        //double partonresult = xs.dHadronMultiplicity_dyd2pt_ktfact_parton(y, pt, sqrts, &N);
        outputFile1 << std::scientific << std::setprecision(9) << pt << " " << partonresult  << endl;
        //cout << pt << " " << partonresult << endl;            
    }
    outputFile1.close();
}











//=================================================main=====================================================//

int main(){
    string outputFileName = "xg1e-4_MV.txt";
    ofstream outputFile(outputFileName);

    TCanvas *c = new TCanvas("c1","c1",500,700);
    c->Divide(4,4);
    double y=0.1;
    double x0 = 0.0001;
    double xbj = x0; //N.X0()*std::exp(-y);
    N.InitializeInterpolation(xbj);

    TGraph* graph = new TGraph();
    graph->SetNameTitle("Dipole Amplitude", "Dipole Amplitude; r [GeV^{-1}]; N(r,x = 1e-2)");

    double minr = N.MinR(); double maxr=N.MaxR();

    int i =1;
    for (double r=minr; r<maxr; r*=1.04725)
    {
        graph->SetPoint(i, r, N.N(r, xbj));
        i++;
    }

    c->cd(1);
    TH1D *h = new TH1D("h","Dipole amplitude;r[GeV^{-1}];N(r,x)",1,1e-2,20);
    h->GetYaxis()->SetRangeUser(1e-4,1);
    h->SetStats(0);
    h->Draw();
    graph->Draw("SAME");


    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("S(k)", "S(k); k [GeV]; S(k)");

    for (double k=0.01; k<maxr; k=k+0.01)
    {
        graph1->SetPoint(i, k, S_k(k,xbj,-1));
        i++;
    }
    c->cd(2);
    TH1D *h1 = new TH1D("h","S(k);k[GeV];S(k)",1,1e-1,15);
    h1->GetYaxis()->SetRangeUser(1e-4,120);
    h1->SetStats(0);
    h1->Draw();
    graph1->Draw("SAME");

    TCanvas *c1 = new TCanvas("c1","c1",500,700);
    
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    for (double k=1; k<10; k=k+0.01)
    {
        graph2->SetPoint(i, k*k, xg(xbj,k));
        outputFile << std::scientific << std::setprecision(9) << k*k << " " << xg(xbj,k)  << endl;
        i++;
    }

    TH1D *h2 = new TH1D("h","xg;Q^{2}[GeV^{2}];xg(x,Q^{2})",1,0,1e+3);
    h2->GetYaxis()->SetRangeUser(0,1e+2);
    h2->SetStats(0);
    h2->Draw();
    graph2->Draw("SAME");

    outputFile.close();
   
    return 0;
}

double x_g(double k){
    double sqrt_s = 7000;
    double y = 0.0;
    return (k/sqrt_s) * TMath::Exp(-y);
  }

void dNOverdydpT(){

    string outputFileName = "update_MVgamma_dNOverdydpT_y=0_0.txt";
    ofstream outputFile(outputFileName);
    double sigma0Over2 = 16.45;
    double sigma_inel = 42;
    double A = pow(2*TMath::Pi(),2);

    for(double k = 1;k<= 10;k=k+0.1){
        double xbj = x_g(k);
        N.InitializeInterpolation(xbj);
        double result = (sigma0Over2/sigma_inel) * 1/A * xg(xbj, k) * S_k(k,xbj,-1);
        outputFile << std::scientific << std::setprecision(9) << k << " " << result  << endl;

    }
    outputFile.close();

}

void drawInvYield(){
    string datafile1 = "update_MVgamma_dNOverdydpT_y=0_0.txt"; 
    string datafile2 = "MVgamma_dNOverdydpT_y=3_3.txt"; 
    string datafile3 = "data_MVgamma_dNOverdydpT_y=3_3.txt"; 

    string datafile4 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/MVgamma_dNdy.txt"; 
    
    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    TGraph* graph3 = new TGraph();
    graph3->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    TGraph* graph4 = new TGraph();
    graph4->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);
    ifstream inputFile3(datafile3);
    ifstream inputFile4(datafile4);

    double k, result;
    int i = 1;
    while (inputFile1 >> k >> result)
    {
        graph1->SetPoint(i, k, result*5.3);
        i++;
    }
    inputFile1.close();


    while (inputFile2 >> k >> result)
    {
        graph2->SetPoint(i, k, result);
        i++;
    }
    inputFile2.close();

    while (inputFile3 >> k >> result)
    {
        graph3->SetPoint(i, k, result*10);
        i++;
    }
    inputFile3.close();

    while (inputFile4 >> k >> result)
    {
        graph4->SetPoint(i, k, result*10);
        i++;
    }
    inputFile4.close();

    graph1->SetLineColorAlpha(kBlack,1); 
    graph1->SetLineStyle(2);      
    graph1->SetLineWidth(6);   
    graph2->SetLineColorAlpha(kBlue,1); 
    graph2->SetLineStyle(1);      
    graph2->SetLineWidth(6);  
    
    graph3->SetLineColorAlpha(kRed,1); 
    graph3->SetLineStyle(1);      
    graph3->SetLineWidth(6); 

    graph4->SetLineColorAlpha(kGreen,1); 
    graph4->SetLineStyle(1);      
    graph4->SetLineWidth(6); 

    TH1D *h = new TH1D("h","dN/dyd^{2}p_{T}; p_{T} [GeV/c]; dN/dyd^{2}p_{T} [1/GeV^{2}]",1,1,10);
    h->GetYaxis()->SetRangeUser(1e-8,1);
    h->SetStats(0);
    h->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    graph3->Draw("same");

    graph4->Draw("same");

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "update : MV^{#gamma},y = 3.3 ", "lp"); 
    legend->AddEntry(graph2, "MV^{#gamma},y = 3.3 ", "lp");
    legend->AddEntry(graph3, "Fig.6 : MV^{#gamma},y = 3.3 ", "lp");
    legend->AddEntry(graph4, "simulation : MV^{#gamma},y = 0.0 ", "lp");
    legend->Draw("SAME");

}

void drawGDF(){
    string datafile1 = "xg1e-2_MV.txt"; 
    string datafile2 = "data_MV_x=1e-2.txt"; 

    string datafile3 = "xg1e-4_MV.txt"; 
    string datafile4 = "data_MV_x=1e-4.txt";

    string datafile5 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/x=1e-4_MV_xg.txt";

   
    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    TGraph* graph3 = new TGraph();
    graph3->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    TGraph* graph4 = new TGraph();
    graph4->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
   
    TGraph* graph5 = new TGraph();
    graph5->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    ifstream inputFile3(datafile3);
    ifstream inputFile4(datafile4);

    ifstream inputFile5(datafile5);
   

    double k, xg;
    int i = 1;
    while (inputFile1 >> k >> xg)
    {
        graph1->SetPoint(i, k, xg);//*pow(1.73,1.73)
        i++;
    }
    inputFile1.close();

    while (inputFile2 >> k >> xg)
    {
        graph2->SetPoint(i, k, xg);
        i++;
    }
    inputFile2.close();

    while (inputFile3 >> k >> xg)
    {
        graph3->SetPoint(i, k, xg);//*pow(1.73,1.73)
        i++;
    }
    inputFile3.close();

    while (inputFile4 >> k >> xg)
    {
        graph4->SetPoint(i, k, xg);
        i++;
    }
    inputFile4.close();

    while (inputFile5 >> k >> xg)
    {
        graph5->SetPoint(i, k, xg);//*pow(1.73,1.73)
        i++;
    }
    inputFile5.close();

  

    graph1->SetLineColorAlpha(kRed,1); 
    graph1->SetLineStyle(2);      
    graph1->SetLineWidth(6);   

    graph2->SetLineColorAlpha(kRed,0.3); 
    graph2->SetLineStyle(2);      
    graph2->SetLineWidth(6); 

    graph3->SetLineColorAlpha(kBlue,1); 
    graph3->SetLineStyle(2);      
    graph3->SetLineWidth(6);   

    graph4->SetLineColorAlpha(kBlue,0.3); 
    graph4->SetLineStyle(2);      
    graph4->SetLineWidth(6); 


    graph5->SetLineColorAlpha(kBlack,1); 
    graph5->SetLineStyle(2);      
    graph5->SetLineWidth(6); 

 

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "Model:MV^{},x = 1e-2,adjoint rep", "l");   
    legend->AddEntry(graph2, "data:MV^{},x = 1e-2,adjoint rep", "l");  

    legend->AddEntry(graph3, "Model:MV^{},x = 1e-4,adjoint rep", "l");   
    legend->AddEntry(graph4, "data:MV^{},x = 1e-4,adjoint rep", "l");  
    //legend->AddEntry(graph5, "data:MV^{},x = 1e-4,adjoint rep", "l");  

    TH1D *h2 = new TH1D("h","Gluon distribution function;Q^{2}[GeV^{2}];xg(x,Q^{2})",1,12,1e+3);
    h2->GetYaxis()->SetRangeUser(5,1e+2);
    h2->SetStats(0);
    h2->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    graph3->Draw("same");
    graph4->Draw("same");
    //graph5->Draw("same");
    legend->Draw("same");
}

void dNdydpT(){
    string datafile1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/MVgamma_dNdy.txt"; 
    string datafile2 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/MVgamma_dNdy.txt"; 

    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    double k, xg;
    int i = 1;
    while (inputFile1 >> k >> xg)
    {
        graph1->SetPoint(i, k, xg);
        i++;
    }
    inputFile1.close();

    while (inputFile2 >> k >> xg)
    {
        graph2->SetPoint(i, k, xg);
        i++;
    }
    inputFile2.close();

    graph1->SetLineColorAlpha(kRed,1); 
    graph1->SetLineStyle(1);      
    graph1->SetLineWidth(6);   

    graph2->SetLineColorAlpha(kBlue,1); 
    graph2->SetLineStyle(2);      
    graph2->SetLineWidth(6); 

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "data my code:MV^{#gamma},adjoint rep", "l");   
    legend->AddEntry(graph2, "data H.Mantysaari code:MV^{#gamma},adjoint rep", "l");  

    TH1D *h = new TH1D("h","d#sigma^{p+p->g}/dyd^{2}p_{T}; p_{T} [GeV/c]; d#sigma^{p+p->g}/dyd^{2}p_{T} [1/GeV^{2}]",1,1,10);
    h->GetYaxis()->SetRangeUser(1e-8,1);
    h->SetStats(0);
    h->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    legend->Draw();


}

void f_phi(){
    string datafile1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/f_phi.txt"; 
    string datafile2 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/f_phi.txt"; 

    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    double k, xg;
    int i = 1;
    while (inputFile1 >> k >> xg)
    {
        graph1->SetPoint(i, k, xg);
        i++;
    }
    inputFile1.close();

    while (inputFile2 >> k >> xg)
    {
        graph2->SetPoint(i, k, xg);
        i++;
    }
    inputFile2.close();

    graph1->SetLineColorAlpha(kRed,1); 
    graph1->SetLineStyle(1);      
    graph1->SetLineWidth(6);   

    graph2->SetLineColorAlpha(kBlue,1); 
    graph2->SetLineStyle(2);      
    graph2->SetLineWidth(6); 

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "data my code:MV^{#gamma},adjoint rep f(#phi)", "l");   
    legend->AddEntry(graph2, "data H.Mantysaari code:MV^{#gamma},adjoint rep f(#phi)", "l");  

    TH1D *h = new TH1D("h","f(#phi); #phi; f(#phi)",1,-TMath::Pi(),TMath::Pi());
    h->GetYaxis()->SetRangeUser(1e-4,1e-2);
    h->SetStats(0);
    h->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    legend->Draw();

}

void f_qqt(){
    string datafile1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/f_qqt.txt"; 
    string datafile2 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/f_qqt.txt"; 

    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    double k, xg;
    int i = 1;
    while (inputFile1 >> k >> xg)
    {
        graph1->SetPoint(i, k, xg);
        i++;
    }
    inputFile1.close();

    while (inputFile2 >> k >> xg)
    {
        graph2->SetPoint(i, k, xg);
        i++;
    }
    inputFile2.close();

    graph1->SetLineColorAlpha(kRed,1); 
    graph1->SetLineStyle(1);      
    graph1->SetLineWidth(6);   

    graph2->SetLineColorAlpha(kBlue,1); 
    graph2->SetLineStyle(2);      
    graph2->SetLineWidth(6); 

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "data my code:MV^{#gamma},adjoint rep q(q_{t})", "l");   
    legend->AddEntry(graph2, "data H.Mantysaari code:MV^{#gamma},adjoint rep q(q_{t})", "l");  

    TH1D *h = new TH1D("h","q(q_{t}); q_{T} GeV/c; q(q_{t})",1,0.1,6);
    h->GetYaxis()->SetRangeUser(1e-8,1);
    h->SetStats(0);
    h->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    legend->Draw();

}

void run(){
    cout<<"successful"<<endl;
}

double SingleInclusive::dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2, double scale )
{
	double x1 = pt*std::exp(-y)/sqrts;
	double x2 = pt*std::exp(y)/sqrts;
	double y1 = std::log(N->X0()/x1);
	double y2;
	if (N2==NULL)
		y2 = std::log(N->X0()/x2);
	else
		y2 = std::log(N2->X0()/x2);
	
	if (y1<0 or y2<0)
	{
		cerr << "Evolution variables y1=" << y1 <<", y2=" << y2 <<" is not possible to compute, too large x. pt=" << pt << ", y=" << y << ", sqrts=" << sqrts <<" " << LINEINFO << endl;
		return 0;
		if (y1<0) y1=0; if (y2<0)y2=0;
		
	}
	
	Inthelper_ktfact par; par.pt=pt; par.N1=N;
    par.x1=x1; par.x2=x2;
	par.N2=N2;
    if (scale<0)  
        par.scale=pt*pt;
    else
        par.scale=scale;
	gsl_function fun; fun.params=&par;
	fun.function=Inthelperf_ktfact_q;
	
	if (N2!=NULL) // Initialize interpolators only once
	{
		N->InitializeInterpolation(x1);
		N2->InitializeInterpolation(x2);
	}
	
	double maxq = std::max(3.5*pt, 35.0);
    if (maxq>80)
        maxq=80;

    
    double result, abserr;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.001, maxq, 0, 0.05,		
		INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		//#pragma omp critical
		//cerr << "kt-factorization q integral failed at " << LINEINFO <<", pt=" << pt <<", x1=" << x1 <<", x2=" << x2 
		//	<< ", result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    result *= 2.0/(qcd.Cf()*SQR(pt));
    double alphas=qcd.Alphas(scale);
    result *= alphas;
    
    return result;
}

void drawdNdydpTPi0(){
    string datafile1 = "/Users/suyoupeng/Downloads/rcbkdipole-master/src/data_MVePi0.txt"; 
    string datafile2 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/MVe_pi0.txt"; 

    TGraph* graph1 = new TGraph();
    graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    TGraph* graph2 = new TGraph();
    graph2->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");

    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    double k, xg;
    int i = 1;
    while (inputFile1 >> k >> xg)
    {
        graph1->SetPoint(i, k, xg);
        i++;
    }
    inputFile1.close();

    while (inputFile2 >> k >> xg)
    {
        graph2->SetPoint(i, k, xg*16.36*16.36*1/70*2*TMath::Pi());//*
        i++;
    }
    inputFile2.close();

    graph1->SetLineColorAlpha(kRed,1); 
    graph1->SetLineStyle(1);      
    graph1->SetLineWidth(6);   

    graph2->SetLineColorAlpha(kBlue,1); 
    graph2->SetLineStyle(2);      
    graph2->SetLineWidth(6); 

    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(graph1, "k_{T} & kkp:data from FIG. 9:MV^{e},adjoint rep ", "l");   
    legend->AddEntry(graph2, "(x (#sigma_{0}/2)^{2}x 1/#sigma_{inel} x 2#pi ) k_{T} & kkp:H.Mantysaari code:MV^{e},adjoint rep ", "l");  

    TH1D *h = new TH1D("h","dN/d^{2}p_{T}dy p+p->#pi^{0}+X #sqrt{s}=7000 GeV ; p_{T} GeV/c; dN/d^{2}p_{T}dy 1/GeV^{2}",1,1,8);
    h->GetYaxis()->SetRangeUser(1e-6,1);
    h->SetStats(0);
    h->Draw();
    graph1->Draw("same");
    graph2->Draw("same");
    legend->Draw();

}


void Dipole_UGD_2D() {
    // Parameters
    double x0 = 0.01;
    double Y_min = 0;
    double Y_step = 0.1;
    double Y_max = 16.2;
    double q_max = 10;
    double q_step = 0.2; // Finer step for better resolution

    // Create containers for data
    std::vector<double> x_values, q_values, ugd_values;

    string outputFileName1 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/UGD_x_pt.txt";
    ofstream outputFile1(outputFileName1);

    // Calculate UGD values
    for (double Y = Y_min; Y <= Y_max; Y += Y_step) {
        double x = x0 / pow(10, Y);
        for (double q = 0; q <= q_max; q += q_step) {
            double ugd = Dipole_UGD(q, x, q*q, 1); // Assuming Dipole_UGD returns a double
            x_values.push_back(x); // Using log10(x) for better visualization
            q_values.push_back(q);
            ugd_values.push_back(ugd);
            outputFile1 << std::scientific << std::setprecision(9) << x << " " << q << " "<< ugd <<endl;
        }
    }
    outputFile1.close();
}

void Dipole_UGD_x(){
    string datafile1 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/UGD_x_pt.txt"; 
    string datafile2 = "/Users/suyoupeng/Downloads/amplitudelib-2/build/bin/auther_UGD_x_pt.txt"; 

    TGraph2D *g2d = new TGraph2D();
    g2d->SetTitle("Dipole UGD (2D Scatter);log_{10}(x); p_{T} [GeV/c];#phi_{p}(x,p_{T})");

    TGraph2D *g2d1 = new TGraph2D();
    g2d1->SetTitle("Dipole UGD (2D Scatter);log_{10}(x); p_{T} [GeV/c];#phi_{p}(x,p_{T})");

    //TGraph* graph1 = new TGraph();
    //graph1->SetNameTitle("xg", "xg; Q^{2} [GeV^{2}]; xg(x,Q^{2})");
    ifstream inputFile1(datafile1);
    ifstream inputFile2(datafile2);

    g2d->SetMarkerStyle(2);  // 设置点的样式
    g2d->SetMarkerColor(kBlue); // 设置点的颜色

    g2d1->SetMarkerStyle(5);  // 设置点的样式
    g2d1->SetMarkerColor(kRed); // 设置点的颜色

    

    double x, p, ugd;  // dummy变量用于存储不用的第二列
    int i = 0;            // 建议从0开始计数（ROOT的SetPoint索引通常从0开始）
    while (inputFile1 >> x >> p >> ugd) {
        g2d->SetPoint(i, log10(x), p,ugd);
        i++;
    }
    inputFile1.close();

    double x1, p1, ugd1; 
    int j = 0;   
    while (inputFile2 >> x1 >> p1 >> ugd1) {
       
            g2d1->SetPoint(j, log10(x1), p1,ugd1);
            j++;
        
           
    }
    inputFile2.close();

    
    
    
    
    g2d->Draw();
    g2d1->Draw("PL same");
    TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89); 
    legend->AddEntry(g2d, "My result:dipole unintegrated gluon distribution ", "lp");   
    legend->AddEntry(g2d1, "H.Mantysaari:dipole unintegrated gluon distribution ", "lp");  
    legend->Draw();

}



