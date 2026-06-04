#ifdef __CLING__
#endif

#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TString.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include "TMath.h"
#include "TLatex.h"

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooProdPdf.h>
#include <RooAbsReal.h>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooMinimizer.h>
#include <RooArgSet.h>
#include <RooCurve.h>
#include <RooHist.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooCrystalBall.h>
#include <RooGenericPdf.h>
#include <RooChebychev.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooExponential.h>
#include <RooConstVar.h>
#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"
#include "RooMath.h"

using namespace RooFit;

// --------------------------- qGaus right-tail ---------------------------
class RooQGausRightTail : public RooAbsPdf {
public:
  RooQGausRightTail() {}
  RooQGausRightTail(const char* name, const char* title,
                    RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau)
    : RooAbsPdf(name, title),
      x("x","x",this,_x),
      mu("mu","mu",this,_mu),
      sigma("sigma","sigma",this,_sigma),
      tau("tau","tau",this,_tau)
  {}

  RooQGausRightTail(const RooQGausRightTail& other, const char* name=nullptr)
    : RooAbsPdf(other, name),
      x("x",this,other.x),
      mu("mu",this,other.mu),
      sigma("sigma",this,other.sigma),
      tau("tau",this,other.tau)
  {}

  TObject* clone(const char* newname) const override { return new RooQGausRightTail(*this,newname); }
  ~RooQGausRightTail() override {}

protected:
  RooRealProxy x, mu, sigma, tau;

  Double_t evaluate() const override {
    const double xx = (double)x;
    const double m  = (double)mu;
    const double s  = std::max(1e-6, (double)sigma);
    const double t  = std::max(1e-6, (double)tau);

    const double x0 = m + t;
    if (xx <= x0) {
      const double u = (xx - m)/s;
      return std::exp(-0.5*u*u);
    } else {
      const double u0 = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double expo = std::exp( -t*(xx - m - t)/(s*s) );
      return cont * expo;
    }
  }
};

// --------------------------- qGaus left-tail ---------------------------
class RooQGausLeftTail : public RooAbsPdf {
public:
  RooQGausLeftTail() {}
  RooQGausLeftTail(const char* name, const char* title,
                   RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau)
    : RooAbsPdf(name, title),
      x("x","x",this,_x),
      mu("mu","mu",this,_mu),
      sigma("sigma","sigma",this,_sigma),
      tau("tau","tau",this,_tau)
  {}

  RooQGausLeftTail(const RooQGausLeftTail& other, const char* name=nullptr)
    : RooAbsPdf(other, name),
      x("x",this,other.x),
      mu("mu",this,other.mu),
      sigma("sigma",this,other.sigma),
      tau("tau",this,other.tau)
  {}

  TObject* clone(const char* newname) const override { return new RooQGausLeftTail(*this,newname); }
  ~RooQGausLeftTail() override {}

protected:
  RooRealProxy x, mu, sigma, tau;

  Double_t evaluate() const override {
    const double xx = (double)x;
    const double m  = (double)mu;
    const double s  = std::max(1e-6, (double)sigma);
    const double t  = std::max(1e-6, (double)tau);

    const double x0 = m - t;
    if (xx >= x0) {
      const double u = (xx - m)/s;
      return std::exp(-0.5*u*u);
    } else {
      const double u0   = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double k    = t/(s*s);
      return cont * std::exp(+k*(xx - x0));
    }
  }
};

// --------------------------- Gaussian core + exponential RIGHT tail ---------------------------
class RooGausExpRightTail : public RooAbsPdf {
public:
  RooGausExpRightTail() {}
  RooGausExpRightTail(const char* name, const char* title,
                      RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau)
    : RooAbsPdf(name, title),
      x("x","x",this,_x),
      mu("mu","mu",this,_mu),
      sigma("sigma","sigma",this,_sigma),
      tau("tau","tau",this,_tau)
  {}

  RooGausExpRightTail(const RooGausExpRightTail& other, const char* name=nullptr)
    : RooAbsPdf(other, name),
      x("x",this,other.x),
      mu("mu",this,other.mu),
      sigma("sigma",this,other.sigma),
      tau("tau",this,other.tau)
  {}

  TObject* clone(const char* newname) const override { return new RooGausExpRightTail(*this,newname); }
  ~RooGausExpRightTail() override {}

protected:
  RooRealProxy x, mu, sigma, tau;

  Double_t evaluate() const override {
    const double xx = (double)x;
    const double m  = (double)mu;
    const double s  = std::max(1e-9, (double)sigma);
    const double t  = std::max(1e-9, (double)tau);

    const double x0 = m + t;
    if (xx <= x0) {
      const double u = (xx - m)/s;
      return std::exp(-0.5*u*u);
    } else {
      const double u0   = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double k    = t/(s*s);
      return cont * std::exp(-k*(xx - x0));
    }
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* /*rangeName*/) const override {
    if (matchArgs(allVars, analVars, x)) return 1;
    return 0;
  }

  Double_t analyticalIntegral(Int_t code, const char* rangeName) const override {
    if (code != 1) return 0.0;

    const double m  = (double)mu;
    const double s  = std::max(1e-9, (double)sigma);
    const double t  = std::max(1e-9, (double)tau);
    const double x0 = m + t;

    const double a = x.min(rangeName);
    const double b = x.max(rangeName);

    auto intG = [&](double lo, double hi) {
      const double inv = 1.0/(std::sqrt(2.0)*s);
      const double Ehi = TMath::Erf((hi - m)*inv);
      const double Elo = TMath::Erf((lo - m)*inv);
      return s*std::sqrt(TMath::Pi()/2.0) * (Ehi - Elo);
    };

    auto intE = [&](double lo, double hi) {
      const double u0   = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double k    = t/(s*s);
      const double eLo  = std::exp(-k*(lo - x0));
      const double eHi  = std::exp(-k*(hi - x0));
      return cont*(1.0/k)*(eLo - eHi);
    };

    if (b <= x0) return intG(a,b);
    if (a >= x0) return intE(a,b);
    return intG(a, x0) + intE(x0, b);
  }
};

// --------------------------- Gaussian core + exponential LEFT tail ---------------------------
class RooGausExpLeftTail : public RooAbsPdf {
public:
  RooGausExpLeftTail() {}
  RooGausExpLeftTail(const char* name, const char* title,
                     RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau)
    : RooAbsPdf(name, title),
      x("x","x",this,_x),
      mu("mu","mu",this,_mu),
      sigma("sigma","sigma",this,_sigma),
      tau("tau","tau",this,_tau)
  {}

  RooGausExpLeftTail(const RooGausExpLeftTail& other, const char* name=nullptr)
    : RooAbsPdf(other, name),
      x("x",this,other.x),
      mu("mu",this,other.mu),
      sigma("sigma",this,other.sigma),
      tau("tau",this,other.tau)
  {}

  TObject* clone(const char* newname) const override { return new RooGausExpLeftTail(*this,newname); }
  ~RooGausExpLeftTail() override {}

protected:
  RooRealProxy x, mu, sigma, tau;

  Double_t evaluate() const override {
    const double xx = (double)x;
    const double m  = (double)mu;
    const double s  = std::max(1e-9, (double)sigma);
    const double t  = std::max(1e-9, (double)tau);

    const double x0 = m - t;
    if (xx >= x0) {
      const double u = (xx - m)/s;
      return std::exp(-0.5*u*u);
    } else {
      const double u0   = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double k    = t/(s*s);
      return cont * std::exp(+k*(xx - x0));
    }
  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* /*rangeName*/) const override {
    if (matchArgs(allVars, analVars, x)) return 1;
    return 0;
  }

  Double_t analyticalIntegral(Int_t code, const char* rangeName) const override {
    if (code != 1) return 0.0;

    const double m  = (double)mu;
    const double s  = std::max(1e-9, (double)sigma);
    const double t  = std::max(1e-9, (double)tau);
    const double x0 = m - t;

    const double a = x.min(rangeName);
    const double b = x.max(rangeName);

    auto intG = [&](double lo, double hi) {
      const double inv = 1.0/(std::sqrt(2.0)*s);
      const double Ehi = TMath::Erf((hi - m)*inv);
      const double Elo = TMath::Erf((lo - m)*inv);
      return s*std::sqrt(TMath::Pi()/2.0) * (Ehi - Elo);
    };

    auto intE = [&](double lo, double hi) {
      const double u0   = t/s;
      const double cont = std::exp(-0.5*u0*u0);
      const double k    = t/(s*s);
      const double eHi  = std::exp(+k*(hi - x0));
      const double eLo  = std::exp(+k*(lo - x0));
      return cont*(1.0/k)*(eHi - eLo);
    };

    if (a >= x0) return intG(a,b);
    if (b <= x0) return intE(a,b);
    return intE(a, x0) + intG(x0, b);
  }
};

// ---------------- helpers ----------------
static void IntegralAndErrorH2Window(const TH2* h2,
                                     double xlow, double xhi,
                                     double ylow, double yhi,
                                     double& val, double& err)
{
  val = 0.0; err = 0.0;
  if (!h2) return;

  auto ax = h2->GetXaxis();
  auto ay = h2->GetYaxis();
  if (!ax || !ay) return;

  const double eps = 1e-12;

  int bx1 = ax->FindBin(xlow);
  int bx2 = ax->FindBin(xhi - eps);
  int by1 = ay->FindBin(ylow);
  int by2 = ay->FindBin(yhi - eps);

  bx1 = std::max(1, std::min(bx1, ax->GetNbins()));
  bx2 = std::max(1, std::min(bx2, ax->GetNbins()));
  by1 = std::max(1, std::min(by1, ay->GetNbins()));
  by2 = std::max(1, std::min(by2, ay->GetNbins()));

  if (bx2 < bx1) std::swap(bx1, bx2);
  if (by2 < by1) std::swap(by1, by2);

  val = h2->IntegralAndError(bx1, bx2, by1, by2, err);
}

struct ValErr { double v{0.0}; double e{0.0}; };

static inline ValErr AddVE(const ValErr& a, const ValErr& b) { return {a.v + b.v, std::sqrt(a.e*a.e + b.e*b.e)}; }
static inline ValErr SubVE(const ValErr& a, const ValErr& b) { return {a.v - b.v, std::sqrt(a.e*a.e + b.e*b.e)}; }
static inline ValErr ScaleVE(const ValErr& a, double s)      { return {s*a.v, std::abs(s)*a.e}; }
static inline ValErr MakeSigFromAandBkg(const ValErr& A, const ValErr& Bkg) { return {A.v - Bkg.v, std::sqrt(A.e*A.e + Bkg.e*Bkg.e)}; }
static inline ValErr MakeFracSigFromAandBkg(const ValErr& A, const ValErr& Bkg) {
  ValErr out;
  if (A.v <= 0.0) return out;
  out.v = 1.0 - Bkg.v / A.v;
  const double dfdA = Bkg.v / (A.v * A.v);
  const double dfdB = -1.0 / A.v;
  out.e = std::sqrt(dfdA*dfdA * A.e*A.e + dfdB*dfdB * Bkg.e*Bkg.e);
  return out;
}
static inline ValErr IntegralVE(const TH2* h2, double x1,double x2,double y1,double y2) {
  ValErr out;
  IntegralAndErrorH2Window(h2, x1,x2,y1,y2, out.v, out.e);
  return out;
}

static int findBinForEdgeLow(const TAxis* ax, double edge) {
  int b = ax->FindBin(edge + 1e-9);
  b = std::max(1, std::min(b, ax->GetNbins()));
  return b;
}
static int findBinForEdgeHigh(const TAxis* ax, double edge) {
  int b = ax->FindBin(edge - 1e-9);
  b = std::max(1, std::min(b, ax->GetNbins()));
  return b;
}

static TGraphErrors* MakeRatioGraph(RooPlot* frame,
                                    const char* dataName,
                                    const char* curveName)
{
  auto* hData  = dynamic_cast<RooHist*>(frame->findObject(dataName));
  auto* cModel = dynamic_cast<RooCurve*>(frame->findObject(curveName));
  if (!hData || !cModel) return nullptr;

  const int n = hData->GetN();
  auto* gr = new TGraphErrors(n);
  gr->SetName(Form("ratio_%s_%s", dataName, curveName));

  for (int i = 0; i < n; ++i) {
    double xx, yy;
    hData->GetPoint(i, xx, yy);
    const double ey = hData->GetErrorY(i);
    const double yfit = cModel->Eval(xx);
    if (yfit <= 0) {
      gr->SetPoint(i, xx, 0.0);
      gr->SetPointError(i, 0.0, 0.0);
      continue;
    }
    gr->SetPoint(i, xx, yy / yfit);
    gr->SetPointError(i, 0.0, ey / yfit);
  }
  return gr;
}

static double GetPropErrSafe(const RooAbsReal* v, const RooFitResult* res)
{
  if (!v || !res) return 0.0;
  double e = 0.0;
  try { e = v->getPropagatedError(*res); } catch (...) { e = 0.0; }
  if (!std::isfinite(e)) e = 0.0;
  return e;
}

static TH1D* MakeRatioHist(const TH1D* num, const TH1D* den, const char* name, const char* title)
{
  if (!num || !den) return nullptr;
  auto* h = (TH1D*)num->Clone(name);
  h->SetTitle(title);
  h->Reset("ICES");
  h->Sumw2();
  for (int b = 1; b <= num->GetNbinsX(); ++b) {
    const double a  = num->GetBinContent(b);
    const double ea = num->GetBinError(b);
    const double c  = den->GetBinContent(b);
    const double ec = den->GetBinError(b);
    if (c <= 0.0) {
      h->SetBinContent(b, 0.0);
      h->SetBinError(b, 0.0);
      continue;
    }
    const double r = a / c;
    double er = 0.0;
    if (a > 0.0) er = std::abs(r) * std::sqrt((ea*ea)/(a*a) + (ec*ec)/(c*c));
    else         er = std::abs(ea / c);
    h->SetBinContent(b, r);
    h->SetBinError(b, er);
  }
  return h;
}

// --------------------------- ABCD drawing helpers ---------------------------
static void DrawSingleBox(double x1, double x2, double y1, double y2,
                          int color, int style = 1, int width = 4)
{
  TLine* l1 = new TLine(x1, y1, x2, y1);
  TLine* l2 = new TLine(x2, y1, x2, y2);
  TLine* l3 = new TLine(x2, y2, x1, y2);
  TLine* l4 = new TLine(x1, y2, x1, y1);
  for (auto* l : {l1, l2, l3, l4}) {
    l->SetLineColor(color);
    l->SetLineStyle(style);
    l->SetLineWidth(width);
    l->Draw("SAME");
  }
}

static void DrawABCDOverlays(double xlow, double xhi, double ylow, double yhi,
                             double xSB1_lo, double xSB1_hi,
                             double xSB2_lo, double xSB2_hi,
                             double ySB1_lo, double ySB1_hi,
                             double ySB2_lo, double ySB2_hi)
{
  DrawSingleBox(xlow, xhi, ylow, yhi, kRed+1, 1, 5);
  DrawSingleBox(xlow, xhi, ySB1_lo, ySB1_hi, kBlue+1, 2, 4);
  DrawSingleBox(xlow, xhi, ySB2_lo, ySB2_hi, kBlue+1, 2, 4);
  DrawSingleBox(xSB1_lo, xSB1_hi, ylow, yhi, kGreen+2, 2, 4);
  DrawSingleBox(xSB2_lo, xSB2_hi, ylow, yhi, kGreen+2, 2, 4);
  DrawSingleBox(xSB1_lo, xSB1_hi, ySB1_lo, ySB1_hi, kMagenta+1, 3, 4);
  DrawSingleBox(xSB1_lo, xSB1_hi, ySB2_lo, ySB2_hi, kMagenta+1, 3, 4);
  DrawSingleBox(xSB2_lo, xSB2_hi, ySB1_lo, ySB1_hi, kMagenta+1, 3, 4);
  DrawSingleBox(xSB2_lo, xSB2_hi, ySB2_lo, ySB2_hi, kMagenta+1, 3, 4);
}

static TCanvas* MakeABCDCanvas(TH2D* h2,
                               const TString& cname,
                               const TString& ctitle,
                               double xlow, double xhi,
                               double ylow, double yhi,
                               double xSB1_lo, double xSB1_hi,
                               double xSB2_lo, double xSB2_hi,
                               double ySB1_lo, double ySB1_hi,
                               double ySB2_lo, double ySB2_hi,
                               bool drawSidebands)
{
  if (!h2) return nullptr;
  TCanvas* c = new TCanvas(cname, ctitle, 900, 780);
  c->cd();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetLogz();
  h2->SetTitle(";m_{#Lambda1} (GeV/c^{2});m_{#Lambda2} (GeV/c^{2})");
  h2->Draw("COLZ");
  DrawSingleBox(xlow, xhi, ylow, yhi, kRed+1, 1, 5);
  if (drawSidebands) {
    DrawABCDOverlays(xlow, xhi, ylow, yhi,
                     xSB1_lo, xSB1_hi, xSB2_lo, xSB2_hi,
                     ySB1_lo, ySB1_hi, ySB2_lo, ySB2_hi);
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.040);
    lat.SetTextFont(42);
    lat.SetTextColor(kRed+1);    lat.DrawLatex(0.73, 0.90, "A: signal box");
    lat.SetTextColor(kBlue+1);   lat.DrawLatex(0.73, 0.84, "B: y-sidebands");
    lat.SetTextColor(kGreen+2);  lat.DrawLatex(0.73, 0.78, "C: x-sidebands");
    lat.SetTextColor(kMagenta+1);lat.DrawLatex(0.73, 0.72, "D: corners");
  }
  return c;
}

// --------------------------- 1D background builder ---------------------------
static std::unique_ptr<RooAbsPdf> buildBkg1D(int bkgOpt,
                                             RooRealVar& v,
                                             const TString& name,
                                             std::vector<std::unique_ptr<RooAbsArg>>& keep)
{
  if (bkgOpt == 0) {
    auto c0 = std::make_unique<RooRealVar>(Form("b2c0_%s", name.Data()), "b2c0", 0.5, 0.0, 10.0);
    auto c1 = std::make_unique<RooRealVar>(Form("b2c1_%s", name.Data()), "b2c1", 0.5, 0.0, 10.0);
    auto c2 = std::make_unique<RooRealVar>(Form("b2c2_%s", name.Data()), "b2c2", 0.5, 0.0, 10.0);
    auto pdf = std::make_unique<RooBernstein>(name, name, v, RooArgList(*c0, *c1, *c2));
    keep.emplace_back(std::move(c0));
    keep.emplace_back(std::move(c1));
    keep.emplace_back(std::move(c2));
    return pdf;
  }
  if (bkgOpt == 1) {
    auto c0 = std::make_unique<RooRealVar>(Form("b3c0_%s", name.Data()), "b3c0", 0.5, 0.0, 10.0);
    auto c1 = std::make_unique<RooRealVar>(Form("b3c1_%s", name.Data()), "b3c1", 0.5, 0.0, 10.0);
    auto c2 = std::make_unique<RooRealVar>(Form("b3c2_%s", name.Data()), "b3c2", 0.5, 0.0, 10.0);
    auto c3 = std::make_unique<RooRealVar>(Form("b3c3_%s", name.Data()), "b3c3", 0.5, 0.0, 10.0);
    auto pdf = std::make_unique<RooBernstein>(name, name, v, RooArgList(*c0, *c1, *c2, *c3));
    keep.emplace_back(std::move(c0));
    keep.emplace_back(std::move(c1));
    keep.emplace_back(std::move(c2));
    keep.emplace_back(std::move(c3));
    return pdf;
  }
  if (bkgOpt == 2) {
    auto p1 = std::make_unique<RooRealVar>(Form("ch1_%s", name.Data()), "ch1", 0.0, -0.5, 0.5);
    auto p2 = std::make_unique<RooRealVar>(Form("ch2_%s", name.Data()), "ch2", 0.0, -0.5, 0.5);
    auto pdf = std::make_unique<RooChebychev>(name, name, v, RooArgList(*p1, *p2));
    keep.emplace_back(std::move(p1));
    keep.emplace_back(std::move(p2));
    return pdf;
  }
  auto p0 = std::make_unique<RooRealVar>(Form("ep0_%s", name.Data()), "ep0", 0.0, -50.0, 50.0);
  auto p1 = std::make_unique<RooRealVar>(Form("ep1_%s", name.Data()), "ep1", 0.0, -50.0, 50.0);
  auto p2 = std::make_unique<RooRealVar>(Form("ep2_%s", name.Data()), "ep2", 0.0, -50.0, 50.0);
  auto pdf = std::make_unique<RooGenericPdf>(name, name, "exp(@0 + @1*@3 + @2*@3*@3)", RooArgList(*p0, *p1, *p2, v));
  keep.emplace_back(std::move(p0));
  keep.emplace_back(std::move(p1));
  keep.emplace_back(std::move(p2));
  return pdf;
}

// --------------------------- 1D signal builder ---------------------------
static RooAbsPdf* buildSig1D(int signalOpt,
                             RooRealVar& v,
                             const TString& base,
                             RooRealVar& /*aTailShared*/,
                             RooRealVar& /*nTailShared*/,
                             std::vector<std::unique_ptr<RooAbsArg>>& owned,
                             RooRealVar*& outMean,
                             RooAbsReal*& outSigmaEff,
                             RooRealVar*& outFTailLike)
{
  auto mean = std::make_unique<RooRealVar>(Form("mean_%s", base.Data()), "mean", 1.115, 1.113, 1.117);
  auto s1 = std::make_unique<RooRealVar>(Form("s1_%s", base.Data()), "s1", 0.0015, 0.0008, 0.010);
  auto s2 = std::make_unique<RooRealVar>(Form("s2_%s", base.Data()), "s2", 0.0022, 0.0010, 0.015);
  auto g1 = std::make_unique<RooGaussian>(Form("g1_%s", base.Data()), "g1", v, *mean, *s1);
  auto g2 = std::make_unique<RooGaussian>(Form("g2_%s", base.Data()), "g2", v, *mean, *s2);
  auto fG = std::make_unique<RooRealVar>(Form("fG_%s", base.Data()), "fG", 0.7, 0.0, 1.0);
  auto var2G = std::make_unique<RooFormulaVar>(Form("var2G_%s", base.Data()), "@0*@1*@1 + (1-@0)*@2*@2", RooArgList(*fG, *s1, *s2));
  auto sig2G = std::make_unique<RooFormulaVar>(Form("sig2G_%s", base.Data()), "sqrt(@0)", RooArgList(*var2G));

  outSigmaEff  = sig2G.get();
  outFTailLike = fG.get();
  if (signalOpt == 3) signalOpt = 2;

  std::unique_ptr<RooAbsPdf> sigPdf;
  if (signalOpt == 0) {
    sigPdf = std::make_unique<RooAddPdf>(Form("sig_%s", base.Data()), "sig", RooArgList(*g1, *g2), RooArgList(*fG));
  } else if (signalOpt == 1) {
    auto sigma = std::make_unique<RooRealVar>(Form("sigma_%s", base.Data()), "sigma", 0.0020, 0.0010, 0.0200);
    auto aL = std::make_unique<RooRealVar>(Form("aL_%s", base.Data()), "aL", 1.2, 0.3, 6.0);
    auto nL = std::make_unique<RooRealVar>(Form("nL_%s", base.Data()), "nL", 8.0, 1.0, 80.0);
    auto aR = std::make_unique<RooRealVar>(Form("aR_%s", base.Data()), "aR", 1.2, 0.3, 6.0);
    auto nR = std::make_unique<RooRealVar>(Form("nR_%s", base.Data()), "nR", 8.0, 1.0, 80.0);
    sigPdf = std::make_unique<RooCrystalBall>(Form("sig_%s", base.Data()), "sig", v, *mean, *sigma, *aL, *nL, *aR, *nR);
    outSigmaEff = sigma.get();
    outFTailLike = nullptr;
    owned.emplace_back(std::move(sigma));
    owned.emplace_back(std::move(aL));
    owned.emplace_back(std::move(nL));
    owned.emplace_back(std::move(aR));
    owned.emplace_back(std::move(nR));
  } else if (signalOpt == 2) {
    auto s3 = std::make_unique<RooRealVar>(Form("s3_%s", base.Data()), "s3", 0.0035, 0.0010, 0.030);
    auto g3 = std::make_unique<RooGaussian>(Form("g3_%s", base.Data()), "g3", v, *mean, *s3);
    auto f1 = std::make_unique<RooRealVar>(Form("f1_%s", base.Data()), "f1", 0.60, 0.0, 1.0);
    auto f2 = std::make_unique<RooRealVar>(Form("f2_%s", base.Data()), "f2", 0.30, 0.0, 1.0);
    auto g12 = std::make_unique<RooAddPdf>(Form("g12_%s", base.Data()), "g12", RooArgList(*g1, *g2), RooArgList(*f1));
    sigPdf = std::make_unique<RooAddPdf>(Form("sig_%s", base.Data()), "sig", RooArgList(*g12, *g3), RooArgList(*f2));
    auto w1 = std::make_unique<RooFormulaVar>(Form("w1_%s", base.Data()), "@0*@1", RooArgList(*f2, *f1));
    auto w2 = std::make_unique<RooFormulaVar>(Form("w2_%s", base.Data()), "@0*(1-@1)", RooArgList(*f2, *f1));
    auto w3 = std::make_unique<RooFormulaVar>(Form("w3_%s", base.Data()), "(1-@0)", RooArgList(*f2));
    auto var3G = std::make_unique<RooFormulaVar>(Form("var3G_%s", base.Data()), "@0*@3*@3 + @1*@4*@4 + @2*@5*@5", RooArgList(*w1, *w2, *w3, *s1, *s2, *s3));
    auto sig3G = std::make_unique<RooFormulaVar>(Form("sig3G_%s", base.Data()), "sqrt(@0)", RooArgList(*var3G));
    outSigmaEff = sig3G.get();
    outFTailLike = f1.get();
    owned.emplace_back(std::move(s3));
    owned.emplace_back(std::move(g3));
    owned.emplace_back(std::move(f1));
    owned.emplace_back(std::move(f2));
    owned.emplace_back(std::move(g12));
    owned.emplace_back(std::move(w1));
    owned.emplace_back(std::move(w2));
    owned.emplace_back(std::move(w3));
    owned.emplace_back(std::move(var3G));
    owned.emplace_back(std::move(sig3G));
  } else if (signalOpt == 4) {
    auto sigmaG = std::make_unique<RooRealVar>(Form("sigmaG_%s", base.Data()), "sigmaG", 0.0015, 0.0005, 0.0080);
    auto kSigma = std::make_unique<RooRealVar>(Form("kSigma_%s", base.Data()), "kSigma", 1.2, 0.7, 3.0);
    auto sigmaCB = std::make_unique<RooFormulaVar>(Form("sigmaCB_%s", base.Data()), "@0*@1", RooArgList(*kSigma, *sigmaG));
    auto gaus = std::make_unique<RooGaussian>(Form("g_%s", base.Data()), "gaus", v, *mean, *sigmaG);
    auto aL = std::make_unique<RooRealVar>(Form("aL_%s", base.Data()), "aL", 1.5, 0.3, 6.0);
    auto nL = std::make_unique<RooRealVar>(Form("nL_%s", base.Data()), "nL", 6.0, 1.0, 80.0);
    auto aR = std::make_unique<RooRealVar>(Form("aR_%s", base.Data()), "aR", 2.0, 0.3, 6.0);
    auto nR = std::make_unique<RooRealVar>(Form("nR_%s", base.Data()), "nR", 6.0, 1.0, 80.0);
    auto dscb = std::make_unique<RooCrystalBall>(Form("cb_%s", base.Data()), "dscb", v, *mean, *sigmaCB, *aL, *nL, *aR, *nR);
    auto fMix = std::make_unique<RooRealVar>(Form("fGmix_%s", base.Data()), "fGmix", 0.5, 0.0, 1.0);
    sigPdf = std::make_unique<RooAddPdf>(Form("sig_%s", base.Data()), "sig", RooArgList(*gaus, *dscb), RooArgList(*fMix));
    outSigmaEff = sigmaG.get();
    outFTailLike = fMix.get();
    owned.emplace_back(std::move(sigmaG));
    owned.emplace_back(std::move(kSigma));
    owned.emplace_back(std::move(sigmaCB));
    owned.emplace_back(std::move(gaus));
    owned.emplace_back(std::move(aL));
    owned.emplace_back(std::move(nL));
    owned.emplace_back(std::move(aR));
    owned.emplace_back(std::move(nR));
    owned.emplace_back(std::move(dscb));
    owned.emplace_back(std::move(fMix));
  } else if (signalOpt == 6) {
    auto sigmaR = std::make_unique<RooRealVar>(Form("sigmaR_%s", base.Data()), "sigmaR", 0.0016, 0.0006, 0.010);
    auto sigmaL = std::make_unique<RooRealVar>(Form("sigmaL_%s", base.Data()), "sigmaL", 0.0016, 0.0006, 0.010);
    auto tauR = std::make_unique<RooRealVar>(Form("tauR_%s", base.Data()), "tauR", 0.0020, 0.0002, 0.020);
    auto tauL = std::make_unique<RooRealVar>(Form("tauL_%s", base.Data()), "tauL", 0.0020, 0.0002, 0.020);
    auto right = std::make_unique<RooGausExpRightTail>(Form("gxR_%s", base.Data()), "gxR", v, *mean, *sigmaR, *tauR);
    auto left  = std::make_unique<RooGausExpLeftTail >(Form("gxL_%s", base.Data()), "gxL", v, *mean, *sigmaL, *tauL);
    auto fR = std::make_unique<RooRealVar>(Form("fR_%s", base.Data()), "fR", 0.5, 0.0, 1.0);
    sigPdf = std::make_unique<RooAddPdf>(Form("sig_%s", base.Data()), "sig", RooArgList(*right, *left), RooArgList(*fR));
    outSigmaEff = sigmaR.get();
    outFTailLike = fR.get();
    owned.emplace_back(std::move(sigmaR));
    owned.emplace_back(std::move(sigmaL));
    owned.emplace_back(std::move(tauR));
    owned.emplace_back(std::move(tauL));
    owned.emplace_back(std::move(right));
    owned.emplace_back(std::move(left));
    owned.emplace_back(std::move(fR));
  } else {
    sigPdf = std::make_unique<RooAddPdf>(Form("sig_%s", base.Data()), "sig", RooArgList(*g1, *g2), RooArgList(*fG));
  }

  outMean = mean.get();
  owned.emplace_back(std::move(mean));
  owned.emplace_back(std::move(s1));
  owned.emplace_back(std::move(s2));
  owned.emplace_back(std::move(g1));
  owned.emplace_back(std::move(g2));
  owned.emplace_back(std::move(fG));
  owned.emplace_back(std::move(var2G));
  owned.emplace_back(std::move(sig2G));
  RooAbsPdf* raw = sigPdf.release();
  owned.emplace_back(std::unique_ptr<RooAbsArg>(raw));
  return raw;
}

static void BuildWindowNoRound(double mean, double sigma, double nSigma, double& low, double& high)
{
  low  = mean - nSigma * sigma;
  high = mean + nSigma * sigma;
}

static void BuildSignalWindowAligned(const TH2D* h2,
                                     double meanX, double sigmaX,
                                     double meanY, double sigmaY,
                                     double nSigmaWin,
                                     double& xlow, double& xhi,
                                     double& ylow, double& yhi,
                                     double& xlow_bin, double& xhi_bin,
                                     double& ylow_bin, double& yhi_bin)
{
  BuildWindowNoRound(meanX, sigmaX, nSigmaWin, xlow, xhi);
  BuildWindowNoRound(meanY, sigmaY, nSigmaWin, ylow, yhi);

  auto axHX = h2->GetXaxis();
  auto axHY = h2->GetYaxis();
  const int bx1_sig = std::max(1, std::min(axHX->FindBin(xlow), axHX->GetNbins()));
  const int bx2_sig = std::max(1, std::min(axHX->FindBin(xhi - 1e-12), axHX->GetNbins()));
  const int by1_sig = std::max(1, std::min(axHY->FindBin(ylow), axHY->GetNbins()));
  const int by2_sig = std::max(1, std::min(axHY->FindBin(yhi - 1e-12), axHY->GetNbins()));

  xlow_bin = axHX->GetBinLowEdge(std::min(bx1_sig, bx2_sig));
  xhi_bin  = axHX->GetBinUpEdge (std::max(bx1_sig, bx2_sig));
  ylow_bin = axHY->GetBinLowEdge(std::min(by1_sig, by2_sig));
  yhi_bin  = axHY->GetBinUpEdge (std::max(by1_sig, by2_sig));
}



//signalOpt = 0: Double Gaussian
//signalOpt = 1: Double-sided CrystalBall(DCSB)
//signalOpt = 2: Triple Gaussian
//signalOpt = 3: Triple Gaussian
//signalOpt = 4: One Gaussian + Double-sided CrystalBall(DSCB)
//signalOpt = 6: GausExpLeft/RightTail

//bkgOpt = 0: 3rd order Bernstein
//bkgOpt = 1: 4th order Bernstein
//bkgOpt = 2: 2nd order Chebychev
//bkgOpt = other: Exponential of 2nd order polynomial

void Analysis_Rbin(TString inFile = "AnalysisResults.root",
                   TString outDir = "fit_Rbins_out",
                   TString taskPath = "lambdaspincorrderived",
                   TString sparseSEName = "hSparseLambdaLambda",
                   TString sparseMEName = "hSparseLambdaLambdaMixed",
                   bool interactive = false,
                   bool savePng = true,
                   bool doME = true,
                   bool meUseSEWindow = true,
                   int signalOpt = 0,//signal function
                   int bkgOpt = 1,//background function
                   bool propagateFracErr = true,
                   double sbLeftLow = 1.104,//side band range low edge for left side band
                   double sbLeftHigh = 1.108,
                   double sbRightLow = 1.122,
                   double sbRightHigh = 1.126)
{
  gROOT->SetBatch(!interactive);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gSystem->Load("libRooFitCore");
  gSystem->Load("libRooFit");

  const double VAR_MIN = 1.1;
  const double VAR_MAX = 1.13;
  const double FIT_MIN = 1.1;
  const double FIT_MAX = 1.13;
  const double nSigmaWin = 2.0;

  const bool useManualWindowReference = false;
  const double manualMeanX  = 1.115;
  const double manualMeanY  = 1.115;
  const double manualSigmaX = 0.00175;
  const double manualSigmaY = 0.00175;

  const int AX_M1  = 0;
  const int AX_M2  = 1;
  const int AX_COS = 2;
  const int AX_R   = 3;
  const int AX_Rap   = 4;
  const int AX_Phi   = 5;

  double rapidityMin = 0.0;
  double rapidityMax = 0.5;
  //const double REdges[] = {0.0, TMath::Pi()/8.0, TMath::Pi()/4.0, TMath::Pi()/2.0, 3.0*TMath::Pi()/4.0, TMath::Pi()};
  //const double REdges[] = {0.0, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.3, 3.1};
  //const double REdges[] = {0.0, 3.1};
  const double REdges[] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};
  const int nRbins = sizeof(REdges)/sizeof(REdges[0]) - 1;// 6 R bins

  const int nCosBins = 10;
  const double cosEdges[nCosBins + 1] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};

  const int colTot = kBlack;
  const int colSS  = kRed + 1;
  const int colSB  = kBlue + 1;
  const int colBS  = kGreen + 2;
  const int colBB  = kMagenta + 1;

  if (!(sbLeftLow < sbLeftHigh && sbRightLow < sbRightHigh &&
        sbLeftLow >= VAR_MIN && sbRightHigh <= VAR_MAX)) {
    std::cerr << "ERROR: invalid fixed sidebands.\n";
    return;
  }

  //open AnalysisResults.root
  TFile* fin = TFile::Open(inFile);
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inFile << "\n";
    if (fin) { fin->Close(); delete fin; }
    return;
  }

  auto getSparse = [&](const TString& name)->THnSparse*{
    TString fullPath = Form("%s/%s", taskPath.Data(), name.Data());
    return dynamic_cast<THnSparse*>(fin->Get(fullPath));
  };

  THnSparse* hsSE = getSparse(sparseSEName);
  if (!hsSE) {
    std::cerr << "ERROR: cannot find THnSparse at: " << Form("%s/%s", taskPath.Data(), sparseSEName.Data()) << "\n";
    fin->Close();
    delete fin;
    return;
  }

  THnSparse* hsME = nullptr;
  if (doME) {
    hsME = getSparse(sparseMEName);
    if (!hsME) {
      std::cerr << "WARNING: ME sparse not found at: " << Form("%s/%s", taskPath.Data(), sparseMEName.Data()) << "\n";
      doME = false;
    }
  }

  TAxis* axM1  = hsSE->GetAxis(AX_M1);
  TAxis* axM2  = hsSE->GetAxis(AX_M2);
  TAxis* axCos = hsSE->GetAxis(AX_COS);
  TAxis* axR   = hsSE->GetAxis(AX_R);
  // TAxis* axRap   = hsSE->GetAxis(AX_Rap);
  // TAxis* axPhi   = hsSE->GetAxis(AX_Phi);


  TAxis* axM1m  = hsME->GetAxis(AX_M1);
  TAxis* axM2m  = hsME->GetAxis(AX_M2);
  TAxis* axCosm = hsME->GetAxis(AX_COS);
  TAxis* axRm   = hsME->GetAxis(AX_R);
  // TAxis* axPhim   = hsME->GetAxis(AX_Phi);
  // TAxis* axRapm   = hsME->GetAxis(AX_Rap);

  
  axM1->SetRange(axM1->FindBin(VAR_MIN + 1e-6), axM1->FindBin(VAR_MAX - 1e-6));
  axM2->SetRange(axM2->FindBin(VAR_MIN + 1e-6), axM2->FindBin(VAR_MAX - 1e-6));
  //axRap->SetRange(axRap->FindBin(rapidityMin + 1e-6), axRap->FindBin(rapidityMax - 1e-6));

  axM1m->SetRange(axM1m->FindBin(VAR_MIN + 1e-6), axM1m->FindBin(VAR_MAX - 1e-6));
  axM2m->SetRange(axM2m->FindBin(VAR_MIN + 1e-6), axM2m->FindBin(VAR_MAX - 1e-6));
  //axRapm ->SetRange(axRapm ->FindBin(rapidityMin + 1e-6),     axRapm ->FindBin(rapidityMax - 1e-6));

  gSystem->mkdir(outDir.Data(), kTRUE);

  TH1D* hFsAvg_vs_R = new TH1D("hFsAvg_vs_R", "Average f_s;#DeltaR;f_s", nRbins, REdges);
    

  for (int iR = 0; iR < nRbins; ++iR) {
    const double rLo = REdges[iR];
    const double rHi = REdges[iR+1];

    //axPhi->SetRange(axPhi->FindBin(rLo + 1e-6), axPhi->FindBin(rHi - 1e-6));
    //axPhim ->SetRange(axPhim ->FindBin(rLo + 1e-6),     axPhim ->FindBin(rHi - 1e-6));

    axR->SetRange(axR->FindBin(rLo + 1e-6), axR->FindBin(rHi - 1e-6));
    axRm ->SetRange(axRm ->FindBin(rLo + 1e-6),     axRm ->FindBin(rHi - 1e-6));
      
    TString outRoot = Form("%s/fit_Rbin%d.root", outDir.Data(), iR);
    TFile* fout = TFile::Open(outRoot, "RECREATE");
    if (!fout || fout->IsZombie()) {
      std::cerr << "ERROR: cannot create output file: " << outRoot << "\n";
      if (fout) { fout->Close(); delete fout; }
      continue;
    }
    fout->cd();

    TH1D* hYieldSS          = new TH1D("hYieldSS_vsCos","SS yield in sigWin;cos#theta*;Yield", nCosBins, cosEdges);
    TH1D* hYieldSS_Model    = new TH1D("hYieldSS_Model_vsCos","SS model yield in sigWin;cos#theta*;Yield", nCosBins, cosEdges);
    TH1D* hYieldSBsum_Model = new TH1D("hYieldSBsum_Model_vsCos","(SB+BS) model yield in sigWin;cos#theta*;Yield", nCosBins, cosEdges);
    TH1D* hYieldBB_Model    = new TH1D("hYieldBB_Model_vsCos","BB model yield in sigWin;cos#theta*;Yield", nCosBins, cosEdges);
    TH1D* hEDM       = new TH1D("hEDM_vsCos", "EDM;cos#theta*;EDM", nCosBins, cosEdges);
    TH1D* hStatus    = new TH1D("hFitStatus_vsCos", "fit status;cos#theta*;status", nCosBins, cosEdges);
    TH1D* hCovQual   = new TH1D("hCovQual_vsCos", "covQual;cos#theta*;covQual", nCosBins, cosEdges);
    TH1D* hMeanX     = new TH1D("hMeanX_vsCos", "#mu_{x};cos#theta*;#mu_{x}", nCosBins, cosEdges);
    TH1D* hMeanY     = new TH1D("hMeanY_vsCos", "#mu_{y};cos#theta*;#mu_{y}", nCosBins, cosEdges);
    TH1D* hSigmaEffX = new TH1D("hSigmaEffX_vsCos", "#sigma^{eff}_{x};cos#theta*;#sigma^{eff}_{x}", nCosBins, cosEdges);
    TH1D* hSigmaEffY = new TH1D("hSigmaEffY_vsCos", "#sigma^{eff}_{y};cos#theta*;#sigma^{eff}_{y}", nCosBins, cosEdges);
    TH1D* hNSS_vsCos = new TH1D("hNSS_vsCos","N_{SS};cos#theta*;N_{SS}", nCosBins, cosEdges);
    TH1D* hNSSwin_vsCos = new TH1D("hNSSwin_vsCos","N_{SS}^{win};cos#theta*;N_{SS}^{win}", nCosBins, cosEdges);
    TH1D* hFracSSwin_vsCos = new TH1D("hFracSSwin_vsCos","f_{SS}^{win};cos#theta*;f_{SS}^{win}", nCosBins, cosEdges);
    TH1D* hDataWin_vsCos = new TH1D("hDataWin_vsCos","Data integral in SS window (SE);cos#theta*;N_{data}^{win}", nCosBins, cosEdges);
    TH1D* hSigWin_BCD_vsCos   = new TH1D("hSigWin_BCD_vsCos","S^{win} from A-(B+C-D);cos#theta*;S^{win}", nCosBins, cosEdges);
    TH1D* hBkgWin_BCD_vsCos   = new TH1D("hBkgWin_BCD_vsCos","B^{win} from (B+C-D);cos#theta*;B^{win}", nCosBins, cosEdges);
    
    //Add the histogram for the background B,C and D separately in SE.
    TH1D* hBkgWin_B_vsCos   = new TH1D("hBkgWin_B_vsCos","B^{win} from B;cos#theta*;B^{win}", nCosBins, cosEdges);
    TH1D* hBkgWin_C_vsCos   = new TH1D("hBkgWin_C_vsCos","B^{win} from C;cos#theta*;B^{win}", nCosBins, cosEdges);
    TH1D* hBkgWin_D_vsCos   = new TH1D("hBkgWin_D_vsCos","B^{win} from D;cos#theta*;B^{win}", nCosBins, cosEdges);

    
    TH1D* hFracSig_BCD_vsCos  = new TH1D("hFracSig_BCD_vsCos","f_{sig}^{win} from B+C-D;cos#theta*;f_{sig}^{win}", nCosBins, cosEdges);
    TH1D* hX_SE_vsCos         = new TH1D("hX_SE_vsCos","X^{SE}=B+C-D;cos#theta*;X^{SE}", nCosBins, cosEdges);
    TH1D* hFTailX_vsCos = new TH1D("hFTailX_vsCos","Tail fraction X;cos#theta*;f_{tail}^{x}", nCosBins, cosEdges);
    TH1D* hFTailY_vsCos = new TH1D("hFTailY_vsCos","Tail fraction Y;cos#theta*;f_{tail}^{y}", nCosBins, cosEdges);
    TH1D* hWinXlow_SE = new TH1D("hWinXlow_SE","xlow(SE);cos#theta*;xlow", nCosBins, cosEdges);
    TH1D* hWinXhi_SE  = new TH1D("hWinXhi_SE","xhi(SE);cos#theta*;xhi", nCosBins, cosEdges);
    TH1D* hWinYlow_SE = new TH1D("hWinYlow_SE","ylow(SE);cos#theta*;ylow", nCosBins, cosEdges);
    TH1D* hWinYhi_SE  = new TH1D("hWinYhi_SE","yhi(SE);cos#theta*;yhi", nCosBins, cosEdges);

    for (auto* h : {hYieldSS,hYieldSS_Model,hYieldSBsum_Model,hYieldBB_Model,
                    hEDM,hStatus,hCovQual,hMeanX,hMeanY,hSigmaEffX,hSigmaEffY,
                    hNSS_vsCos,hNSSwin_vsCos,hFracSSwin_vsCos,hDataWin_vsCos,
                    hSigWin_BCD_vsCos,hBkgWin_BCD_vsCos,hBkgWin_B_vsCos,hBkgWin_C_vsCos,hBkgWin_D_vsCos,hFracSig_BCD_vsCos,hX_SE_vsCos,
                    hFTailX_vsCos,hFTailY_vsCos,hWinXlow_SE,hWinXhi_SE,hWinYlow_SE,hWinYhi_SE}) h->Sumw2();

    TH1D *hYieldSS_ME=nullptr, *hNSS_ME=nullptr, *hFracSSwin_ME=nullptr, *hDataWin_vsCos_ME=nullptr;
    TH1D *hSigWin_BCD_vsCos_ME=nullptr, *hBkgWin_BCD_vsCos_ME=nullptr, *hFracSig_BCD_vsCos_ME=nullptr;
    TH1D *hYieldSS_Model_ME=nullptr, *hYieldSBsum_Model_ME=nullptr, *hYieldBB_Model_ME=nullptr;
    TH1D *hMeanX_ME=nullptr, *hMeanY_ME=nullptr, *hSigmaEffX_ME=nullptr, *hSigmaEffY_ME=nullptr;
    TH1D *hWinXlow_used_ME=nullptr, *hWinXhi_used_ME=nullptr, *hWinYlow_used_ME=nullptr, *hWinYhi_used_ME=nullptr;
    //Define the histograms for the background B,C,D separately in ME
    TH1D *hBkgWin_B_vsCos_ME=nullptr, *hBkgWin_C_vsCos_ME=nullptr, *hBkgWin_D_vsCos_ME=nullptr;

    if (doME && hsME) {
      hYieldSS_ME          = new TH1D("hYieldSS_vsCos_ME","SS yield in sigWin (ME);cos#theta*;Yield", nCosBins, cosEdges);
      hYieldSS_Model_ME    = new TH1D("hYieldSS_model_vsCos_ME","SS model yield in sigWin (ME);cos#theta*;Yield", nCosBins, cosEdges);
      hYieldSBsum_Model_ME = new TH1D("hYieldSBsum_Model_vsCos_ME","(SB+BS) model yield in sigWin (ME);cos#theta*;Yield", nCosBins, cosEdges);
      hYieldBB_Model_ME    = new TH1D("hYieldBB_Model_vsCos_ME","BB model yield in sigWin (ME);cos#theta*;Yield", nCosBins, cosEdges);
      hNSS_ME     = new TH1D("hNSS_vsCos_ME","N_{SS} (ME);cos#theta*;N_{SS}", nCosBins, cosEdges);
      hFracSSwin_ME = new TH1D("hFracSSwin_vsCos_ME","f_{SS}^{win} (ME);cos#theta*;f_{SS}^{win}", nCosBins, cosEdges);
      hDataWin_vsCos_ME = new TH1D("hDataWin_vsCos_ME","Data integral in SS window (ME);cos#theta*;N_{data}^{win}", nCosBins, cosEdges);
      hSigWin_BCD_vsCos_ME   = new TH1D("hSigWin_BCD_vsCos_ME","S^{win} (ME) from A-(B+C-D);cos#theta*;S^{win}", nCosBins, cosEdges);
      hBkgWin_BCD_vsCos_ME   = new TH1D("hBkgWin_BCD_vsCos_ME","B^{win} (ME) from (B+C-D);cos#theta*;B^{win}", nCosBins, cosEdges);
      
      //Add the histogram for the background B,C and D separately in ME
      hBkgWin_B_vsCos_ME   = new TH1D("hBkgWin_B_vsCos_ME","B^{win} from B;cos#theta*;B^{win}", nCosBins, cosEdges);
      hBkgWin_C_vsCos_ME   = new TH1D("hBkgWin_C_vsCos_ME","B^{win} from C;cos#theta*;B^{win}", nCosBins, cosEdges);
      hBkgWin_D_vsCos_ME   = new TH1D("hBkgWin_D_vsCos_ME","B^{win} from D;cos#theta*;B^{win}", nCosBins, cosEdges);
      
      
      
      hFracSig_BCD_vsCos_ME  = new TH1D("hFracSig_BCD_vsCos_ME","f_{sig}^{win} (ME) from B+C-D;cos#theta*;f_{sig}^{win}", nCosBins, cosEdges);
      hMeanX_ME = new TH1D("hMeanX_vsCos_ME","#mu_{x} (ME);cos#theta*;#mu_{x}", nCosBins, cosEdges);
      hMeanY_ME = new TH1D("hMeanY_vsCos_ME","#mu_{y} (ME);cos#theta*;#mu_{y}", nCosBins, cosEdges);
      hSigmaEffX_ME = new TH1D("hSigmaEffX_vsCos_ME","#sigma^{eff}_{x} (ME);cos#theta*;#sigma^{eff}_{x}", nCosBins, cosEdges);
      hSigmaEffY_ME = new TH1D("hSigmaEffY_vsCos_ME","#sigma^{eff}_{y} (ME);cos#theta*;#sigma^{eff}_{y}", nCosBins, cosEdges);
      hWinXlow_used_ME = new TH1D("hWinXlow_used_ME","xlow used for ME;cos#theta*;xlow", nCosBins, cosEdges);
      hWinXhi_used_ME  = new TH1D("hWinXhi_used_ME","xhi used for ME;cos#theta*;xhi", nCosBins, cosEdges);
      hWinYlow_used_ME = new TH1D("hWinYlow_used_ME","ylow used for ME;cos#theta*;ylow", nCosBins, cosEdges);
      hWinYhi_used_ME  = new TH1D("hWinYhi_used_ME","yhi used for ME;cos#theta*;yhi", nCosBins, cosEdges);
      for (auto* h : {hYieldSS_ME, hYieldSS_Model_ME, hYieldSBsum_Model_ME, hYieldBB_Model_ME,
                      hNSS_ME, hFracSSwin_ME, hDataWin_vsCos_ME,
                      hSigWin_BCD_vsCos_ME, hBkgWin_BCD_vsCos_ME, hFracSig_BCD_vsCos_ME,
                      hMeanX_ME, hMeanY_ME, hSigmaEffX_ME, hSigmaEffY_ME,
                      hWinXlow_used_ME, hWinXhi_used_ME, hWinYlow_used_ME, hWinYhi_used_ME}) h->Sumw2();
    }
    for (int k = 0; k < nCosBins; ++k) {

      const double cLo = cosEdges[k];
      const double cHi = cosEdges[k + 1];

      const int bLo = findBinForEdgeLow(axCos, cLo);
      const int bHi = findBinForEdgeHigh(axCos, cHi);
      axCos->SetRange(bLo, bHi);
      axCos->SetBit(TAxis::kAxisRange);

      TH2D* h2 = dynamic_cast<TH2D*>(hsSE->Projection(AX_M2, AX_M1, "E"));
      if (!h2) { axCos->SetRange(0,0); continue; }
      h2->SetDirectory(nullptr);
      h2->Sumw2();

      const double sumW = h2->Integral();
      if (sumW <= 0) { delete h2; axCos->SetRange(0,0); continue; }

      const TString tag = Form("_R%d_c%02d", iR, k+1);

      RooRealVar x(Form("x%s", tag.Data()), "m_{#Lambda1}", VAR_MIN, VAR_MAX);
      RooRealVar y(Form("y%s", tag.Data()), "m_{#Lambda2}", VAR_MIN, VAR_MAX);
      x.setRange("fit", FIT_MIN, FIT_MAX);
      y.setRange("fit", FIT_MIN, FIT_MAX);

      RooDataHist data2D(Form("data2D%s", tag.Data()), "data2D", RooArgList(x, y), h2);

      RooRealVar aTail(Form("aTail%s", tag.Data()), "aTail", 3.0, 2.5, 8.0);
      RooRealVar nTail(Form("nTail%s", tag.Data()), "nTail", 10.0, 2.0, 20.0);

      std::vector<std::unique_ptr<RooAbsArg>> ownedX, ownedY;
      RooRealVar *meanX=nullptr, *meanY=nullptr;
      RooAbsReal *sigEffX=nullptr, *sigEffY=nullptr;
      RooRealVar *fTailX=nullptr, *fTailY=nullptr;

      RooAbsPdf* sig1_ptr = buildSig1D(signalOpt, x, Form("X%s",tag.Data()), aTail, nTail, ownedX, meanX, sigEffX, fTailX);
      RooAbsPdf* sig2_ptr = buildSig1D(signalOpt, y, Form("Y%s",tag.Data()), aTail, nTail, ownedY, meanY, sigEffY, fTailY);

      std::vector<std::unique_ptr<RooAbsArg>> keepBx, keepBy;
      auto bkg1_ptr = buildBkg1D(bkgOpt, x, Form("bkg1%s",tag.Data()), keepBx);
      auto bkg2_ptr = buildBkg1D(bkgOpt, y, Form("bkg2%s",tag.Data()), keepBy);

      RooProdPdf pdfSS(Form("pdfSS%s", tag.Data()), "SS", RooArgList(*sig1_ptr, *sig2_ptr));
      RooProdPdf pdfSB(Form("pdfSB%s", tag.Data()), "SB", RooArgList(*sig1_ptr, *bkg2_ptr));
      RooProdPdf pdfBS(Form("pdfBS%s", tag.Data()), "BS", RooArgList(*bkg1_ptr, *sig2_ptr));
      RooProdPdf pdfBB(Form("pdfBB%s", tag.Data()), "BB", RooArgList(*bkg1_ptr, *bkg2_ptr));

      const double Nguess = std::max(1.0, sumW);
      RooRealVar N_SS(Form("N_SS%s", tag.Data()), "N_SS", 0.80*Nguess, 0.0, 20.0*Nguess);
      RooRealVar N_SB(Form("N_SB%s", tag.Data()), "N_SB", 0.07*Nguess, 0.0, 20.0*Nguess);
      RooRealVar N_BS(Form("N_BS%s", tag.Data()), "N_BS", 0.07*Nguess, 0.0, 20.0*Nguess);
      RooRealVar N_BB(Form("N_BB%s", tag.Data()), "N_BB", 0.06*Nguess, 0.0, 20.0*Nguess);

      RooAddPdf model(Form("model%s", tag.Data()), "model",
                      RooArgList(pdfSS,pdfSB,pdfBS,pdfBB),
                      RooArgList(N_SS,N_SB,N_BS,N_BB));

      std::unique_ptr<RooAbsReal> nll(model.createNLL(data2D, Extended(true), Range("fit"), Offset(true)));
      if (!nll) { delete h2; axCos->SetRange(0,0); continue; }

      RooMinimizer minim(*nll);
      minim.setMinimizerType("Minuit2");
      minim.setStrategy(2);
      minim.setPrintLevel(interactive ? 1 : -1);
      minim.setMaxFunctionCalls(200000);
      minim.setMaxIterations(200000);
      minim.setEps(1e-4);
      minim.setOffsetting(true);
      minim.optimizeConst(2);
      minim.minimize("Minuit2","simplex");
      minim.minimize("Minuit2","migrad");
      minim.hesse();

      std::unique_ptr<RooFitResult> res(minim.save());
      if (!res) { delete h2; axCos->SetRange(0,0); continue; }

      const double mx = meanX ? meanX->getVal() : 1.115;
      const double my = meanY ? meanY->getVal() : 1.115;
      const double sx = sigEffX ? sigEffX->getVal() : 0.002;
      const double sy = sigEffY ? sigEffY->getVal() : 0.002;

      const double exSig = GetPropErrSafe(sigEffX, res.get());
      const double eySig = GetPropErrSafe(sigEffY, res.get());

      double xlow = 0.0, xhi = 0.0;
      double ylow = 0.0, yhi = 0.0;
      double xlow_bin = 0.0, xhi_bin = 0.0;
      double ylow_bin = 0.0, yhi_bin = 0.0;

      const double winMeanX  = useManualWindowReference ? manualMeanX  : mx;
      const double winMeanY  = useManualWindowReference ? manualMeanY  : my;
      const double winSigmaX = useManualWindowReference ? manualSigmaX : sx;
      const double winSigmaY = useManualWindowReference ? manualSigmaY : sy;

      BuildSignalWindowAligned(h2,
                               winMeanX, winSigmaX,
                               winMeanY, winSigmaY,
                               nSigmaWin,
                               xlow, xhi, ylow, yhi,
                               xlow_bin, xhi_bin, ylow_bin, yhi_bin);

      x.setRange("sigWin", xlow_bin, xhi_bin);
      y.setRange("sigWin", ylow_bin, yhi_bin);

      hWinXlow_SE->SetBinContent(k+1, xlow_bin);
      hWinXhi_SE ->SetBinContent(k+1, xhi_bin);
      hWinYlow_SE->SetBinContent(k+1, ylow_bin);
      hWinYhi_SE ->SetBinContent(k+1, yhi_bin);

      RooArgSet normSet(x,y);

      hNSS_vsCos->SetBinContent(k+1, N_SS.getVal());
      hNSS_vsCos->SetBinError  (k+1, N_SS.getError());

      auto intSS_win = std::unique_ptr<RooAbsReal>(pdfSS.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
      auto intSB_win = std::unique_ptr<RooAbsReal>(pdfSB.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
      auto intBS_win = std::unique_ptr<RooAbsReal>(pdfBS.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
      auto intBB_win = std::unique_ptr<RooAbsReal>(pdfBB.createIntegral(normSet, NormSet(normSet), Range("sigWin")));

      const double YSS_win = N_SS.getVal() * (intSS_win ? intSS_win->getVal() : 0.0);
      double eYSS_win = 0.0;
      if (intSS_win && N_SS.getVal() > 0 && intSS_win->getVal() > 0) {
        eYSS_win = YSS_win * N_SS.getError() / std::max(1e-12, N_SS.getVal());
      }

      hYieldSS_Model->SetBinContent(k+1, YSS_win);
      hYieldSS_Model->SetBinError  (k+1, eYSS_win);
      hNSSwin_vsCos->SetBinContent(k+1, YSS_win);
      hNSSwin_vsCos->SetBinError  (k+1, eYSS_win);

      const double YSB_win = N_SB.getVal() * (intSB_win ? intSB_win->getVal() : 0.0);
      const double YBS_win = N_BS.getVal() * (intBS_win ? intBS_win->getVal() : 0.0);
      const double YBB_win = N_BB.getVal() * (intBB_win ? intBB_win->getVal() : 0.0);

      double eYSB_win = 0.0;
      double eYBS_win = 0.0;
      double eYBB_win = 0.0;

      if (intSB_win && N_SB.getVal() > 0 && intSB_win->getVal() > 0) eYSB_win = YSB_win * N_SB.getError() / std::max(1e-12, N_SB.getVal());
      if (intBS_win && N_BS.getVal() > 0 && intBS_win->getVal() > 0) eYBS_win = YBS_win * N_BS.getError() / std::max(1e-12, N_BS.getVal());
      if (intBB_win && N_BB.getVal() > 0 && intBB_win->getVal() > 0) eYBB_win = YBB_win * N_BB.getError() / std::max(1e-12, N_BB.getVal());

      hYieldSBsum_Model->SetBinContent(k+1, YSB_win + YBS_win);
      hYieldSBsum_Model->SetBinError  (k+1, std::sqrt(eYSB_win*eYSB_win + eYBS_win*eYBS_win));
      hYieldBB_Model->SetBinContent(k+1, YBB_win);
      hYieldBB_Model->SetBinError  (k+1, eYBB_win);

      RooFormulaVar fracSSwin(Form("fracSSwin%s",tag.Data()),
                              "(@0*@1)/(@0*@1 + @2*@3 + @4*@5 + @6*@7)",
                              RooArgList(N_SS, *intSS_win, N_SB, *intSB_win, N_BS, *intBS_win, N_BB, *intBB_win));

      const double fss  = fracSSwin.getVal();
      const double efss = fracSSwin.getPropagatedError(*res);

      hFracSSwin_vsCos->SetBinContent(k+1, fss);
      hFracSSwin_vsCos->SetBinError  (k+1, efss);

      double Ndata_win = 0.0, eNdata_win = 0.0;
      IntegralAndErrorH2Window(h2, xlow_bin, xhi_bin, ylow_bin, yhi_bin, Ndata_win, eNdata_win);
      hDataWin_vsCos->SetBinContent(k+1, Ndata_win);
      hDataWin_vsCos->SetBinError  (k+1, eNdata_win);

      {
        const double xSB1_lo = sbLeftLow;
        const double xSB1_hi = sbLeftHigh;
        const double xSB2_lo = sbRightLow;
        const double xSB2_hi = sbRightHigh;

        const double ySB1_lo = sbLeftLow;
        const double ySB1_hi = sbLeftHigh;
        const double ySB2_lo = sbRightLow;
        const double ySB2_hi = sbRightHigh;

        const bool xSBinside = (xSB1_lo >= VAR_MIN) && (xSB2_hi <= VAR_MAX);
        const bool ySBinside = (ySB1_lo >= VAR_MIN) && (ySB2_hi <= VAR_MAX);

        if (xSBinside && ySBinside) {
          const double sigWidthX = xhi_bin - xlow_bin;
          const double sigWidthY = yhi_bin - ylow_bin;

          const double totSbWidthX = (xSB1_hi - xSB1_lo) + (xSB2_hi - xSB2_lo);
          const double totSbWidthY = (ySB1_hi - ySB1_lo) + (ySB2_hi - ySB2_lo);

          const double rX = sigWidthX / std::max(1e-12, totSbWidthX);
          const double rY = sigWidthY / std::max(1e-12, totSbWidthY);

          ValErr Ase{Ndata_win, eNdata_win};

          ValErr Braw = AddVE(
            IntegralVE(h2, xlow_bin, xhi_bin, ySB1_lo, ySB1_hi),
            IntegralVE(h2, xlow_bin, xhi_bin, ySB2_lo, ySB2_hi)
          );
          ValErr Craw = AddVE(
            IntegralVE(h2, xSB1_lo, xSB1_hi, ylow_bin, yhi_bin),
            IntegralVE(h2, xSB2_lo, xSB2_hi, ylow_bin, yhi_bin)
          );
          ValErr Draw = AddVE(
            AddVE(IntegralVE(h2, xSB1_lo, xSB1_hi, ySB1_lo, ySB1_hi),
                  IntegralVE(h2, xSB1_lo, xSB1_hi, ySB2_lo, ySB2_hi)),
            AddVE(IntegralVE(h2, xSB2_lo, xSB2_hi, ySB1_lo, ySB1_hi),
                  IntegralVE(h2, xSB2_lo, xSB2_hi, ySB2_lo, ySB2_hi))
          );

          ValErr B = ScaleVE(Braw, rY);
          ValErr C = ScaleVE(Craw, rX);
          ValErr D = ScaleVE(Draw, rX*rY);

          ValErr Xse = SubVE(AddVE(B, C), D);
          ValErr Bkg = Xse;
          ValErr Sig = MakeSigFromAandBkg(Ase, Bkg);
          ValErr FracSig = MakeFracSigFromAandBkg(Ase, Bkg);

          hX_SE_vsCos->SetBinContent(k+1, Xse.v);
          hX_SE_vsCos->SetBinError  (k+1, Xse.e);

          hBkgWin_BCD_vsCos->SetBinContent(k+1, Bkg.v);
          hBkgWin_BCD_vsCos->SetBinError  (k+1, Bkg.e);
          //For the mixed events, the B,C,D should considered sepreately
          //sepreately store the B and  C and D for the ME in the ME histograms
          /*
          hBkgWin_B_vsCos->SetBinContent(k+1, B.v);
          hBkgWin_B_vsCos->SetBinError  (k+1, B.e);
          hBkgWin_C_vsCos->SetBinContent(k+1, C.v);
          hBkgWin_C_vsCos->SetBinError  (k+1, C.e);
          hBkgWin_D_vsCos->SetBinContent(k+1, D.v);
          hBkgWin_D_vsCos->SetBinError  (k+1, D.e);
          */
          hBkgWin_B_vsCos->SetBinContent(k+1, Braw.v);
          hBkgWin_B_vsCos->SetBinError  (k+1, Braw.e);
          hBkgWin_C_vsCos->SetBinContent(k+1, Craw.v);
          hBkgWin_C_vsCos->SetBinError  (k+1, Craw.e);
          hBkgWin_D_vsCos->SetBinContent(k+1, Draw.v);
          hBkgWin_D_vsCos->SetBinError  (k+1, Draw.e);
            

          //

          hSigWin_BCD_vsCos->SetBinContent(k+1, Sig.v);
          hSigWin_BCD_vsCos->SetBinError  (k+1, Sig.e);

          hFracSig_BCD_vsCos->SetBinContent(k+1, FracSig.v);
          hFracSig_BCD_vsCos->SetBinError  (k+1, FracSig.e);

          if (savePng) {
            TH2D* h2ABCD = (TH2D*)h2->Clone(Form("h2ABCD_SE_R%d_c%02d", iR, k+1));
            h2ABCD->SetDirectory(nullptr);

            TCanvas* cABCD = MakeABCDCanvas(
              h2ABCD,
              Form("cABCD_SE_R%d_c%02d", iR, k+1),
              Form("SE ABCD regions: Rbin %d, cos bin %d", iR, k+1),
              xlow_bin, xhi_bin, ylow_bin, yhi_bin,
              xSB1_lo, xSB1_hi, xSB2_lo, xSB2_hi,
              ySB1_lo, ySB1_hi, ySB2_lo, ySB2_hi,
              true
            );

            if (cABCD) {
              TString png2D = Form("%s/ABCD_SE_Rbin%d_cos%02d.png", outDir.Data(), iR, k+1);
              cABCD->SaveAs(png2D);
              fout->cd();
              h2ABCD->Write();
              cABCD->Write();
              delete cABCD;
            }
            delete h2ABCD;
          }
        }
      }

      const double Ysig = fss * Ndata_win;
      double eYsig = std::abs(Ndata_win) * efss;
      if (propagateFracErr) {
        eYsig = std::sqrt((Ndata_win*efss)*(Ndata_win*efss) + (fss*eNdata_win)*(fss*eNdata_win));
      }

      hYieldSS->SetBinContent(k+1, Ysig);
      hYieldSS->SetBinError  (k+1, eYsig);

      hEDM->SetBinContent(k+1, res->edm());
      hStatus->SetBinContent(k+1, res->status());
      hCovQual->SetBinContent(k+1, res->covQual());
      hMeanX->SetBinContent(k+1, mx); if (meanX) hMeanX->SetBinError(k+1, meanX->getError());
      hMeanY->SetBinContent(k+1, my); if (meanY) hMeanY->SetBinError(k+1, meanY->getError());
      hSigmaEffX->SetBinContent(k+1, sx); hSigmaEffX->SetBinError(k+1, exSig);
      hSigmaEffY->SetBinContent(k+1, sy); hSigmaEffY->SetBinError(k+1, eySig);

      hFTailX_vsCos->SetBinContent(k+1, (fTailX ? fTailX->getVal() : 0.0));
      hFTailY_vsCos->SetBinContent(k+1, (fTailY ? fTailY->getVal() : 0.0));
      if (fTailX) hFTailX_vsCos->SetBinError(k+1, fTailX->getError());
      if (fTailY) hFTailY_vsCos->SetBinError(k+1, fTailY->getError());

      if (savePng) {
        RooPlot* fx = x.frame(Range("fit"));
        RooPlot* fy = y.frame(Range("fit"));

        data2D.plotOn(fx, Name("datax"), CutRange("fit"));
        data2D.plotOn(fy, Name("datay"), CutRange("fit"));

        const double nSS  = N_SS.getVal();
        const double nSB  = N_SB.getVal();
        const double nBS  = N_BS.getVal();
        const double nBB  = N_BB.getVal();
        const double nTot = nSS + nSB + nBS + nBB;

        model.plotOn(fx, ProjWData(data2D), Range("fit"), NormRange("fit"),
                     Normalization(nTot, RooAbsReal::NumEvent),
                     LineColor(colTot), LineWidth(2), Name("totx"));
        model.plotOn(fy, ProjWData(data2D), Range("fit"), NormRange("fit"),
                     Normalization(nTot, RooAbsReal::NumEvent),
                     LineColor(colTot), LineWidth(2), Name("toty"));

        sig1_ptr->plotOn(fx, Range("fit"), NormRange("fit"),
                         Normalization(nSS, RooAbsReal::NumEvent),
                         LineColor(colSS), LineStyle(kSolid), LineWidth(2), Name("ssx"));
        sig1_ptr->plotOn(fx, Range("fit"), NormRange("fit"),
                         Normalization(nSB, RooAbsReal::NumEvent),
                         LineColor(colSB), LineStyle(kDashed), LineWidth(2), Name("sbx"));
        bkg1_ptr->plotOn(fx, Range("fit"), NormRange("fit"),
                         Normalization(nBS, RooAbsReal::NumEvent),
                         LineColor(colBS), LineStyle(kSolid), LineWidth(2), Name("bsx"));
        bkg1_ptr->plotOn(fx, Range("fit"), NormRange("fit"),
                         Normalization(nBB, RooAbsReal::NumEvent),
                         LineColor(colBB), LineStyle(kDashDotted), LineWidth(2), Name("bbx"));

        sig2_ptr->plotOn(fy, Range("fit"), NormRange("fit"),
                         Normalization(nSS, RooAbsReal::NumEvent),
                         LineColor(colSS), LineStyle(kSolid), LineWidth(2), Name("ssy"));
        bkg2_ptr->plotOn(fy, Range("fit"), NormRange("fit"),
                         Normalization(nSB, RooAbsReal::NumEvent),
                         LineColor(colSB), LineStyle(kDashed), LineWidth(2), Name("sby"));
        sig2_ptr->plotOn(fy, Range("fit"), NormRange("fit"),
                         Normalization(nBS, RooAbsReal::NumEvent),
                         LineColor(colBS), LineStyle(kSolid), LineWidth(2), Name("bsy"));
        bkg2_ptr->plotOn(fy, Range("fit"), NormRange("fit"),
                         Normalization(nBB, RooAbsReal::NumEvent),
                         LineColor(colBB), LineStyle(kDashDotted), LineWidth(2), Name("bby"));

        TGraphErrors* grRx = MakeRatioGraph(fx, "datax", "totx");
        TGraphErrors* grRy = MakeRatioGraph(fy, "datay", "toty");

        TCanvas* c = new TCanvas(Form("c%s", tag.Data()), "diagnostic", 1800, 700);
        c->Divide(2,1);

        c->cd(1);
        TPad* p1_top = new TPad(Form("p1_top%s",tag.Data()),"",0,0.30,1,1);
        TPad* p1_bot = new TPad(Form("p1_bot%s",tag.Data()),"",0,0.00,1,0.30);
        p1_top->SetBottomMargin(0.0);
        p1_bot->SetTopMargin(0.0);
        p1_bot->SetBottomMargin(0.25);
        p1_top->SetLogy();
        p1_top->Draw();
        p1_bot->Draw();
        p1_top->cd();
        fx->GetXaxis()->SetTitle("m_{#Lambda1} (GeV/c^{2})");
        fx->Draw();

        {
          TLegend* leg = new TLegend(0.50,0.45,0.88,0.88);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->AddEntry(fx->findObject("totx"), "Total", "l");
          leg->AddEntry(fx->findObject("ssx"),  "SS", "l");
          leg->AddEntry(fx->findObject("sbx"),  "SB", "l");
          leg->AddEntry(fx->findObject("bsx"),  "BS", "l");
          leg->AddEntry(fx->findObject("bbx"),  "BB", "l");
          leg->AddEntry((TObject*)0, Form("EDM=%.3g  status=%d  cov=%d", res->edm(), res->status(), res->covQual()), "");
          leg->AddEntry((TObject*)0, Form("#mu_{x}=%.6f  #sigma^{eff}_{x}=%.6f", mx, sx), "");
          leg->AddEntry((TObject*)0, Form("#mu_{y}=%.6f  #sigma^{eff}_{y}=%.6f", my, sy), "");
          leg->AddEntry((TObject*)0, Form("Y_{SS}^{win}(#pm%.0f#sigma)=%.0f", nSigmaWin, YSS_win), "");
          leg->Draw();
        }
      

        p1_bot->cd();
        p1_bot->SetGridy();
        auto* fr1 = p1_bot->DrawFrame(FIT_MIN, 0.8, FIT_MAX, 1.2);
        fr1->GetYaxis()->SetTitle("Data/Fit");
        fr1->GetXaxis()->SetTitle("m_{#Lambda1} (GeV/c^{2})");
        fr1->GetXaxis()->SetTitleSize(0.1);
        fr1->GetXaxis()->SetLabelSize(0.08);
        fr1->GetYaxis()->SetTitleSize(0.1);
        fr1->GetYaxis()->SetLabelSize(0.08);
       
        if (grRx) grRx->Draw("P SAME");
        {
          TLine* lx1 = new TLine(xlow_bin,0.9,xlow_bin,1.1);
          TLine* lx2 = new TLine(xhi_bin ,0.9,xhi_bin ,1.1);
          lx1->SetLineStyle(kDashed);
          lx2->SetLineStyle(kDashed);
          lx1->Draw("SAME");
          lx2->Draw("SAME");
        }

        c->cd(2);
        TPad* p2_top = new TPad(Form("p2_top%s",tag.Data()),"",0,0.30,1,1);
        TPad* p2_bot = new TPad(Form("p2_bot%s",tag.Data()),"",0,0.00,1,0.30);
        p2_top->SetBottomMargin(0.0);
        p2_bot->SetTopMargin(0.0);
        p2_bot->SetBottomMargin(0.25);
        p2_top->SetLogy();
        p2_top->Draw();
        p2_bot->Draw();

              p2_top->cd();
        fy->GetXaxis()->SetTitle("m_{#Lambda2} (GeV/c^{2})");
       
        fy->Draw();

        {
          TLegend* leg = new TLegend(0.50,0.45,0.88,0.88);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->AddEntry(fy->findObject("toty"), "Total", "l");
          leg->AddEntry(fy->findObject("ssy"),  "SS", "l");
          leg->AddEntry(fy->findObject("sby"),  "SB", "l");
          leg->AddEntry(fy->findObject("bsy"),  "BS", "l");
          leg->AddEntry(fy->findObject("bby"),  "BB", "l");
          leg->AddEntry((TObject*)0, Form("EDM=%.3g  status=%d  cov=%d", res->edm(), res->status(), res->covQual()), "");
          leg->AddEntry((TObject*)0, Form("#mu_{x}=%.6f  #sigma^{eff}_{x}=%.6f", mx, sx), "");
          leg->AddEntry((TObject*)0, Form("#mu_{y}=%.6f  #sigma^{eff}_{y}=%.6f", my, sy), "");
          leg->AddEntry((TObject*)0, Form("Y_{SS}^{win}(#pm%.0f#sigma)=%.0f", nSigmaWin, YSS_win), "");
          leg->Draw();
        }

        p2_bot->cd();
        p2_bot->SetGridy();
        auto* fr2 = p2_bot->DrawFrame(FIT_MIN, 0.8, FIT_MAX, 1.2);
        fr2->GetYaxis()->SetTitle("Data/Fit");
        fr2->GetXaxis()->SetTitle("m_{#Lambda2} (GeV/c^{2})");
        fr2->GetXaxis()->SetTitleSize(0.1);
        fr2->GetYaxis()->SetTitleSize(0.1);
        fr2->GetXaxis()->SetLabelSize(0.08);
        fr2->GetYaxis()->SetLabelSize(0.08);
      
        if (grRy) grRy->Draw("P SAME");
        {
          TLine* ly1 = new TLine(ylow_bin,0.9,ylow_bin,1.1);
          TLine* ly2 = new TLine(yhi_bin ,0.9,yhi_bin ,1.1);
          ly1->SetLineStyle(kDashed);
          ly2->SetLineStyle(kDashed);
          ly1->Draw("SAME");
          ly2->Draw("SAME");
        }

        TString png = Form("%s/diag_Rbin%d_cos%02d.png", outDir.Data(), iR, k+1);
        c->SaveAs(png);
        fout->cd();
        c->Write();

        delete fx;
        delete fy;
        delete grRx;
        delete grRy;
        delete c;
      }

      delete h2;
      axCos->SetRange(0,0);
    }
    if (doME && hsME) {
      for (int k = 0; k < nCosBins; ++k) {

        const double cLo = cosEdges[k];
        const double cHi = cosEdges[k + 1];

        const int bLo = findBinForEdgeLow(axCosm, cLo);
        const int bHi = findBinForEdgeHigh(axCosm, cHi);
        axCosm->SetRange(bLo, bHi);
        axCosm->SetBit(TAxis::kAxisRange);

        TH2D* h2 = dynamic_cast<TH2D*>(hsME->Projection(AX_M2, AX_M1, "E"));
        if (!h2) { axCosm->SetRange(0,0); continue; }
        h2->SetDirectory(nullptr);
        h2->Sumw2();

        const double sumW = h2->Integral();
        if (sumW <= 0) { delete h2; axCosm->SetRange(0,0); continue; }

        const TString tag = Form("_ME_R%d_c%02d", iR, k+1);

        RooRealVar x(Form("x%s", tag.Data()), "m_{#Lambda1}", VAR_MIN, VAR_MAX);
        RooRealVar y(Form("y%s", tag.Data()), "m_{#Lambda2}", VAR_MIN, VAR_MAX);
        x.setRange("fit", FIT_MIN, FIT_MAX);
        y.setRange("fit", FIT_MIN, FIT_MAX);

        RooDataHist data2D(Form("data2D%s", tag.Data()), "data2D", RooArgList(x, y), h2);

        RooRealVar aTail(Form("aTail%s", tag.Data()), "aTail", 3.0, 2.5, 8.0);
        RooRealVar nTail(Form("nTail%s", tag.Data()), "nTail", 10.0, 2.0, 20.0);

        std::vector<std::unique_ptr<RooAbsArg>> ownedX, ownedY;
        RooRealVar *meanX=nullptr, *meanY=nullptr;
        RooAbsReal *sigEffX=nullptr, *sigEffY=nullptr;
        RooRealVar *fTailX=nullptr, *fTailY=nullptr;

        RooAbsPdf* sig1_ptr = buildSig1D(signalOpt, x, Form("X%s",tag.Data()), aTail, nTail, ownedX, meanX, sigEffX, fTailX);
        RooAbsPdf* sig2_ptr = buildSig1D(signalOpt, y, Form("Y%s",tag.Data()), aTail, nTail, ownedY, meanY, sigEffY, fTailY);

        std::vector<std::unique_ptr<RooAbsArg>> keepBx, keepBy;
        auto bkg1_ptr = buildBkg1D(bkgOpt, x, Form("bkg1%s",tag.Data()), keepBx);
        auto bkg2_ptr = buildBkg1D(bkgOpt, y, Form("bkg2%s",tag.Data()), keepBy);

        RooProdPdf pdfSS(Form("pdfSS%s", tag.Data()), "SS", RooArgList(*sig1_ptr, *sig2_ptr));
        RooProdPdf pdfSB(Form("pdfSB%s", tag.Data()), "SB", RooArgList(*sig1_ptr, *bkg2_ptr));
        RooProdPdf pdfBS(Form("pdfBS%s", tag.Data()), "BS", RooArgList(*bkg1_ptr, *sig2_ptr));
        RooProdPdf pdfBB(Form("pdfBB%s", tag.Data()), "BB", RooArgList(*bkg1_ptr, *bkg2_ptr));

        const double Nguess = std::max(1.0, sumW);
        RooRealVar N_SS(Form("N_SS%s", tag.Data()), "N_SS", 0.80*Nguess, 0.0, 20.0*Nguess);
        RooRealVar N_SB(Form("N_SB%s", tag.Data()), "N_SB", 0.07*Nguess, 0.0, 20.0*Nguess);
        RooRealVar N_BS(Form("N_BS%s", tag.Data()), "N_BS", 0.07*Nguess, 0.0, 20.0*Nguess);
        RooRealVar N_BB(Form("N_BB%s", tag.Data()), "N_BB", 0.06*Nguess, 0.0, 20.0*Nguess);

        RooAddPdf model(Form("model%s", tag.Data()), "model",
                        RooArgList(pdfSS,pdfSB,pdfBS,pdfBB),
                        RooArgList(N_SS,N_SB,N_BS,N_BB));

        std::unique_ptr<RooAbsReal> nll(model.createNLL(data2D, Extended(true), Range("fit"), Offset(true)));
        if (!nll) { delete h2; axCosm->SetRange(0,0); continue; }

        RooMinimizer minim(*nll);
        minim.setMinimizerType("Minuit2");
        minim.setStrategy(2);
        minim.setPrintLevel(interactive ? 1 : -1);
        minim.setMaxFunctionCalls(200000);
        minim.setMaxIterations(200000);
        minim.setEps(1e-4);
        minim.setOffsetting(true);
        minim.optimizeConst(2);
        minim.minimize("Minuit2","simplex");
        minim.minimize("Minuit2","migrad");
        minim.hesse();

        std::unique_ptr<RooFitResult> res(minim.save());
        if (!res) { delete h2; axCosm->SetRange(0,0); continue; }

        const double mx_fit = meanX ? meanX->getVal() : 1.115;
        const double my_fit = meanY ? meanY->getVal() : 1.115;
        const double sx_fit = sigEffX ? sigEffX->getVal() : 0.002;
        const double sy_fit = sigEffY ? sigEffY->getVal() : 0.002;

        const double exSig_fit = GetPropErrSafe(sigEffX, res.get());
        const double eySig_fit = GetPropErrSafe(sigEffY, res.get());

        const double xlow_se = hWinXlow_SE ? hWinXlow_SE->GetBinContent(k+1) : 0.0;
        const double xhi_se  = hWinXhi_SE  ? hWinXhi_SE ->GetBinContent(k+1) : 0.0;
        const double ylow_se = hWinYlow_SE ? hWinYlow_SE->GetBinContent(k+1) : 0.0;
        const double yhi_se  = hWinYhi_SE  ? hWinYhi_SE ->GetBinContent(k+1) : 0.0;

        const bool validSEwin =
          (xhi_se > xlow_se) && (yhi_se > ylow_se) &&
          (xlow_se >= VAR_MIN) && (xhi_se <= VAR_MAX) &&
          (ylow_se >= VAR_MIN) && (yhi_se <= VAR_MAX);

        double xlow = 0.0, xhi = 0.0;
        double ylow = 0.0, yhi = 0.0;
        double xlow_bin = 0.0, xhi_bin = 0.0;
        double ylow_bin = 0.0, yhi_bin = 0.0;

        if (useManualWindowReference) {
          BuildSignalWindowAligned(h2,
                                   manualMeanX, manualSigmaX,
                                   manualMeanY, manualSigmaY,
                                   nSigmaWin,
                                   xlow, xhi, ylow, yhi,
                                   xlow_bin, xhi_bin, ylow_bin, yhi_bin);
        } else if (meUseSEWindow && validSEwin) {
          xlow = xlow_se; xhi = xhi_se;
          ylow = ylow_se; yhi = yhi_se;

          auto axHX = h2->GetXaxis();
          auto axHY = h2->GetYaxis();

          const int bx1_sig = std::max(1, std::min(axHX->FindBin(xlow), axHX->GetNbins()));
          const int bx2_sig = std::max(1, std::min(axHX->FindBin(xhi - 1e-12), axHX->GetNbins()));
          const int by1_sig = std::max(1, std::min(axHY->FindBin(ylow), axHY->GetNbins()));
          const int by2_sig = std::max(1, std::min(axHY->FindBin(yhi - 1e-12), axHY->GetNbins()));

          xlow_bin = axHX->GetBinLowEdge(std::min(bx1_sig, bx2_sig));
          xhi_bin  = axHX->GetBinUpEdge (std::max(bx1_sig, bx2_sig));
          ylow_bin = axHY->GetBinLowEdge(std::min(by1_sig, by2_sig));
          yhi_bin  = axHY->GetBinUpEdge (std::max(by1_sig, by2_sig));
        } else {
          BuildSignalWindowAligned(h2,
                                   mx_fit, sx_fit,
                                   my_fit, sy_fit,
                                   nSigmaWin,
                                   xlow, xhi, ylow, yhi,
                                   xlow_bin, xhi_bin, ylow_bin, yhi_bin);
        }

        x.setRange("sigWin", xlow_bin, xhi_bin);
        y.setRange("sigWin", ylow_bin, yhi_bin);

        if (hWinXlow_used_ME) hWinXlow_used_ME->SetBinContent(k+1, xlow_bin);
        if (hWinXhi_used_ME)  hWinXhi_used_ME ->SetBinContent(k+1, xhi_bin);
        if (hWinYlow_used_ME) hWinYlow_used_ME->SetBinContent(k+1, ylow_bin);
        if (hWinYhi_used_ME)  hWinYhi_used_ME ->SetBinContent(k+1, yhi_bin);

        RooArgSet normSet(x,y);

        if (hNSS_ME) {
          hNSS_ME->SetBinContent(k+1, N_SS.getVal());
          hNSS_ME->SetBinError  (k+1, N_SS.getError());
        }

        auto intSS_win = std::unique_ptr<RooAbsReal>(pdfSS.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
        auto intSB_win = std::unique_ptr<RooAbsReal>(pdfSB.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
        auto intBS_win = std::unique_ptr<RooAbsReal>(pdfBS.createIntegral(normSet, NormSet(normSet), Range("sigWin")));
        auto intBB_win = std::unique_ptr<RooAbsReal>(pdfBB.createIntegral(normSet, NormSet(normSet), Range("sigWin")));

        const double YSS_win = N_SS.getVal() * (intSS_win ? intSS_win->getVal() : 0.0);
        double eYSS_win = 0.0;
        if (intSS_win && N_SS.getVal() > 0 && intSS_win->getVal() > 0) {
          eYSS_win = YSS_win * N_SS.getError() / std::max(1e-12, N_SS.getVal());
        }

        if (hYieldSS_Model_ME) {
          hYieldSS_Model_ME->SetBinContent(k+1, YSS_win);
          hYieldSS_Model_ME->SetBinError  (k+1, eYSS_win);
        }

        const double YSB_win = N_SB.getVal() * (intSB_win ? intSB_win->getVal() : 0.0);
        const double YBS_win = N_BS.getVal() * (intBS_win ? intBS_win->getVal() : 0.0);
        const double YBB_win = N_BB.getVal() * (intBB_win ? intBB_win->getVal() : 0.0);

        double eYSB_win = 0.0;
        double eYBS_win = 0.0;
        double eYBB_win = 0.0;

        if (intSB_win && N_SB.getVal() > 0 && intSB_win->getVal() > 0) eYSB_win = YSB_win * N_SB.getError() / std::max(1e-12, N_SB.getVal());
        if (intBS_win && N_BS.getVal() > 0 && intBS_win->getVal() > 0) eYBS_win = YBS_win * N_BS.getError() / std::max(1e-12, N_BS.getVal());
        if (intBB_win && N_BB.getVal() > 0 && intBB_win->getVal() > 0) eYBB_win = YBB_win * N_BB.getError() / std::max(1e-12, N_BB.getVal());

        if (hYieldSBsum_Model_ME) {
          hYieldSBsum_Model_ME->SetBinContent(k+1, YSB_win + YBS_win);
          hYieldSBsum_Model_ME->SetBinError  (k+1, std::sqrt(eYSB_win*eYSB_win + eYBS_win*eYBS_win));
        }
        if (hYieldBB_Model_ME) {
          hYieldBB_Model_ME->SetBinContent(k+1, YBB_win);
          hYieldBB_Model_ME->SetBinError  (k+1, eYBB_win);
        }

        RooFormulaVar fracSSwin(Form("fracSSwin%s",tag.Data()),
                                "(@0*@1)/(@0*@1 + @2*@3 + @4*@5 + @6*@7)",
                                RooArgList(N_SS, *intSS_win, N_SB, *intSB_win, N_BS, *intBS_win, N_BB, *intBB_win));

        const double fss  = fracSSwin.getVal();
        const double efss = fracSSwin.getPropagatedError(*res);

        if (hFracSSwin_ME) {
          hFracSSwin_ME->SetBinContent(k+1, fss);
          hFracSSwin_ME->SetBinError  (k+1, efss);
        }

        double Ndata_win = 0.0, eNdata_win = 0.0;
        IntegralAndErrorH2Window(h2, xlow_bin, xhi_bin, ylow_bin, yhi_bin, Ndata_win, eNdata_win);
        if (hDataWin_vsCos_ME) {
          hDataWin_vsCos_ME->SetBinContent(k+1, Ndata_win);
          hDataWin_vsCos_ME->SetBinError  (k+1, eNdata_win);
        }

        {
          const double xSB1_lo = sbLeftLow;
          const double xSB1_hi = sbLeftHigh;
          const double xSB2_lo = sbRightLow;
          const double xSB2_hi = sbRightHigh;

          const double ySB1_lo = sbLeftLow;
          const double ySB1_hi = sbLeftHigh;
          const double ySB2_lo = sbRightLow;
          const double ySB2_hi = sbRightHigh;

          const bool xSBinside = (xSB1_lo >= VAR_MIN) && (xSB2_hi <= VAR_MAX);
          const bool ySBinside = (ySB1_lo >= VAR_MIN) && (ySB2_hi <= VAR_MAX);

          if (xSBinside && ySBinside) {
            const double sigWidthX = xhi_bin - xlow_bin;
            const double sigWidthY = yhi_bin - ylow_bin;

            const double totSbWidthX = (xSB1_hi - xSB1_lo) + (xSB2_hi - xSB2_lo);
            const double totSbWidthY = (ySB1_hi - ySB1_lo) + (ySB2_hi - ySB2_lo);

            const double rX = sigWidthX / std::max(1e-12, totSbWidthX);
            const double rY = sigWidthY / std::max(1e-12, totSbWidthY);

            ValErr Ame{Ndata_win, eNdata_win};

            ValErr Braw = AddVE(
              IntegralVE(h2, xlow_bin, xhi_bin, ySB1_lo, ySB1_hi),
              IntegralVE(h2, xlow_bin, xhi_bin, ySB2_lo, ySB2_hi)
            );
            ValErr Craw = AddVE(
              IntegralVE(h2, xSB1_lo, xSB1_hi, ylow_bin, yhi_bin),
              IntegralVE(h2, xSB2_lo, xSB2_hi, ylow_bin, yhi_bin)
            );
            ValErr Draw = AddVE(
              AddVE(IntegralVE(h2, xSB1_lo, xSB1_hi, ySB1_lo, ySB1_hi),
                    IntegralVE(h2, xSB1_lo, xSB1_hi, ySB2_lo, ySB2_hi)),
              AddVE(IntegralVE(h2, xSB2_lo, xSB2_hi, ySB1_lo, ySB1_hi),
                    IntegralVE(h2, xSB2_lo, xSB2_hi, ySB2_lo, ySB2_hi))
            );

            ValErr B = ScaleVE(Braw, rY);
            ValErr C = ScaleVE(Craw, rX);
            ValErr D = ScaleVE(Draw, rX*rY);

            ValErr Bkg = SubVE(AddVE(B, C), D);
            ValErr Sig = MakeSigFromAandBkg(Ame, Bkg);
            ValErr FracSig = MakeFracSigFromAandBkg(Ame, Bkg);

            if (hBkgWin_BCD_vsCos_ME) {
              hBkgWin_BCD_vsCos_ME->SetBinContent(k+1, Bkg.v);
              hBkgWin_BCD_vsCos_ME->SetBinError  (k+1, Bkg.e);
                /*
              hBkgWin_B_vsCos_ME->SetBinContent(k+1, B.v);
              hBkgWin_B_vsCos_ME->SetBinError  (k+1, B.e);
              hBkgWin_C_vsCos_ME->SetBinContent(k+1, C.v);
              hBkgWin_C_vsCos_ME->SetBinError  (k+1, C.e);
              hBkgWin_D_vsCos_ME->SetBinContent(k+1, D.v);
              hBkgWin_D_vsCos_ME->SetBinError  (k+1, D.e);
              */
             
              hBkgWin_B_vsCos_ME->SetBinContent(k+1, Braw.v);
              hBkgWin_B_vsCos_ME->SetBinError  (k+1, Braw.e);
              hBkgWin_C_vsCos_ME->SetBinContent(k+1, Craw.v);
              hBkgWin_C_vsCos_ME->SetBinError  (k+1, Craw.e);
              hBkgWin_D_vsCos_ME->SetBinContent(k+1, Draw.v);
              hBkgWin_D_vsCos_ME->SetBinError  (k+1, Draw.e);
              

            }
            if (hSigWin_BCD_vsCos_ME) {
              hSigWin_BCD_vsCos_ME->SetBinContent(k+1, Sig.v);
              hSigWin_BCD_vsCos_ME->SetBinError  (k+1, Sig.e);
            }
            if (hFracSig_BCD_vsCos_ME) {
              hFracSig_BCD_vsCos_ME->SetBinContent(k+1, FracSig.v);
              hFracSig_BCD_vsCos_ME->SetBinError  (k+1, FracSig.e);
            }

            if (savePng) {
              TH2D* h2ABCD = (TH2D*)h2->Clone(Form("h2ABCD_ME_R%d_c%02d", iR, k+1));
              h2ABCD->SetDirectory(nullptr);

              TCanvas* cABCD = MakeABCDCanvas(
                h2ABCD,
                Form("cABCD_ME_R%d_c%02d", iR, k+1),
                Form("ME ABCD regions: Rbin %d, cos bin %d", iR, k+1),
                xlow_bin, xhi_bin, ylow_bin, yhi_bin,
                xSB1_lo, xSB1_hi, xSB2_lo, xSB2_hi,
                ySB1_lo, ySB1_hi, ySB2_lo, ySB2_hi,
                true
              );

              if (cABCD) {
                TString png2D = Form("%s/ABCD_ME_Rbin%d_cos%02d.png", outDir.Data(), iR, k+1);
                cABCD->SaveAs(png2D);
                fout->cd();
                h2ABCD->Write();
                cABCD->Write();
                delete cABCD;
              }
              delete h2ABCD;
            }
          }
        }

        const double Ysig = fss * Ndata_win;
        double eYsig = std::abs(Ndata_win) * efss;
        if (propagateFracErr) {
          eYsig = std::sqrt((Ndata_win*efss)*(Ndata_win*efss) + (fss*eNdata_win)*(fss*eNdata_win));
        }

        if (hYieldSS_ME) {
          hYieldSS_ME->SetBinContent(k+1, Ysig);
          hYieldSS_ME->SetBinError  (k+1, eYsig);
        }

        if (hMeanX_ME) {
          hMeanX_ME->SetBinContent(k+1, mx_fit);
          if (meanX) hMeanX_ME->SetBinError(k+1, meanX->getError());
        }
        if (hMeanY_ME) {
          hMeanY_ME->SetBinContent(k+1, my_fit);
          if (meanY) hMeanY_ME->SetBinError(k+1, meanY->getError());
        }
        if (hSigmaEffX_ME) {
          hSigmaEffX_ME->SetBinContent(k+1, sx_fit);
          hSigmaEffX_ME->SetBinError(k+1, exSig_fit);
        }
        if (hSigmaEffY_ME) {
          hSigmaEffY_ME->SetBinContent(k+1, sy_fit);
          hSigmaEffY_ME->SetBinError(k+1, eySig_fit);
        }

        if (savePng) {
          RooPlot* fx = x.frame(Range("fit"));
          RooPlot* fy = y.frame(Range("fit"));

          data2D.plotOn(fx, Name("datax"), CutRange("fit"));
          data2D.plotOn(fy, Name("datay"), CutRange("fit"));

          const double nSS  = N_SS.getVal();
          const double nSB  = N_SB.getVal();
          const double nBS  = N_BS.getVal();
          const double nBB  = N_BB.getVal();
          const double nTot = nSS+nSB+nBS+nBB;

          model.plotOn(fx, ProjWData(data2D), Range("fit"), NormRange("fit"),
                       Normalization(nTot, RooAbsReal::NumEvent),
                       LineColor(colTot), LineWidth(2), Name("totx"));
          model.plotOn(fy, ProjWData(data2D), Range("fit"), NormRange("fit"),
                       Normalization(nTot, RooAbsReal::NumEvent),
                       LineColor(colTot), LineWidth(2), Name("toty"));

          sig1_ptr->plotOn(fx, Range("fit"), NormRange("fit"), Normalization(nSS, RooAbsReal::NumEvent),
                           LineColor(colSS), LineStyle(kSolid), LineWidth(2), Name("ssx"));
          sig1_ptr->plotOn(fx, Range("fit"), NormRange("fit"), Normalization(nSB, RooAbsReal::NumEvent),
                           LineColor(colSB), LineStyle(kDashed), LineWidth(2), Name("sbx"));
          bkg1_ptr->plotOn(fx, Range("fit"), NormRange("fit"), Normalization(nBS, RooAbsReal::NumEvent),
                           LineColor(colBS), LineStyle(kSolid), LineWidth(2), Name("bsx"));
          bkg1_ptr->plotOn(fx, Range("fit"), NormRange("fit"), Normalization(nBB, RooAbsReal::NumEvent),
                           LineColor(colBB), LineStyle(kDashDotted), LineWidth(2), Name("bbx"));

          sig2_ptr->plotOn(fy, Range("fit"), NormRange("fit"), Normalization(nSS, RooAbsReal::NumEvent),
                           LineColor(colSS), LineStyle(kSolid), LineWidth(2), Name("ssy"));
          bkg2_ptr->plotOn(fy, Range("fit"), NormRange("fit"), Normalization(nSB, RooAbsReal::NumEvent),
                           LineColor(colSB), LineStyle(kDashed), LineWidth(2), Name("sby"));
          sig2_ptr->plotOn(fy, Range("fit"), NormRange("fit"), Normalization(nBS, RooAbsReal::NumEvent),
                           LineColor(colBS), LineStyle(kSolid), LineWidth(2), Name("bsy"));
          bkg2_ptr->plotOn(fy, Range("fit"), NormRange("fit"), Normalization(nBB, RooAbsReal::NumEvent),
                           LineColor(colBB), LineStyle(kDashDotted), LineWidth(2), Name("bby"));

          TGraphErrors* grRx = MakeRatioGraph(fx, "datax", "totx");
          TGraphErrors* grRy = MakeRatioGraph(fy, "datay", "toty");

          TCanvas* c = new TCanvas(Form("c%s", tag.Data()), "diagnostic", 1800, 700);
          c->Divide(2,1);

          c->cd(1);
          TPad* p1_top = new TPad(Form("p1_top%s",tag.Data()),"",0,0.30,1,1);
          TPad* p1_bot = new TPad(Form("p1_bot%s",tag.Data()),"",0,0.00,1,0.30);
          p1_top->SetBottomMargin(0.0);
          p1_bot->SetTopMargin(0.0);
          p1_bot->SetBottomMargin(0.25);
          p1_top->SetLogy();
          p1_top->Draw();
          p1_bot->Draw();

                 p1_top->cd();
          fx->GetXaxis()->SetTitle("m_{#Lambda1} (GeV/c^{2})");
          fx->Draw();

          {
            TLegend* leg = new TLegend(0.50,0.45,0.88,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(fx->findObject("totx"), "Total", "l");
            leg->AddEntry(fx->findObject("ssx"),  "SS", "l");
            leg->AddEntry(fx->findObject("sbx"),  "SB", "l");
            leg->AddEntry(fx->findObject("bsx"),  "BS", "l");
            leg->AddEntry(fx->findObject("bbx"),  "BB", "l");
            leg->AddEntry((TObject*)0, Form("EDM=%.3g  status=%d  cov=%d", res->edm(), res->status(), res->covQual()), "");
            leg->AddEntry((TObject*)0, Form("#mu_{x}=%.6f  #sigma^{eff}_{x}=%.6f", mx_fit, sx_fit), "");
            leg->AddEntry((TObject*)0, Form("#mu_{y}=%.6f  #sigma^{eff}_{y}=%.6f", my_fit, sy_fit), "");
            leg->AddEntry((TObject*)0, Form("Y_{SS}^{win}(#pm%.0f#sigma)=%.0f", nSigmaWin, YSS_win), "");
            leg->Draw();
          }

          p1_bot->cd();
          p1_bot->SetGridy();
          auto* fr1 = p1_bot->DrawFrame(FIT_MIN, 0.8, FIT_MAX, 1.2);
          fr1->GetYaxis()->SetTitle("Data/Fit");
          fr1->GetXaxis()->SetTitle("m_{#Lambda1} (GeV/c^{2})");
          fr1->GetXaxis()->SetTitleSize(0.1);
          fr1->GetYaxis()->SetTitleSize(0.1);
          fr1->GetXaxis()->SetLabelSize(0.08);
          fr1->GetYaxis()->SetLabelSize(0.08);
          if (grRx) grRx->Draw("P SAME");
          {
            TLine* lx1 = new TLine(xlow_bin,0.9,xlow_bin,1.1);
            TLine* lx2 = new TLine(xhi_bin ,0.9,xhi_bin ,1.1);
            lx1->SetLineStyle(kDashed);
            lx2->SetLineStyle(kDashed);
            lx1->Draw("SAME");
            lx2->Draw("SAME");
          }

          c->cd(2);
          TPad* p2_top = new TPad(Form("p2_top%s",tag.Data()),"",0,0.30,1,1);
          TPad* p2_bot = new TPad(Form("p2_bot%s",tag.Data()),"",0,0.00,1,0.30);
          p2_top->SetBottomMargin(0.0);
          p2_bot->SetTopMargin(0.0);
          p2_bot->SetBottomMargin(0.25);
          p2_top->SetLogy();
          p2_top->Draw();
          p2_bot->Draw();
          p2_top->cd();
          fy->GetXaxis()->SetTitle("m_{#Lambda2} (GeV/c^{2})");
          fy->Draw();

          {
            TLegend* leg = new TLegend(0.50,0.45,0.88,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(fy->findObject("toty"), "Total", "l");
            leg->AddEntry(fy->findObject("ssy"),  "SS", "l");
            leg->AddEntry(fy->findObject("sby"),  "SB", "l");
            leg->AddEntry(fy->findObject("bsy"),  "BS", "l");
            leg->AddEntry(fy->findObject("bby"),  "BB", "l");
            leg->AddEntry((TObject*)0, Form("EDM=%.3g  status=%d  cov=%d", res->edm(), res->status(), res->covQual()), "");
            leg->AddEntry((TObject*)0, Form("#mu_{x}=%.6f  #sigma^{eff}_{x}=%.6f", mx_fit, sx_fit), "");
            leg->AddEntry((TObject*)0, Form("#mu_{y}=%.6f  #sigma^{eff}_{y}=%.6f", my_fit, sy_fit), "");
            leg->AddEntry((TObject*)0, Form("Y_{SS}^{win}(#pm%.0f#sigma)=%.0f", nSigmaWin, YSS_win), "");
            leg->Draw();
          }

          p2_bot->cd();
          p2_bot->SetGridy();
          auto* fr2 = p2_bot->DrawFrame(FIT_MIN, 0.8, FIT_MAX, 1.2);
          fr2->GetYaxis()->SetTitle("Data/Fit");
          fr2->GetXaxis()->SetTitle("m_{#Lambda2} (GeV/c^{2})");
          fr2->GetXaxis()->SetTitleSize(0.1);
          fr2->GetYaxis()->SetTitleSize(0.1);
          fr2->GetXaxis()->SetLabelSize(0.08);
          fr2->GetYaxis()->SetLabelSize(0.08);

          if (grRy) grRy->Draw("P SAME");
          {
            TLine* ly1 = new TLine(ylow_bin,0.9,ylow_bin,1.1);
            TLine* ly2 = new TLine(yhi_bin ,0.9,yhi_bin ,1.1);
            ly1->SetLineStyle(kDashed);
            ly2->SetLineStyle(kDashed);
            ly1->Draw("SAME");
            ly2->Draw("SAME");
          }

          TString png = Form("%s/diagME_Rbin%d_cos%02d.png", outDir.Data(), iR, k+1);
          c->SaveAs(png);
          fout->cd();
          c->Write();

          delete fx;
          delete fy;
          delete grRx;
          delete grRy;
          delete c;
        }

        delete h2;
        axCosm->SetRange(0,0);
      }
    }


    // calculate the average fs.    
    //============================================================================================================================
     
  
    // 计算加权平均 fs
    double sum_total = 0.0;
    double sum_signal = 0.0;
    double sum_total_err2 = 0.0;
    double sum_signal_err2 = 0.0;
    
    for (int k = 0; k < nCosBins; ++k) {
        double N_total = hDataWin_vsCos->GetBinContent(k+1);
        double err_total = hDataWin_vsCos->GetBinError(k+1);
        double N_signal = hSigWin_BCD_vsCos->GetBinContent(k+1);
        double err_signal = hSigWin_BCD_vsCos->GetBinError(k+1);
        
        // 只累加正数
        if (N_total > 0) {
            sum_total += N_total;
            sum_total_err2 += err_total * err_total;
        }
        if (N_signal > 0) {
            sum_signal += N_signal;
            sum_signal_err2 += err_signal * err_signal;
        }
    }
    
    double fs_avg = (sum_total > 0) ? sum_signal / sum_total : 0.0;
    double err_fs_avg = 0.0;
    if (sum_total > 0 && sum_signal > 0) {
        double rel_err_total = sqrt(sum_total_err2) / sum_total;
        double rel_err_signal = sqrt(sum_signal_err2) / sum_signal;
        err_fs_avg = fs_avg * sqrt(rel_err_signal*rel_err_signal + rel_err_total*rel_err_total);
    }
    
    // 存储到直方图（在循环外创建）
   
    hFsAvg_vs_R->SetBinContent(iR+1, fs_avg);
    hFsAvg_vs_R->SetBinError(iR+1, err_fs_avg);
    

    //============================================================================================================================




    TH1D* hRatioYield_to_NSSwin = MakeRatioHist(hYieldSS, hNSSwin_vsCos, "hRatioYield_to_NSSwin", "Yield / N_{SS}^{win};cos#theta*;Ratio");
    TH1D* hRatioABCD_to_NSSwin  = MakeRatioHist(hSigWin_BCD_vsCos, hNSSwin_vsCos, "hRatioABCD_to_NSSwin", "Yield_{ABCD} / N_{SS}^{win};cos#theta*;Ratio");
    {
      TCanvas* cYieldCmp = new TCanvas(Form("cYieldCompare_R%d", iR), "Yield comparison", 2000, 800);
      cYieldCmp->Divide(2,1);
      cYieldCmp->cd(1); gPad->SetGridy();
      hNSSwin_vsCos->SetMarkerStyle(20); hNSSwin_vsCos->SetMarkerColor(kBlack); hNSSwin_vsCos->SetLineColor(kBlack); hNSSwin_vsCos->SetMarkerSize(1.1);
      hYieldSS->SetMarkerStyle(24); hYieldSS->SetMarkerColor(kRed+1); hYieldSS->SetLineColor(kRed+1); hYieldSS->SetMarkerSize(1.1);
      hYieldSS_Model->SetMarkerStyle(21); hYieldSS_Model->SetMarkerColor(kBlue+1); hYieldSS_Model->SetLineColor(kBlue+1); hYieldSS_Model->SetMarkerSize(1.1);
      hSigWin_BCD_vsCos->SetMarkerStyle(25); hSigWin_BCD_vsCos->SetMarkerColor(kGreen+2); hSigWin_BCD_vsCos->SetLineColor(kGreen+2); hSigWin_BCD_vsCos->SetMarkerSize(1.1);
      hNSSwin_vsCos->SetTitle("Yield comparison in signal window;cos#theta*;Yield");
      hNSSwin_vsCos->Draw("E1"); hYieldSS->Draw("E1 SAME"); hYieldSS_Model->Draw("E1 SAME"); hSigWin_BCD_vsCos->Draw("E1 SAME");
      {
        TLegend* leg = new TLegend(0.45,0.62,0.88,0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
        leg->AddEntry(hNSSwin_vsCos, "N_{SS}^{win}", "lep");
        leg->AddEntry(hYieldSS, "Yield = f_{SS}^{win} #times N_{data}^{win}", "lep");
        leg->AddEntry(hYieldSS_Model, "Yield_{model}", "lep");
        leg->AddEntry(hSigWin_BCD_vsCos, "Yield_{ABCD}", "lep");
        leg->Draw();
      }
      cYieldCmp->cd(2); gPad->SetGridy();
      hRatioYield_to_NSSwin->SetMarkerStyle(24); hRatioYield_to_NSSwin->SetMarkerColor(kRed+1); hRatioYield_to_NSSwin->SetLineColor(kRed+1); hRatioYield_to_NSSwin->SetMarkerSize(1.1);
      hRatioABCD_to_NSSwin->SetMarkerStyle(25); hRatioABCD_to_NSSwin->SetMarkerColor(kGreen+2); hRatioABCD_to_NSSwin->SetLineColor(kGreen+2); hRatioABCD_to_NSSwin->SetMarkerSize(1.1);
      hRatioYield_to_NSSwin->SetTitle("Ratios to N_{SS}^{win};cos#theta*;Ratio");
      hRatioYield_to_NSSwin->GetYaxis()->SetRangeUser(0.95, 1.12);
      hRatioYield_to_NSSwin->Draw("E1"); hRatioABCD_to_NSSwin->Draw("E1 SAME");
      { TLine* l = new TLine(cosEdges[0], 1.0, cosEdges[nCosBins], 1.0); l->SetLineStyle(kDashed); l->SetLineWidth(2); l->SetLineColor(kBlue+1); l->Draw("SAME"); }
      {
        TLegend* leg = new TLegend(0.25,0.68,0.88,0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
        leg->AddEntry(hRatioYield_to_NSSwin, "Yield / N_{SS}^{win}", "lep");
        leg->AddEntry(hRatioABCD_to_NSSwin, "Yield_{ABCD} / N_{SS}^{win}", "lep");
        leg->AddEntry((TObject*)0, "Yield_{model} / N_{SS}^{win} = 1 by construction", "");
        leg->Draw();
      }
      TLatex* tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.032);
      tex->DrawLatex(0.15, 0.15, Form("Rbin:[%0.2f, %0.2f]", REdges[iR], REdges[iR+1]));
      fout->cd(); cYieldCmp->Write();
      if (savePng) cYieldCmp->SaveAs(Form("%s/yieldCompare_Rbin%d.png", outDir.Data(), iR));
      if (savePng) cYieldCmp->SaveAs(Form("%s/yieldCompare_Rbin%d.pdf", outDir.Data(), iR));
      delete cYieldCmp;
    }

    TH1D* hFracRatio_SEoverME = nullptr;
    TF1*  fFracRatio_pol1 = nullptr;
    if (doME && hsME && hFracSSwin_ME) {
      hFracRatio_SEoverME = (TH1D*)hFracSSwin_vsCos->Clone("hFracSSwin_ratio_SEoverME");
      hFracRatio_SEoverME->SetTitle("Ratio f_{SS}^{win}: SE/ME;cos#theta*;ratio");
      hFracRatio_SEoverME->Reset("ICES");
      hFracRatio_SEoverME->Sumw2();
      for (int b=1;b<=nCosBins;++b){
        const double a  = hFracSSwin_vsCos->GetBinContent(b);
        const double ea = hFracSSwin_vsCos->GetBinError(b);
        const double c  = hFracSSwin_ME->GetBinContent(b);
        const double ec = hFracSSwin_ME->GetBinError(b);
        if (c<=0.0) { hFracRatio_SEoverME->SetBinContent(b, 0.0); hFracRatio_SEoverME->SetBinError(b, 0.0); continue; }
        const double r  = a/c;
        const double er = std::sqrt( (ea*ea)/(c*c) + (a*a*ec*ec)/(c*c*c*c) );
        hFracRatio_SEoverME->SetBinContent(b, r);
        hFracRatio_SEoverME->SetBinError(b, er);
      }
      TCanvas* cFracCmp = new TCanvas(Form("cFracCompare_R%d", iR), "SE vs ME fraction", 1200, 500);
      cFracCmp->Divide(2,1);
      cFracCmp->cd(1); gPad->SetGridy(); hFracSSwin_vsCos->SetMarkerStyle(20); hFracSSwin_ME->SetMarkerStyle(24); hFracSSwin_vsCos->Draw("E1"); hFracSSwin_ME->Draw("E1 SAME");
      cFracCmp->cd(2); gPad->SetGridy(); hFracRatio_SEoverME->SetMarkerStyle(21); hFracRatio_SEoverME->Draw("E1");
      fFracRatio_pol1 = new TF1(Form("fFracRatio_pol1_R%d", iR),"pol1", cosEdges[0], cosEdges[nCosBins]);
      int nNonZero = 0; for (int b=1; b<=hFracRatio_SEoverME->GetNbinsX(); ++b) if (hFracRatio_SEoverME->GetBinError(b) > 0.0) nNonZero++;
      if (nNonZero >= 2) { hFracRatio_SEoverME->Fit(fFracRatio_pol1, "Q0"); fFracRatio_pol1->Draw("SAME"); }
      fout->cd(); hFracRatio_SEoverME->Write(); if (nNonZero >= 2) fFracRatio_pol1->Write(); cFracCmp->Write();
      if (savePng) cFracCmp->SaveAs(Form("%s/fracCompare_Rbin%d.png", outDir.Data(), iR));
      delete cFracCmp;

      TCanvas* cFrac_SE_ME = new TCanvas(Form("cFracCompare_SE_ME_R%d", iR), "SE vs ME fraction", 700, 500);
      gPad->SetGridy(); hFracSSwin_vsCos->SetMarkerStyle(20); hFracSSwin_ME->SetMarkerStyle(24); 
      hFracSSwin_vsCos->GetYaxis()->SetRangeUser(0.8, 1.0);
      hFracSSwin_vsCos->Draw("E1"); hFracSSwin_ME->Draw("E1 SAME");

      TLegend* leg = new TLegend(0.45,0.62,0.88,0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
      leg->AddEntry(hFracSSwin_vsCos, "f_{SS}^{win} SE", "lep");
      leg->AddEntry(hFracSSwin_ME, "f_{SS}^{win} ME", "lep");
      leg->Draw("same");
      TLatex* tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.032);
      tex->DrawLatex(0.15, 0.85, Form("Rbin:[%0.2f, %0.2f]", REdges[iR], REdges[iR+1]));
    
      if (savePng) cFrac_SE_ME->SaveAs(Form("%s/fracCompare_SE_ME_Rbin%d.png", outDir.Data(), iR));
      if (savePng) cFrac_SE_ME->SaveAs(Form("%s/fracCompare_SE_ME_Rbin%d.pdf", outDir.Data(), iR));





    }

    fout->cd();
    hYieldSS->Write(); hYieldSS_Model->Write(); hYieldSBsum_Model->Write(); hYieldBB_Model->Write();
    hEDM->Write(); hStatus->Write(); hCovQual->Write(); hMeanX->Write(); hMeanY->Write(); hSigmaEffX->Write(); hSigmaEffY->Write();
    hNSS_vsCos->Write(); hNSSwin_vsCos->Write(); hFracSSwin_vsCos->Write(); hDataWin_vsCos->Write(); hSigWin_BCD_vsCos->Write(); hBkgWin_BCD_vsCos->Write(); hFracSig_BCD_vsCos->Write(); hX_SE_vsCos->Write(); hFTailX_vsCos->Write(); hFTailY_vsCos->Write(); hWinXlow_SE->Write(); hWinXhi_SE->Write(); hWinYlow_SE->Write(); hWinYhi_SE->Write();
    if (hRatioYield_to_NSSwin) hRatioYield_to_NSSwin->Write();
    if (hRatioABCD_to_NSSwin)  hRatioABCD_to_NSSwin->Write();
    if (doME && hsME) {
      hYieldSS_ME->Write(); hYieldSS_Model_ME->Write(); hYieldSBsum_Model_ME->Write(); hYieldBB_Model_ME->Write();
      hNSS_ME->Write(); hFracSSwin_ME->Write(); hDataWin_vsCos_ME->Write();
      if (hSigWin_BCD_vsCos_ME) hSigWin_BCD_vsCos_ME->Write();
      if (hBkgWin_BCD_vsCos_ME) hBkgWin_BCD_vsCos_ME->Write();

      if (hBkgWin_B_vsCos_ME) hBkgWin_B_vsCos_ME->Write();
      if (hBkgWin_C_vsCos_ME)  hBkgWin_C_vsCos_ME->Write();
      if (hBkgWin_D_vsCos_ME)  hBkgWin_D_vsCos_ME->Write();
     
      if (hFracSig_BCD_vsCos_ME) hFracSig_BCD_vsCos_ME->Write();
      hMeanX_ME->Write(); hMeanY_ME->Write(); hSigmaEffX_ME->Write(); hSigmaEffY_ME->Write();
      hWinXlow_used_ME->Write(); hWinXhi_used_ME->Write(); hWinYlow_used_ME->Write(); hWinYhi_used_ME->Write();
    }

    
    fout->Write();

    delete hYieldSS; delete hYieldSS_Model; delete hYieldSBsum_Model; delete hYieldBB_Model;
    delete hEDM; delete hStatus; delete hCovQual; delete hDataWin_vsCos; delete hMeanX; delete hMeanY; delete hSigmaEffX; delete hSigmaEffY;
    delete hNSS_vsCos; delete hNSSwin_vsCos; delete hFracSSwin_vsCos; delete hSigWin_BCD_vsCos; delete hBkgWin_BCD_vsCos; delete hFracSig_BCD_vsCos; delete hX_SE_vsCos; delete hFTailX_vsCos; delete hFTailY_vsCos; delete hWinXlow_SE; delete hWinXhi_SE; delete hWinYlow_SE; delete hWinYhi_SE;
    if (hRatioYield_to_NSSwin) delete hRatioYield_to_NSSwin;
    if (hRatioABCD_to_NSSwin)  delete hRatioABCD_to_NSSwin;
    if (hFracRatio_SEoverME) delete hFracRatio_SEoverME;
    if (fFracRatio_pol1) delete fFracRatio_pol1;
    if (doME && hsME) {
      delete hYieldSS_ME; delete hYieldSS_Model_ME; delete hYieldSBsum_Model_ME; delete hYieldBB_Model_ME;
      delete hNSS_ME; delete hFracSSwin_ME; delete hDataWin_vsCos_ME;
      if (hSigWin_BCD_vsCos_ME) delete hSigWin_BCD_vsCos_ME;
      if (hBkgWin_BCD_vsCos_ME) delete hBkgWin_BCD_vsCos_ME;

      if (hBkgWin_B_vsCos_ME)  delete hBkgWin_B_vsCos_ME;
      if (hBkgWin_C_vsCos_ME)  delete hBkgWin_C_vsCos_ME;
      if (hBkgWin_D_vsCos_ME)  delete hBkgWin_D_vsCos_ME;



      if (hFracSig_BCD_vsCos_ME) delete hFracSig_BCD_vsCos_ME;
      delete hMeanX_ME; delete hMeanY_ME; delete hSigmaEffX_ME; delete hSigmaEffY_ME;
      delete hWinXlow_used_ME; delete hWinXhi_used_ME; delete hWinYlow_used_ME; delete hWinYhi_used_ME;
    }
    fout->Close(); delete fout;
  }

  if (hFsAvg_vs_R) {
      TString summaryRoot = Form("%s/fs_avg.root", outDir.Data());
      TFile* fsum = TFile::Open(summaryRoot, "RECREATE");
      if (fsum && !fsum->IsZombie()) {
          hFsAvg_vs_R->Write();
          fsum->Close();
          delete fsum;
          std::cout << "Saved average f_s to " << summaryRoot << std::endl;
      }
  }

  fin->Close(); delete fin;
  std::cout << "Done. Outputs in: " << outDir << "\n"
            << "ROOT:  " << outDir << "/fit_Rbin*.root\n"
            << "PNGs:  " << outDir << "/diag_Rbin*_cos*.png\n"
            << "Yield compare PNGs: " << outDir << "/yieldCompare_Rbin*.png\n";
  if (doME) std::cout << "ME PNGs: " << outDir << "/diagME_Rbin*_cos*.png\n";
}

