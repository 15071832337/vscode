#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TPad.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <vector>

// ============================================================
// Fit model
//   baseline * (1 + CP * cos(theta*))
// where CP = alpha1 * alpha2 * P
// ============================================================
double Signal1d(double* x, double* par)
{
  return par[0] * (1.0 + par[1] * x[0]); // par[0]=baseline, par[1]=C*P
  //return  (par[0]+ par[1] * x[0]); 
}

struct FitOut {
  std::unique_ptr<TF1> f;
  double baseline = 0, baselineErr = 0;
  double CP = 0, CPerr = 0;
  double chi2ndf = 0;
  double P = 0, Perr = 0;
  int fitStatus = -999;
};

static FitOut FitAndReturn(TH1D* h, double Cscale, const char* fname)
{
  FitOut out;
  if (!h) return out;

  const int nb = h->GetNbinsX();
  int nuse = 0;
  double sum = 0.0;

  for (int b = 1; b <= nb; ++b) {
    const double y = h->GetBinContent(b);
    const double e = h->GetBinError(b);
    if (std::isfinite(y) && std::isfinite(e) && e > 0.0) {
      sum += y;
      ++nuse;
    }
  }
  if (nuse < 2) return out;

  const double mean = (nuse > 0) ? sum / nuse : 1.0;

  out.f = std::make_unique<TF1>(fname, Signal1d, -1.0, 1.0, 2);
  out.f->SetParameter(0, mean);
  out.f->SetParameter(1, 0.0);

  out.f->SetParLimits(0, 0.0, 10.0 * std::max(1.0, mean));
  out.f->SetParLimits(1, -1.0, 1.0);

  out.fitStatus = h->Fit(out.f.get(), "REMS");

  out.baseline    = out.f->GetParameter(0);
  out.baselineErr = out.f->GetParError(0);
  out.CP          = out.f->GetParameter(1);
  out.CPerr       = out.f->GetParError(1);
  out.chi2ndf     = (out.f->GetNDF() > 0) ? out.f->GetChisquare() / out.f->GetNDF() : 0.0;

  if (Cscale != 0.0) {
    
    out.P    = out.CP / (Cscale);
    out.Perr = out.CPerr / std::abs(Cscale);
    
    /*
    out.P    = out.CP / (Cscale*out.baseline);
    double rel_err_p1 = (out.CP != 0.0) ? out.CPerr / std::abs(out.CP) : 0.0;
    double rel_err_p0 = out.baselineErr / out.baseline;
    double rel_err_Cscale = 0.0;  // 如果 Cscale 有误差可以加进来
    
    double rel_err_P = sqrt(rel_err_p1*rel_err_p1 + rel_err_p0*rel_err_p0 + rel_err_Cscale*rel_err_Cscale);
    out.Perr = std::abs(out.P) * rel_err_P;
    */



  }

  return out;
}

static void DrawPanel(TH1D* h, const FitOut& fo, const char* titleLine, double yMin, double yMax)
{
  if (!h) {
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.06);
    lat.DrawLatex(0.20, 0.55, "Missing histogram");
    return;
  }

  h->SetTitle(Form("%s;cos#theta*;SE/ME(norm)", titleLine));
  h->GetYaxis()->SetRangeUser(yMin, yMax);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetXaxis()->SetTitleOffset(1.05);

  h->Draw("E1");
  if (fo.f) fo.f->Draw("SAME");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.050);
  lat.DrawLatex(0.16, 0.86, Form("P = %.4f #pm %.4f", fo.P, fo.Perr));
  lat.DrawLatex(0.16, 0.79, Form("CP = %.4f #pm %.4f", fo.CP, fo.CPerr));
  lat.DrawLatex(0.16, 0.72, Form("Base = %.4f #pm %.4f", fo.baseline, fo.baselineErr));
  lat.DrawLatex(0.16, 0.65, Form("#chi^{2}/NDF = %.2f", fo.chi2ndf));
}

static void UpdateRangeFrom(const TH1D* h, double& ymin, double& ymax)
{
  if (!h) return;
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double y = h->GetBinContent(b);
    const double e = h->GetBinError(b);
    if (!std::isfinite(y) || !std::isfinite(e)) continue;
    ymin = std::min(ymin, y - e);
    ymax = std::max(ymax, y + e);
  }
}

// ============================================================
// Acceptance correction
//   1) scale ME to same integral as SE
//   2) build SE/ME bin-by-bin
// ============================================================
static TH1D* MakeAcceptanceCorrected(const TH1D* hSE_in, const TH1D* hME_in, const char* newname)
{
  if (!hSE_in || !hME_in) return nullptr;

  auto hSE = std::unique_ptr<TH1D>((TH1D*)hSE_in->Clone(Form("%s_tmpSE", newname)));
  auto hME = std::unique_ptr<TH1D>((TH1D*)hME_in->Clone(Form("%s_tmpME", newname)));
  hSE->SetDirectory(nullptr);
  hME->SetDirectory(nullptr);

  const int nb = hSE->GetNbinsX();
  if (hME->GetNbinsX() != nb) {
    std::cerr << "ERROR: SE/ME bin mismatch for " << newname << "\n";
    return nullptr;
  }

  const double intSE = hSE->Integral(1, nb);
  const double intME = hME->Integral(1, nb);
  if (intSE <= 0.0 || intME <= 0.0) {
    std::cerr << "ERROR: non-positive SE/ME integral for " << newname << "\n";
    return nullptr;
  }

  //const double norm =  1.0/intME;
  //const double norm =  intSE/intME;
  
  //hSE->Scale(1.0/intSE);
  hME->Scale(1.0/intME);
  //hME->Scale(norm);

  TH1D* hCorr = (TH1D*)hSE->Clone(newname);
  hCorr->SetDirectory(nullptr);
  hCorr->Reset("ICES");

  for (int b = 1; b <= nb; ++b) {
    const double se  = hSE->GetBinContent(b);
    const double seE = hSE->GetBinError(b);
    const double me  = hME->GetBinContent(b);
    const double meE = hME->GetBinError(b);

    if (!std::isfinite(se) || !std::isfinite(seE) ||
        !std::isfinite(me) || !std::isfinite(meE) || me <= 0.0) {
      hCorr->SetBinContent(b, 0.0);
      hCorr->SetBinError(b, 0.0);
      continue;
    }

    const double val = se / me;

    double err = 0.0;
    if (se > 0.0 && me > 0.0) {
      const double rel2 =
        (seE > 0.0 ? (seE / se) * (seE / se) : 0.0) +
        (meE > 0.0 ? (meE / me) * (meE / me) : 0.0);
      err = std::abs(val) * std::sqrt(rel2);
    } else if (me > 0.0) {
      err = (seE > 0.0) ? std::abs(seE / me) : 0.0;
    }

    hCorr->SetBinContent(b, val);
    hCorr->SetBinError(b, err);
  }
  


  return hCorr;
}

// ============================================================
// Sum helper
// ============================================================
static TH1D* MakeSummedHist(const TH1D* h1, const TH1D* h2, const char* newname)
{
  if (!h1 || !h2) return nullptr;
  if (h1->GetNbinsX() != h2->GetNbinsX()) return nullptr;

  TH1D* out = (TH1D*)h1->Clone(newname);
  out->SetDirectory(nullptr);
  out->Add(h2);
  return out;
}
// ============================================================
// B+C-D
// ============================================================
static TH1D* MakeBPlusCMinusDHist(const TH1D* hB, const TH1D* hC, const TH1D* hD, const char* newname)
{
  if (!hB || !hC || !hD) return nullptr;
  if (hB->GetNbinsX() != hC->GetNbinsX() || hB->GetNbinsX() != hD->GetNbinsX()) return nullptr;
  TH1D* out = (TH1D*)hB->Clone(newname);
  out->SetDirectory(nullptr);
  int nb = hB->GetNbinsX();
  for (int b = 1; b <= nb; ++b) {

    double B = hB->GetBinContent(b);
    double eB = hB->GetBinError(b);
    double C = hC->GetBinContent(b);
    double eC = hC->GetBinError(b);
    double D = hD->GetBinContent(b);
    double eD = hD->GetBinError(b);
    
    // B + C - D
    double val = B + C - D;
    double err = sqrt(eB*eB + eC*eC + eD*eD);

    out->SetBinContent(b, val);
    out->SetBinError(b, err);

  }
  return out;
}
//The weighted method for B+C-D
static TH1D* MakeWeightedAverageHist(const TH1D* h1, const TH1D* h2, const TH1D* h3, const char* newname)
{
  if (!h1 || !h2 || !h3) return nullptr;
  if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsX() != h3->GetNbinsX()) return nullptr;
  
  TH1D* out = (TH1D*)h1->Clone(newname);
  out->SetDirectory(nullptr);
  out->Reset("ICES");
  
  int nb = h1->GetNbinsX();
  for (int b = 1; b <= nb; ++b) {
    double y1 = h1->GetBinContent(b);
    double e1 = h1->GetBinError(b);
    double y2 = h2->GetBinContent(b);
    double e2 = h2->GetBinError(b);
    double y3 = h3->GetBinContent(b);
    double e3 = h3->GetBinError(b);
    
    // 权重 = 1/误差²
    double w1 = (e1 > 0 && std::isfinite(y1)) ? 1.0/(e1*e1) : 0.0;
    double w2 = (e2 > 0 && std::isfinite(y2)) ? 1.0/(e2*e2) : 0.0;
    double w3 = (e3 > 0 && std::isfinite(y3)) ? 1.0/(e3*e3) : 0.0;
    double wsum = w1 + w2 + w3;
    
    if (wsum > 0) {
      double val = (w1*y1 + w2*y2 + w3*y3) / wsum;
      double err = sqrt(1.0 / wsum);
      out->SetBinContent(b, val);
      out->SetBinError(b, err);
    } else {
      out->SetBinContent(b, 0.0);
      out->SetBinError(b, 0.0);
    }
  }
  return out;
}
// only consider B+C, without D
static TH1D* MakeBPlusCHist(const TH1D* hB, const TH1D* hC, const char* newname)
{
  if (!hB || !hC) return nullptr;
  if (hB->GetNbinsX() != hC->GetNbinsX()) return nullptr;
  
  TH1D* out = (TH1D*)hB->Clone(newname);
  out->SetDirectory(nullptr);
  out->Reset("ICES");
  
  int nb = hB->GetNbinsX();
  for (int b = 1; b <= nb; ++b) {
    double B = hB->GetBinContent(b);
    double eB = hB->GetBinError(b);
    double C = hC->GetBinContent(b);
    double eC = hC->GetBinError(b);
    
    double val = B + C;
    double err = sqrt(eB*eB + eC*eC);
    
    out->SetBinContent(b, val);
    out->SetBinError(b, err);
  }
  return out;
}
// ============================================================
// IO helpers
// ============================================================
static TH1D* LoadHistTryNames(const TString& folder, int rbin,
                              const std::vector<TString>& tryNames,
                              const TString& newName)
{
  TString fpath = Form("%s/fit_Rbin%d.root", folder.Data(), rbin);
  TFile* f = TFile::Open(fpath);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fpath << "\n";
    if (f) { f->Close(); delete f; }
    return nullptr;
  }

  TH1D* hFound = nullptr;
  for (const auto& hn : tryNames) {
    TH1D* h = dynamic_cast<TH1D*>(f->Get(hn));
    if (h) {
      hFound = h;
      break;
    }
  }

  if (!hFound) {
    std::cerr << "ERROR: cannot find any of [";
    for (size_t i = 0; i < tryNames.size(); ++i) {
      std::cerr << tryNames[i] << (i + 1 < tryNames.size() ? ", " : "");
    }
    std::cerr << "] in " << fpath << "\n";
    f->Close();
    delete f;
    return nullptr;
  }

  TH1D* out = (TH1D*)hFound->Clone(newName);
  out->SetDirectory(nullptr);

  f->Close();
  delete f;
  return out;
}
// SS
static TH1D* LoadNSS_SE_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldSS_vsCos" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hSigWin_BCD_vsCos" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_ME_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldSS_vsCos_ME" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hSigWin_BCD_vsCos_ME" }, newName);
  }
  return nullptr;
}
// SB
static TH1D* LoadNSS_SE_SB_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldSBsum_Model_vsCos" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_BCD_vsCos" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_ME_SB_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldSBsum_Model_vsCos_ME" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_BCD_vsCos_ME" }, newName);
  }
  return nullptr;
}
// BB
static TH1D* LoadNSS_SE_BB_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldBB_Model_vsCos" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_BCD_vsCos" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_ME_BB_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS){
    return LoadHistTryNames(folder, rbin, { "hYieldBB_Model_vsCos_ME" }, newName);
  }else if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_BCD_vsCos_ME" }, newName);
  }
  return nullptr;
}

//Add the B C D separately for the SE weitth sand band metthod
static TH1D* LoadNSS_SE_B_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_B_vsCos" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_SE_C_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_C_vsCos" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_SE_D_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_D_vsCos" }, newName);
  }
  return nullptr;
}
//Add the B C D separately for the ME weitth sand band metthod
static TH1D* LoadNSS_ME_B_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_B_vsCos_ME" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_ME_C_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_C_vsCos_ME" }, newName);
  }
  return nullptr;
}

static TH1D* LoadNSS_ME_D_singleFile(const TString& folder, int rbin, const TString& newName, bool do_YieldSS, bool do_YieldSS_BCD)
{
  if(do_YieldSS_BCD){
    return LoadHistTryNames(folder, rbin, { "hBkgWin_D_vsCos_ME" }, newName);
  }
  return nullptr;
}

//Load the  data in signal window

static TH1D* LoadData_SE__singleFile(const TString& folder, int rbin, const TString& newName)
{
  return LoadHistTryNames(folder, rbin, { "hDataWin_vsCos" }, newName);
}

static TH1D* LoadData_ME__singleFile(const TString& folder, int rbin, const TString& newName)
{
  return LoadHistTryNames(folder, rbin, { "hDataWin_vsCos_ME" }, newName);
}




// ============================================================
// MAIN
// ============================================================
// 1.efault options: ABCD method: signal:Double Gaussian ,background: Bern2,Bern3,Chebyshev2
// 2.fit method:combined fit method for signal and background:Double Gaussian ,background: Bern2,Bern3,Chebyshev2
void Result2(int centlow = 0, int centhigh = 100, int nRbins = 6, TString inputfold = "Output_pp",
             const double* rEdges = nullptr, int sigstat = 0, int bkgstat = 2, bool do_YieldSS =false, bool do_YieldSS_BCD = true)
{
  gStyle->SetOptStat(0);

  //const double defaultREdges[] = {0.0, 1.5, 3.1};
  //const double defaultREdges[] = {0.0, 3.1};
  const double defaultREdges[] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};
  
  const double* edges = rEdges ? rEdges : defaultREdges;

  // alpha values:
  // Lambda      : +0.746\pm 0.008
  // AntiLambda  : -0.758\pm 0.005

  // Lambda      : +0.747\pm 0.009
  // AntiLambda  : -0.757\pm 0.004

  const double scaleconst1[] = { 0.747,  0.747, -0.757, -0.757 };
  const double scaleconst2[] = { 0.747, -0.757,  0.747, -0.757 };

  std::vector<int> sigOpts = {sigstat};
  std::vector<int> bkgOpts = {bkgstat};

  const TString tag_LL   = "hSparseLambdaLambda";
  const TString tag_LAL  = "hSparseLambdaAntiLambda";
  const TString tag_ALL  = "hSparseAntiLambdaLambda";
  const TString tag_ALAL = "hSparseAntiLambdaAntiLambda";

  TString fitMethod = Form("%s", do_YieldSS ? "YieldSS" : "do_YieldSS_BCD");

  const TString outRootAll = Form("Results_%s_Rdep_%s_sig%d_bkg%d.root",
                                  fitMethod.Data(),inputfold.Data(), sigstat, bkgstat);

  std::unique_ptr<TFile> fAll(TFile::Open(outRootAll, "RECREATE"));
  if (!fAll || fAll->IsZombie()) {
    std::cerr << "ERROR: cannot create " << outRootAll << "\n";
    return;
  }

  const TString pngDir = Form("Results_%s_Rdep_pp_sig%d_bkg%d_PNGs_cent%dto%d_ALLCOMB", fitMethod.Data(), sigstat, bkgstat, centlow, centhigh);
  gSystem->mkdir(pngDir, kTRUE);

  for (int sig : sigOpts) {
    for (int bkg : bkgOpts) {

      const TString comb = Form("sig%d_bkg%d", sig, bkg);
      std::cout << "\n==================== " << comb << " ====================\n";

      const TString dir_LL   = Form("%s_sig%d_bkg%d_%s", inputfold.Data(), sig, bkg, tag_LL.Data());
      const TString dir_LAL  = Form("%s_sig%d_bkg%d_%s", inputfold.Data(), sig, bkg, tag_LAL.Data());
      const TString dir_ALL  = Form("%s_sig%d_bkg%d_%s", inputfold.Data(), sig, bkg, tag_ALL.Data());
      const TString dir_ALAL = Form("%s_sig%d_bkg%d_%s", inputfold.Data(), sig, bkg, tag_ALAL.Data());

      auto hP_LL_vsR   = std::make_unique<TH1D>(Form("hP_LL_vsR_%s",   comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_ALAL_vsR = std::make_unique<TH1D>(Form("hP_ALAL_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_SS_vsR   = std::make_unique<TH1D>(Form("hP_SS_vsR_%s",   comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_US_vsR   = std::make_unique<TH1D>(Form("hP_US_vsR_%s",   comb.Data()), "P vs R;R;P", nRbins, edges);

      auto hP_US_SB_vsR   = std::make_unique<TH1D>(Form("hP_US_SB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_SS_SB_vsR   = std::make_unique<TH1D>(Form("hP_SS_SB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);

      auto hP_US_BB_vsR   = std::make_unique<TH1D>(Form("hP_US_BB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_SS_BB_vsR   = std::make_unique<TH1D>(Form("hP_SS_BB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);

      auto hP_SS_BPlusCMinusD_vsR   = std::make_unique<TH1D>(Form("hP_SS_BPlusCMinusD_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_US_BPlusCMinusD_vsR   = std::make_unique<TH1D>(Form("hP_US_BPlusCMinusD_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);  
      
      auto hP_SS_BPlusC_vsR   = std::make_unique<TH1D>(Form("hP_SS_BPlusC_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_US_BPlusC_vsR   = std::make_unique<TH1D>(Form("hP_US_BPlusC_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);

      //
      auto hP_SS_SPlusB_vsR   = std::make_unique<TH1D>(Form("hP_SS_SPlusB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_US_SPlusB_vsR   = std::make_unique<TH1D>(Form("hP_US_SPlusB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      
      auto hP_LL_SPlusB_vsR   = std::make_unique<TH1D>(Form("hP_LL_SPlusB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);
      auto hP_ALAL_SPlusB_vsR = std::make_unique<TH1D>(Form("hP_ALAL_SPlusB_vsR_%s", comb.Data()), "P vs R;R;P", nRbins, edges);


      hP_LL_vsR->SetDirectory(nullptr);
      hP_ALAL_vsR->SetDirectory(nullptr);
      hP_SS_vsR->SetDirectory(nullptr);
      hP_US_vsR->SetDirectory(nullptr);
      hP_US_SB_vsR->SetDirectory(nullptr);
      hP_SS_SB_vsR->SetDirectory(nullptr);

      hP_US_BB_vsR->SetDirectory(nullptr);
      hP_SS_BB_vsR->SetDirectory(nullptr);

      hP_SS_BPlusCMinusD_vsR->SetDirectory(nullptr);
      hP_US_BPlusCMinusD_vsR->SetDirectory(nullptr);

      hP_SS_BPlusC_vsR->SetDirectory(nullptr);
      hP_US_BPlusC_vsR->SetDirectory(nullptr);

      hP_SS_SPlusB_vsR->SetDirectory(nullptr);
      hP_US_SPlusB_vsR->SetDirectory(nullptr);

      hP_LL_SPlusB_vsR->SetDirectory(nullptr);
      hP_ALAL_SPlusB_vsR->SetDirectory(nullptr);

      fAll->cd();
      TDirectory* dComb = fAll->mkdir(comb);
      if (!dComb) {
        std::cerr << "ERROR: cannot create directory " << comb << " in output file\n";
        continue;
      }
      dComb->cd();

      {
        TLatex meta;
        meta.SetName("meta");
        meta.SetTitle(Form("cent=%d-%d, %s", centlow, centhigh, comb.Data()));
        meta.Write();
      }

      for (int rbin = 0; rbin < nRbins; ++rbin) {

        // ---------------- Load SE and ME ----------------
        auto hLL_SE   = std::unique_ptr<TH1D>(LoadNSS_SE_singleFile(dir_LL,   rbin, Form("hLL_SE_%s_r%d",   comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLL_ME   = std::unique_ptr<TH1D>(LoadNSS_ME_singleFile(dir_LL,   rbin, Form("hLL_ME_%s_r%d",   comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));

        auto hLAL_SE  = std::unique_ptr<TH1D>(LoadNSS_SE_singleFile(dir_LAL,  rbin, Form("hLAL_SE_%s_r%d",  comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLAL_ME  = std::unique_ptr<TH1D>(LoadNSS_ME_singleFile(dir_LAL,  rbin, Form("hLAL_ME_%s_r%d",  comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));

        auto hALL_SE  = std::unique_ptr<TH1D>(LoadNSS_SE_singleFile(dir_ALL,  rbin, Form("hALL_SE_%s_r%d",  comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALL_ME  = std::unique_ptr<TH1D>(LoadNSS_ME_singleFile(dir_ALL,  rbin, Form("hALL_ME_%s_r%d",  comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));

        auto hALAL_SE = std::unique_ptr<TH1D>(LoadNSS_SE_singleFile(dir_ALAL, rbin, Form("hALAL_SE_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALAL_ME = std::unique_ptr<TH1D>(LoadNSS_ME_singleFile(dir_ALAL, rbin, Form("hALAL_ME_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));

      // SB-US
        auto hLAL_SE_SB = std::unique_ptr<TH1D>(LoadNSS_SE_SB_singleFile(dir_LAL, rbin, Form("hLAL_SE_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLAL_ME_SB = std::unique_ptr<TH1D>(LoadNSS_ME_SB_singleFile(dir_LAL, rbin, Form("hLAL_ME_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        
        auto hALL_SE_SB = std::unique_ptr<TH1D>(LoadNSS_SE_SB_singleFile(dir_ALL, rbin, Form("hALL_SE_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALL_ME_SB = std::unique_ptr<TH1D>(LoadNSS_ME_SB_singleFile(dir_ALL, rbin, Form("hALL_ME_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));  
       
      // SB-SS

        auto hLL_SE_SB = std::unique_ptr<TH1D>(LoadNSS_SE_SB_singleFile(dir_LL, rbin, Form("hLL_SE_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLL_ME_SB = std::unique_ptr<TH1D>(LoadNSS_ME_SB_singleFile(dir_LL, rbin, Form("hLL_ME_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        
        auto hALAL_SE_SB = std::unique_ptr<TH1D>(LoadNSS_SE_SB_singleFile(dir_ALAL, rbin, Form("hALAL_SE_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALAL_ME_SB = std::unique_ptr<TH1D>(LoadNSS_ME_SB_singleFile(dir_ALAL, rbin, Form("hALAL_ME_SB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));  
      
      // BB-US
        auto hLAL_SE_BB = std::unique_ptr<TH1D>(LoadNSS_SE_BB_singleFile(dir_LAL, rbin, Form("hLAL_SE_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLAL_ME_BB = std::unique_ptr<TH1D>(LoadNSS_ME_BB_singleFile(dir_LAL, rbin, Form("hLAL_ME_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        
        auto hALL_SE_BB = std::unique_ptr<TH1D>(LoadNSS_SE_BB_singleFile(dir_ALL, rbin, Form("hALL_SE_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALL_ME_BB = std::unique_ptr<TH1D>(LoadNSS_ME_BB_singleFile(dir_ALL, rbin, Form("hALL_ME_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));  
      // BB-SS
        auto hLL_SE_BB = std::unique_ptr<TH1D>(LoadNSS_SE_BB_singleFile(dir_LL, rbin, Form("hLL_SE_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hLL_ME_BB = std::unique_ptr<TH1D>(LoadNSS_ME_BB_singleFile(dir_LL, rbin, Form("hLL_ME_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        
        auto hALAL_SE_BB = std::unique_ptr<TH1D>(LoadNSS_SE_BB_singleFile(dir_ALAL, rbin, Form("hALAL_SE_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        auto hALAL_ME_BB = std::unique_ptr<TH1D>(LoadNSS_ME_BB_singleFile(dir_ALAL, rbin, Form("hALAL_ME_BB_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));  
        // Add the B C D separately for the SE in sand band method

        std::unique_ptr<TH1D> hLL_SE_B = nullptr, hLL_SE_C = nullptr, hLL_SE_D = nullptr;
        std::unique_ptr<TH1D> hALAL_SE_B = nullptr, hALAL_SE_C = nullptr, hALAL_SE_D = nullptr;
        std::unique_ptr<TH1D> hLAL_SE_B = nullptr, hLAL_SE_C = nullptr, hLAL_SE_D = nullptr;
        std::unique_ptr<TH1D> hALL_SE_B = nullptr, hALL_SE_C = nullptr, hALL_SE_D = nullptr;

        std::unique_ptr<TH1D> hLL_ME_B = nullptr, hLL_ME_C = nullptr, hLL_ME_D = nullptr;
        std::unique_ptr<TH1D> hALAL_ME_B = nullptr, hALAL_ME_C = nullptr, hALAL_ME_D = nullptr;
        std::unique_ptr<TH1D> hLAL_ME_B = nullptr, hLAL_ME_C = nullptr, hLAL_ME_D = nullptr;
        std::unique_ptr<TH1D> hALL_ME_B = nullptr, hALL_ME_C = nullptr, hALL_ME_D = nullptr;
        if (do_YieldSS_BCD){
          //Like-sign
          hLL_SE_B = std::unique_ptr<TH1D>(LoadNSS_SE_B_singleFile(dir_LL, rbin, Form("hLL_SE_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLL_SE_C = std::unique_ptr<TH1D>(LoadNSS_SE_C_singleFile(dir_LL, rbin, Form("hLL_SE_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLL_SE_D = std::unique_ptr<TH1D>(LoadNSS_SE_D_singleFile(dir_LL, rbin, Form("hLL_SE_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_SE_B = std::unique_ptr<TH1D>(LoadNSS_SE_B_singleFile(dir_ALAL, rbin, Form("hALAL_SE_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_SE_C = std::unique_ptr<TH1D>(LoadNSS_SE_C_singleFile(dir_ALAL, rbin, Form("hALAL_SE_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_SE_D = std::unique_ptr<TH1D>(LoadNSS_SE_D_singleFile(dir_ALAL, rbin, Form("hALAL_SE_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          //Unlike-sign
          hLAL_SE_B = std::unique_ptr<TH1D>(LoadNSS_SE_B_singleFile(dir_LAL, rbin, Form("hLAL_SE_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLAL_SE_C = std::unique_ptr<TH1D>(LoadNSS_SE_C_singleFile(dir_LAL, rbin, Form("hLAL_SE_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLAL_SE_D = std::unique_ptr<TH1D>(LoadNSS_SE_D_singleFile(dir_LAL, rbin, Form("hLAL_SE_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_SE_B = std::unique_ptr<TH1D>(LoadNSS_SE_B_singleFile(dir_ALL, rbin, Form("hALL_SE_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_SE_C = std::unique_ptr<TH1D>(LoadNSS_SE_C_singleFile(dir_ALL, rbin, Form("hALL_SE_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_SE_D = std::unique_ptr<TH1D>(LoadNSS_SE_D_singleFile(dir_ALL, rbin, Form("hALL_SE_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD)); 


          // Add the B C D separately for the ME in sand band method
          //Like-sign
          hLL_ME_B = std::unique_ptr<TH1D>(LoadNSS_ME_B_singleFile(dir_LL, rbin, Form("hLL_ME_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLL_ME_C = std::unique_ptr<TH1D>(LoadNSS_ME_C_singleFile(dir_LL, rbin, Form("hLL_ME_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLL_ME_D = std::unique_ptr<TH1D>(LoadNSS_ME_D_singleFile(dir_LL, rbin, Form("hLL_ME_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_ME_B = std::unique_ptr<TH1D>(LoadNSS_ME_B_singleFile(dir_ALAL, rbin, Form("hALAL_ME_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_ME_C = std::unique_ptr<TH1D>(LoadNSS_ME_C_singleFile(dir_ALAL, rbin, Form("hALAL_ME_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALAL_ME_D = std::unique_ptr<TH1D>(LoadNSS_ME_D_singleFile(dir_ALAL, rbin, Form("hALAL_ME_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          //Unlike-sign
          hLAL_ME_B = std::unique_ptr<TH1D>(LoadNSS_ME_B_singleFile(dir_LAL, rbin, Form("hLAL_ME_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLAL_ME_C = std::unique_ptr<TH1D>(LoadNSS_ME_C_singleFile(dir_LAL, rbin, Form("hLAL_ME_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hLAL_ME_D = std::unique_ptr<TH1D>(LoadNSS_ME_D_singleFile(dir_LAL, rbin, Form("hLAL_ME_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_ME_B = std::unique_ptr<TH1D>(LoadNSS_ME_B_singleFile(dir_ALL, rbin, Form("hALL_ME_B_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_ME_C = std::unique_ptr<TH1D>(LoadNSS_ME_C_singleFile(dir_ALL, rbin, Form("hALL_ME_C_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
          hALL_ME_D = std::unique_ptr<TH1D>(LoadNSS_ME_D_singleFile(dir_ALL, rbin, Form("hALL_ME_D_%s_r%d", comb.Data(), rbin),do_YieldSS, do_YieldSS_BCD));
        }

        //extract the P_S+B
        //SE
        auto hLLSPlusB_SE = std::unique_ptr<TH1D>(LoadData_SE__singleFile(dir_LL, rbin, Form("hLLSPlusB_SE_%s_r%d", comb.Data(), rbin)));
        auto hLALSPlusB_SE = std::unique_ptr<TH1D>(LoadData_SE__singleFile(dir_LAL, rbin, Form("hLALSPlusB_SE_%s_r%d", comb.Data(), rbin)));
        auto hALLSPlusB_SE = std::unique_ptr<TH1D>(LoadData_SE__singleFile(dir_ALL, rbin, Form("hALLSPlusB_SE_%s_r%d", comb.Data(), rbin)));
        auto hALALSPlusB_SE = std::unique_ptr<TH1D>(LoadData_SE__singleFile(dir_ALAL, rbin, Form("hALALSPlusB_SE_%s_r%d", comb.Data(), rbin)));
        //ME
        auto hLLSPlusB_ME = std::unique_ptr<TH1D>(LoadData_ME__singleFile(dir_LL, rbin, Form("hLLSPlusB_ME_%s_r%d", comb.Data(), rbin)));
        auto hLALSPlusB_ME = std::unique_ptr<TH1D>(LoadData_ME__singleFile(dir_LAL, rbin, Form("hLALSPlusB_ME_%s_r%d", comb.Data(), rbin)));
        auto hALLSPlusB_ME = std::unique_ptr<TH1D>(LoadData_ME__singleFile(dir_ALL, rbin, Form("hALLSPlusB_ME_%s_r%d", comb.Data(), rbin)));
        auto hALALSPlusB_ME = std::unique_ptr<TH1D>(LoadData_ME__singleFile(dir_ALAL, rbin, Form("hALALSPlusB_ME_%s_r%d", comb.Data(), rbin)));

        // ---------------- Acceptance corrected individual same-sign ----------------
        
        auto hCorr_LL = std::unique_ptr<TH1D>(
          (hLL_SE && hLL_ME)
            ? MakeAcceptanceCorrected(hLL_SE.get(), hLL_ME.get(), Form("hCorr_LL_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hCorr_ALAL = std::unique_ptr<TH1D>(
          (hALAL_SE && hALAL_ME)
            ? MakeAcceptanceCorrected(hALAL_SE.get(), hALAL_ME.get(), Form("hCorr_ALAL_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        // ---------------- Same-sign average: sum SE and sum ME first, then divide once ----------------
        auto hSE_SS = std::unique_ptr<TH1D>(
          (hLL_SE && hALAL_SE)
            ? MakeSummedHist(hLL_SE.get(), hALAL_SE.get(), Form("hSE_SS_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hME_SS = std::unique_ptr<TH1D>(
          (hLL_ME && hALAL_ME)
            ? MakeSummedHist(hLL_ME.get(), hALAL_ME.get(), Form("hME_SS_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hCorr_SS = std::unique_ptr<TH1D>(
          (hSE_SS && hME_SS)
            ? MakeAcceptanceCorrected(hSE_SS.get(), hME_SS.get(), Form("hCorr_SS_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        // SB
        auto hSE_SS_SB = std::unique_ptr<TH1D>(
          (hLL_SE_SB && hALAL_SE_SB)
            ? MakeSummedHist(hLL_SE_SB.get(), hALAL_SE_SB.get(), Form("hSE_SS_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_SS_SB = std::unique_ptr<TH1D>(
          (hLL_ME_SB && hALAL_ME_SB)
            ? MakeSummedHist(hLL_ME_SB.get(), hALAL_ME_SB.get(), Form("hME_SS_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_SS_SB = std::unique_ptr<TH1D>(
          (hSE_SS_SB && hME_SS_SB)
            ? MakeAcceptanceCorrected(hSE_SS_SB.get(), hME_SS_SB.get(), Form("hCorr_SS_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        // BB
        auto hSE_SS_BB = std::unique_ptr<TH1D>(
          (hLL_SE_BB && hALAL_SE_BB)
            ? MakeSummedHist(hLL_SE_BB.get(), hALAL_SE_BB.get(), Form("hSE_SS_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_SS_BB = std::unique_ptr<TH1D>(
          (hLL_ME_BB && hALAL_ME_BB)
            ? MakeSummedHist(hLL_ME_BB.get(), hALAL_ME_BB.get(), Form("hME_SS_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_SS_BB = std::unique_ptr<TH1D>(
          (hSE_SS_BB && hME_SS_BB)
            ? MakeAcceptanceCorrected(hSE_SS_BB.get(), hME_SS_BB.get(), Form("hCorr_SS_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        //B C and D separately for the SS///////////////////////////////////////////////////////////////
        auto hSE_SS_B = std::unique_ptr<TH1D>(
          (hLL_SE_B && hALAL_SE_B)
            ? MakeSummedHist(hLL_SE_B.get(), hALAL_SE_B.get(), Form("hSE_SS_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_SS_B = std::unique_ptr<TH1D>(
          (hLL_ME_B && hALAL_ME_B)
            ? MakeSummedHist(hLL_ME_B.get(), hALAL_ME_B.get(), Form("hME_SS_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_SS_B = std::unique_ptr<TH1D>(
          (hSE_SS_B && hME_SS_B)
            ? MakeAcceptanceCorrected(hSE_SS_B.get(), hME_SS_B.get(), Form("hCorr_SS_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hSE_SS_C = std::unique_ptr<TH1D>(
          (hLL_SE_C && hALAL_SE_C)
            ? MakeSummedHist(hLL_SE_C.get(), hALAL_SE_C.get(), Form("hSE_SS_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_SS_C = std::unique_ptr<TH1D>(
          (hLL_ME_C && hALAL_ME_C)
            ? MakeSummedHist(hLL_ME_C.get(), hALAL_ME_C.get(), Form("hME_SS_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_SS_C = std::unique_ptr<TH1D>(
          (hSE_SS_C && hME_SS_C)
            ? MakeAcceptanceCorrected(hSE_SS_C.get(), hME_SS_C.get(), Form("hCorr_SS_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hSE_SS_D = std::unique_ptr<TH1D>(
          (hLL_SE_D && hALAL_SE_D)
            ? MakeSummedHist(hLL_SE_D.get(), hALAL_SE_D.get(), Form("hSE_SS_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_SS_D = std::unique_ptr<TH1D>(
          (hLL_ME_D && hALAL_ME_D)
            ? MakeSummedHist(hLL_ME_D.get(), hALAL_ME_D.get(), Form("hME_SS_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_SS_D = std::unique_ptr<TH1D>(
          (hSE_SS_D && hME_SS_D)
            ? MakeAcceptanceCorrected(hSE_SS_D.get(), hME_SS_D.get(), Form("hCorr_SS_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );


        // ---------------- Unlike-sign average: sum SE and sum ME first, then divide once ----------------
        auto hSE_US = std::unique_ptr<TH1D>(
          (hLAL_SE && hALL_SE)
            ? MakeSummedHist(hLAL_SE.get(), hALL_SE.get(), Form("hSE_US_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hME_US = std::unique_ptr<TH1D>(
          (hLAL_ME && hALL_ME)
            ? MakeSummedHist(hLAL_ME.get(), hALL_ME.get(), Form("hME_US_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hCorr_US = std::unique_ptr<TH1D>(
          (hSE_US && hME_US)
            ? MakeAcceptanceCorrected(hSE_US.get(), hME_US.get(), Form("hCorr_US_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        // SB
        auto hSE_US_SB = std::unique_ptr<TH1D>(
          (hLAL_SE_SB && hALL_SE_SB)
            ? MakeSummedHist(hLAL_SE_SB.get(), hALL_SE_SB.get(), Form("hSE_US_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_US_SB = std::unique_ptr<TH1D>(
          (hLAL_ME_SB && hALL_ME_SB)
            ? MakeSummedHist(hLAL_ME_SB.get(), hALL_ME_SB.get(), Form("hME_US_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_SB = std::unique_ptr<TH1D>(
          (hSE_US_SB && hME_US_SB)
            ? MakeAcceptanceCorrected(hSE_US_SB.get(), hME_US_SB.get(), Form("hCorr_US_SB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        // BB
        auto hSE_US_BB = std::unique_ptr<TH1D>(
          (hLAL_SE_BB && hALL_SE_BB)
            ? MakeSummedHist(hLAL_SE_BB.get(), hALL_SE_BB.get(), Form("hSE_US_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_US_BB = std::unique_ptr<TH1D>(
          (hLAL_ME_BB && hALL_ME_BB)
            ? MakeSummedHist(hLAL_ME_BB.get(), hALL_ME_BB.get(), Form("hME_US_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_BB = std::unique_ptr<TH1D>(
          (hSE_US_BB && hME_US_BB)
            ? MakeAcceptanceCorrected(hSE_US_BB.get(), hME_US_BB.get(), Form("hCorr_US_BB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        //B C and D separately for the US
        auto hSE_US_B = std::unique_ptr<TH1D>(
          (hLAL_SE_B && hALL_SE_B)
            ? MakeSummedHist(hLAL_SE_B.get(), hALL_SE_B.get(), Form("hSE_US_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_US_B = std::unique_ptr<TH1D>(
          (hLAL_ME_B && hALL_ME_B)
            ? MakeSummedHist(hLAL_ME_B.get(), hALL_ME_B.get(), Form("hME_US_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_B = std::unique_ptr<TH1D>(
          (hSE_US_B && hME_US_B)
            ? MakeAcceptanceCorrected(hSE_US_B.get(), hME_US_B.get(), Form("hCorr_US_B_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hSE_US_C = std::unique_ptr<TH1D>(
          (hLAL_SE_C && hALL_SE_C)
            ? MakeSummedHist(hLAL_SE_C.get(), hALL_SE_C.get(), Form("hSE_US_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_US_C = std::unique_ptr<TH1D>(
          (hLAL_ME_C && hALL_ME_C)
            ? MakeSummedHist(hLAL_ME_C.get(), hALL_ME_C.get(), Form("hME_US_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_C = std::unique_ptr<TH1D>(
          (hSE_US_C && hME_US_C)
            ? MakeAcceptanceCorrected(hSE_US_C.get(), hME_US_C.get(), Form("hCorr_US_C_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hSE_US_D = std::unique_ptr<TH1D>(
          (hLAL_SE_D && hALL_SE_D)
            ? MakeSummedHist(hLAL_SE_D.get(), hALL_SE_D.get(), Form("hSE_US_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hME_US_D = std::unique_ptr<TH1D>(
          (hLAL_ME_D && hALL_ME_D)
            ? MakeSummedHist(hLAL_ME_D.get(), hALL_ME_D.get(), Form("hME_US_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_D = std::unique_ptr<TH1D>(
          (hSE_US_D && hME_US_D)
            ? MakeAcceptanceCorrected(hSE_US_D.get(), hME_US_D.get(), Form("hCorr_US_D_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        // make B+C-D separately for the SS and US //////////////////////////////////////////////
        //Like-sign 
        auto hCorr_SS_BPlusCMinusD = std::unique_ptr<TH1D>(
          (hCorr_SS_B && hCorr_SS_C && hCorr_SS_D)
            ? MakeBPlusCMinusDHist(hCorr_SS_B.get(), hCorr_SS_C.get(), hCorr_SS_D.get(), Form("hCorr_SS_BPlusCMinusD_%s_rbin%d", comb.Data(), rbin))  
            : nullptr
        );
        //Unlike-sign
        auto hCorr_US_BPlusCMinusD = std::unique_ptr<TH1D>(
          (hCorr_US_B && hCorr_US_C && hCorr_US_D)
            ? MakeBPlusCMinusDHist(hCorr_US_B.get(), hCorr_US_C.get(), hCorr_US_D.get(), Form("hCorr_US_BPlusCMinusD_%s_rbin%d", comb.Data(), rbin))  
            : nullptr
        );  

        //MakeBPlusCHist
        auto hCorr_SS_BPlusC = std::unique_ptr<TH1D>(
          (hCorr_SS_B && hCorr_SS_C)
            ? MakeBPlusCHist(hCorr_SS_B.get(), hCorr_SS_C.get(), Form("hCorr_SS_BPlusC_%s_rbin%d", comb.Data(), rbin))  
            : nullptr
        );
        auto hCorr_US_BPlusC = std::unique_ptr<TH1D>(
          (hCorr_US_B && hCorr_US_C)
            ? MakeBPlusCHist(hCorr_US_B.get(), hCorr_US_C.get(), Form("hCorr_US_BPlusC_%s_rbin%d", comb.Data(), rbin))  
            : nullptr
        );

        //S Plus B in signal region for the SE and ME separately
        auto hSS_SPlusB_SE = std::unique_ptr<TH1D>(
          (hLLSPlusB_SE && hALALSPlusB_SE)
            ? MakeSummedHist(hLLSPlusB_SE.get(), hALALSPlusB_SE.get(), Form("hSS_SPlusB_SE_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hSS_SPlusB_ME = std::unique_ptr<TH1D>(
          (hLLSPlusB_ME && hALALSPlusB_ME)
            ? MakeSummedHist(hLLSPlusB_ME.get(), hALALSPlusB_ME.get(), Form("hSS_SPlusB_ME_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hUS_SPlusB_SE = std::unique_ptr<TH1D>(
          (hLALSPlusB_SE && hALLSPlusB_SE)
            ? MakeSummedHist(hLALSPlusB_SE.get(), hALLSPlusB_SE.get(), Form("hUS_SPlusB_SE_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hUS_SPlusB_ME = std::unique_ptr<TH1D>(
          (hLALSPlusB_ME && hALLSPlusB_ME)
            ? MakeSummedHist(hLALSPlusB_ME.get(), hALLSPlusB_ME.get(), Form("hUS_SPlusB_ME_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        // do correction for the S+B in the signal region for the SS and US
        auto hCorr_SS_SPlusB = std::unique_ptr<TH1D>(
          (hSS_SPlusB_SE && hSS_SPlusB_ME)
            ? MakeAcceptanceCorrected(hSS_SPlusB_SE.get(), hSS_SPlusB_ME.get(), Form("hCorr_SS_SPlusB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );
        auto hCorr_US_SPlusB = std::unique_ptr<TH1D>(
          (hUS_SPlusB_SE && hUS_SPlusB_ME)
            ? MakeAcceptanceCorrected(hUS_SPlusB_SE.get(), hUS_SPlusB_ME.get(), Form("hCorr_US_SPlusB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );

        auto hCorr_LL_SPlusB = std::unique_ptr<TH1D>(
          (hLLSPlusB_SE && hLLSPlusB_ME)
            ? MakeAcceptanceCorrected(hLLSPlusB_SE.get(), hLLSPlusB_ME.get(), Form("hCorr_LL_SPlusB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );  

        auto hCorr_ALAL_SPlusB = std::unique_ptr<TH1D>(
          (hALALSPlusB_SE && hALALSPlusB_ME)
            ? MakeAcceptanceCorrected(hALALSPlusB_SE.get(), hALALSPlusB_ME.get(), Form("hCorr_ALAL_SPlusB_%s_rbin%d", comb.Data(), rbin))
            : nullptr
        );







        // ---------------- Fit -> P extraction ----------------
        FitOut fLL, fALAL, fSS, fUS, fUS_SB,fSS_SB,fUS_BB,fSS_BB,fSS_BPlusCMinusD,fUS_BPlusCMinusD, fSS_BPlusC, fUS_BPlusC,
        fSS_SPlusB, fUS_SPlusB,fLL_SPlusB, fALAL_SPlusB;

        const double C_LL   = scaleconst1[0] * scaleconst2[0];
        const double C_ALAL = scaleconst1[3] * scaleconst2[3];
        const double C_LAL  = scaleconst1[1] * scaleconst2[1];
        const double C_ALL  = scaleconst1[2] * scaleconst2[2];
        const double C_SS   = 0.5 * (C_LL + C_ALAL);
        const double C_US   = 0.5 * (C_LAL + C_ALL);

        if (hCorr_LL) {
          fLL = FitAndReturn(hCorr_LL.get(), C_LL, Form("fit_LL_%s_rbin%d", comb.Data(), rbin));
          hP_LL_vsR->SetBinContent(rbin + 1, fLL.P);
          hP_LL_vsR->SetBinError  (rbin + 1, fLL.Perr);
        }

        if (hCorr_ALAL) {
          fALAL = FitAndReturn(hCorr_ALAL.get(), C_ALAL, Form("fit_ALAL_%s_rbin%d", comb.Data(), rbin));
          hP_ALAL_vsR->SetBinContent(rbin + 1, fALAL.P);
          hP_ALAL_vsR->SetBinError  (rbin + 1, fALAL.Perr);
        }

        if (hCorr_SS) {
          fSS = FitAndReturn(hCorr_SS.get(), C_SS, Form("fit_SS_%s_rbin%d", comb.Data(), rbin));
          hP_SS_vsR->SetBinContent(rbin + 1, fSS.P);
          hP_SS_vsR->SetBinError  (rbin + 1, fSS.Perr);
        }

        if (hCorr_US) {
          fUS = FitAndReturn(hCorr_US.get(), C_US, Form("fit_US_%s_rbin%d", comb.Data(), rbin));
          hP_US_vsR->SetBinContent(rbin + 1, fUS.P);
          hP_US_vsR->SetBinError  (rbin + 1, fUS.Perr);
        }
        // SB
        if (hCorr_US_SB) {
          fUS_SB = FitAndReturn(hCorr_US_SB.get(), 1.0, Form("fit_US_SB_%s_rbin%d", comb.Data(), rbin));
          hP_US_SB_vsR->SetBinContent(rbin + 1, fUS_SB.P);
          hP_US_SB_vsR->SetBinError  (rbin + 1, fUS_SB.Perr);
        }
        if (hCorr_SS_SB) {
          fSS_SB = FitAndReturn(hCorr_SS_SB.get(), 1.0, Form("fit_SS_SB_%s_rbin%d", comb.Data(), rbin));
          hP_SS_SB_vsR->SetBinContent(rbin + 1, fSS_SB.P);
          hP_SS_SB_vsR->SetBinError  (rbin + 1, fSS_SB.Perr);
        }
        // BB
        if (hCorr_US_BB) {
          fUS_BB = FitAndReturn(hCorr_US_BB.get(), 1.0, Form("fit_US_BB_%s_rbin%d", comb.Data(), rbin));
          hP_US_BB_vsR->SetBinContent(rbin + 1, fUS_BB.P);
          hP_US_BB_vsR->SetBinError(rbin + 1, fUS_BB.Perr);
        }
        if (hCorr_SS_BB) {
          fSS_BB = FitAndReturn(hCorr_SS_BB.get(), 1.0, Form("fit_SS_BB_%s_rbin%d", comb.Data(), rbin));
          hP_SS_BB_vsR->SetBinContent(rbin + 1, fSS_BB.P);
          hP_SS_BB_vsR->SetBinError(rbin + 1, fSS_BB.Perr);
        }

          // B C D separately for the SS and US
          //SS
        if (hCorr_SS_BPlusCMinusD) {
          fSS_BPlusCMinusD = FitAndReturn(hCorr_SS_BPlusCMinusD.get(), 1.0, Form("fit_SS_BPlusCMinusD_%s_rbin%d", comb.Data(), rbin));
          hP_SS_BPlusCMinusD_vsR->SetBinContent(rbin + 1, fSS_BPlusCMinusD.P);
          hP_SS_BPlusCMinusD_vsR->SetBinError(rbin + 1, fSS_BPlusCMinusD.Perr);
        }
        //US
        if (hCorr_US_BPlusCMinusD) {
          fUS_BPlusCMinusD = FitAndReturn(hCorr_US_BPlusCMinusD.get(), 1.0, Form("fit_US_BPlusCMinusD_%s_rbin%d", comb.Data(), rbin));
          hP_US_BPlusCMinusD_vsR->SetBinContent(rbin + 1, fUS_BPlusCMinusD.P);
          hP_US_BPlusCMinusD_vsR->SetBinError  (rbin + 1, fUS_BPlusCMinusD.Perr);
        }
        // B+C
        if (hCorr_SS_BPlusC) {
          fSS_BPlusC = FitAndReturn(hCorr_SS_BPlusC.get(), 1.0, Form("fit_SS_BPlusC_%s_rbin%d", comb.Data(), rbin));
          hP_SS_BPlusC_vsR->SetBinContent(rbin + 1, fSS_BPlusC.P);
          hP_SS_BPlusC_vsR->SetBinError(rbin + 1, fSS_BPlusC.Perr);
        }
        if (hCorr_US_BPlusC) {
          fUS_BPlusC = FitAndReturn(hCorr_US_BPlusC.get(), 1.0, Form("fit_US_BPlusC_%s_rbin%d", comb.Data(), rbin));
          hP_US_BPlusC_vsR->SetBinContent(rbin + 1, fUS_BPlusC.P);
          hP_US_BPlusC_vsR->SetBinError(rbin + 1, fUS_BPlusC.Perr);
        }

        // S+B in the signal region for the SS and US
        if (hCorr_SS_SPlusB) {
          fSS_SPlusB = FitAndReturn(hCorr_SS_SPlusB.get(), C_SS, Form("fit_SS_SPlusB_%s_rbin%d", comb.Data(), rbin));
          hP_SS_SPlusB_vsR->SetBinContent(rbin + 1, fSS_SPlusB.P);
          hP_SS_SPlusB_vsR->SetBinError(rbin + 1, fSS_SPlusB.Perr);
        }
        if (hCorr_US_SPlusB) {
          fUS_SPlusB = FitAndReturn(hCorr_US_SPlusB.get(), C_US, Form("fit_US_SPlusB_%s_rbin%d", comb.Data(), rbin));
          hP_US_SPlusB_vsR->SetBinContent(rbin + 1, fUS_SPlusB.P);
          hP_US_SPlusB_vsR->SetBinError(rbin + 1, fUS_SPlusB.Perr);
        }
        //LL and ALAL separately for the S+B in the signal region
        if (hCorr_LL_SPlusB) {
          fLL_SPlusB = FitAndReturn(hCorr_LL_SPlusB.get(), C_LL, Form("fit_LL_SPlusB_%s_rbin%d", comb.Data(), rbin));
          hP_LL_SPlusB_vsR->SetBinContent(rbin + 1, fLL_SPlusB.P);
          hP_LL_SPlusB_vsR->SetBinError(rbin + 1, fLL_SPlusB.Perr);
        }
        if (hCorr_ALAL_SPlusB) {
          fALAL_SPlusB = FitAndReturn(hCorr_ALAL_SPlusB.get(), C_ALAL, Form("fit_ALAL_SPlusB_%s_rbin%d", comb.Data(), rbin));
          hP_ALAL_SPlusB_vsR->SetBinContent(rbin + 1, fALAL_SPlusB.P);
          hP_ALAL_SPlusB_vsR->SetBinError(rbin + 1, fALAL_SPlusB.Perr);
        }
      
        // ---------------- Common y-range for final 2-panel plot ----------------
        double ymin = +1e9, ymax = -1e9;
        UpdateRangeFrom(hCorr_SS.get(), ymin, ymax);
        UpdateRangeFrom(hCorr_US.get(), ymin, ymax);

        if (!(ymin < ymax)) { ymin = 0.5; ymax = 1.5; }
        const double dy = std::max(1e-6, ymax - ymin);
        ymin -= 0.12 * dy;
        ymax += 0.20 * dy;

        // ---------------- Final 2-panel plot ----------------
        auto c = std::make_unique<TCanvas>(
          Form("c2panel_%s_cent%dto%d_rbin%d", comb.Data(), centlow, centhigh, rbin),
          "2-panel corrected angular distributions", 1000, 500
        );
        c->Divide(2, 1, 0.001, 0.001);

        c->cd(1);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04);
        DrawPanel(hCorr_SS.get(), fSS, Form("Same-sign average (%s, Rbin %d)", comb.Data(), rbin), ymin, ymax);

        c->cd(2);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04);
        DrawPanel(hCorr_US.get(), fUS, Form("Unlike-sign average (%s, Rbin %d)", comb.Data(), rbin), ymin, ymax);

        const TString png2 = Form("%s/AngDist_2panel_%s_cent%dto%d_rbin%d.png",
                                  pngDir.Data(), comb.Data(), centlow, centhigh, rbin);

        const TString pdf2 = Form("%s/AngDist_2panel_%s_cent%dto%d_rbin%d.pdf",
                                  pngDir.Data(), comb.Data(), centlow, centhigh, rbin);
                                  
        c->SaveAs(png2);
        c->SaveAs(pdf2);

        // ---------------- Write per-rbin objects ----------------
        dComb->cd();
        TDirectory* dR = dComb->mkdir(Form("rbin%d", rbin));
        if (dR) {
          dR->cd();

          // keep all existing info
          if (hLL_SE)      hLL_SE->Write();
          if (hLL_ME)      hLL_ME->Write();
          if (hALAL_SE)    hALAL_SE->Write();
          if (hALAL_ME)    hALAL_ME->Write();
          if (hSE_US)      hSE_US->Write();
          if (hME_US)      hME_US->Write();

          // SB - US
          if (hLAL_SE_SB) hLAL_SE_SB->Write();
          if (hLAL_ME_SB) hLAL_ME_SB->Write();
          if (hALL_SE_SB) hALL_SE_SB->Write();
          if (hALL_ME_SB) hALL_ME_SB->Write();
          if (hCorr_US_SB) hCorr_US_SB->Write();
          if (fUS_SB.f)     fUS_SB.f->Write();
          // BB - US
          if (hLAL_SE_BB) hLAL_SE_BB->Write();
          if (hLAL_ME_BB) hLAL_ME_BB->Write();
          if (hALL_SE_BB) hALL_SE_BB->Write();
          if (hALL_ME_BB) hALL_ME_BB->Write();
          if (hCorr_US_BB) hCorr_US_BB->Write();
          if (fUS_BB.f)     fUS_BB.f->Write();
          // SB - SS
          if (hLL_SE_SB) hLL_SE_SB->Write();
          if (hLL_ME_SB) hLL_ME_SB->Write();
          if (hALAL_SE_SB) hALAL_SE_SB->Write();
          if (hALAL_ME_SB) hALAL_ME_SB->Write();
          if (hCorr_SS_SB) hCorr_SS_SB->Write();
          if (fSS_SB.f)     fSS_SB.f->Write();
          // BB - SS
          if (hLL_SE_BB) hLL_SE_BB->Write();
          if (hLL_ME_BB) hLL_ME_BB->Write();  
          if (hALAL_SE_BB) hALAL_SE_BB->Write();
          if (hALAL_ME_BB) hALAL_ME_BB->Write();
          if (hCorr_SS_BB) hCorr_SS_BB->Write();
          if (fSS_BB.f)     fSS_BB.f->Write();

          // S+B in the signal region for the SS and US
          if(hLLSPlusB_SE) hLLSPlusB_SE->Write();
          if(hALALSPlusB_SE) hALALSPlusB_SE->Write();
          if(hLALSPlusB_SE) hLALSPlusB_SE->Write();
          if(hALLSPlusB_SE) hALLSPlusB_SE->Write();

          if(hLLSPlusB_ME) hLLSPlusB_ME->Write();
          if(hALALSPlusB_ME) hALALSPlusB_ME->Write();
          if(hLALSPlusB_ME) hLALSPlusB_ME->Write();
          if(hALLSPlusB_ME) hALLSPlusB_ME->Write();

          if(hSS_SPlusB_SE) hSS_SPlusB_SE->Write();
          if(hSS_SPlusB_ME) hSS_SPlusB_ME->Write();

          if(hUS_SPlusB_SE) hUS_SPlusB_SE->Write();
          if(hUS_SPlusB_ME) hUS_SPlusB_ME->Write();
            
          if(hCorr_US_SPlusB) hCorr_US_SPlusB->Write();
          if(hCorr_SS_SPlusB) hCorr_SS_SPlusB->Write();


          if(fSS_SPlusB.f)fSS_SPlusB.f->Write();
          if(fUS_SPlusB.f)fUS_SPlusB.f->Write();

          if(hCorr_LL_SPlusB) hCorr_LL_SPlusB->Write();
          if(hCorr_ALAL_SPlusB) hCorr_ALAL_SPlusB->Write();

       


          // BB - SS B,C,D
          if (do_YieldSS_BCD){
            // SE
            if (hLL_SE_B) hLL_SE_B->Write();
            if (hLL_SE_C) hLL_SE_C->Write();
            if (hLL_SE_D) hLL_SE_D->Write();
            if (hALAL_SE_B) hALAL_SE_B->Write();
            if (hALAL_SE_C) hALAL_SE_C->Write();
            if (hALAL_SE_D) hALAL_SE_D->Write();
            if (hLAL_SE_B) hLAL_SE_B->Write();
            if (hLAL_SE_C) hLAL_SE_C->Write();
            if (hLAL_SE_D) hLAL_SE_D->Write();
            if (hALL_SE_B) hALL_SE_B->Write();
            if (hALL_SE_C) hALL_SE_C->Write();
            if (hALL_SE_D) hALL_SE_D->Write();
            //ME
            if (hLL_ME_B) hLL_ME_B->Write();
            if (hLL_ME_C) hLL_ME_C->Write();
            if (hLL_ME_D) hLL_ME_D->Write();
            if (hALAL_ME_B) hALAL_ME_B->Write();
            if (hALAL_ME_C) hALAL_ME_C->Write();
            if (hALAL_ME_D) hALAL_ME_D->Write();
            if (hLAL_ME_B) hLAL_ME_B->Write();
            if (hLAL_ME_C) hLAL_ME_C->Write();
            if (hLAL_ME_D) hLAL_ME_D->Write();
            if (hALL_ME_B) hALL_ME_B->Write();
            if (hALL_ME_C) hALL_ME_C->Write();
            if (hALL_ME_D) hALL_ME_D->Write();

            if(hSE_SS_B)  hSE_SS_B->Write();
            if(hME_SS_B)  hME_SS_B->Write();
            if(hCorr_SS_B)  hCorr_SS_B->Write();
            if(hSE_SS_C)  hSE_SS_C->Write();
            if(hME_SS_C)  hME_SS_C->Write();
            if(hCorr_SS_C)  hCorr_SS_C->Write();
            if(hSE_SS_D)  hSE_SS_D->Write();
            if(hME_SS_D)  hME_SS_D->Write();
            if(hCorr_SS_D)  hCorr_SS_D->Write();  
            
            if(hSE_US_B)  hSE_US_B->Write();
            if(hME_US_B)  hME_US_B->Write();
            if(hCorr_US_B)  hCorr_US_B->Write();
            if(hSE_US_C)  hSE_US_C->Write();
            if(hME_US_C)  hME_US_C->Write();
            if(hCorr_US_C)  hCorr_US_C->Write();
            if(hSE_US_D)  hSE_US_D->Write();
            if(hME_US_D)  hME_US_D->Write();
            if(hCorr_US_D)  hCorr_US_D->Write();

            //B+C-D
            if(hCorr_SS_BPlusCMinusD) hCorr_SS_BPlusCMinusD->Write();
            if(hCorr_US_BPlusCMinusD) hCorr_US_BPlusCMinusD->Write();

            if (hCorr_SS_BPlusC) hCorr_SS_BPlusC->Write();  
            if (hCorr_US_BPlusC) hCorr_US_BPlusC->Write();



          }

          if (hCorr_LL)    hCorr_LL->Write();
          if (hCorr_ALAL)  hCorr_ALAL->Write();
          if (hCorr_US)    hCorr_US->Write();
         

          if (fLL.f)       fLL.f->Write();
          if (fALAL.f)     fALAL.f->Write();
          if (fUS.f)       fUS.f->Write();

          // additional same-sign average objects
          if (hSE_SS)      hSE_SS->Write();
          if (hME_SS)      hME_SS->Write();
          if (hCorr_SS)    hCorr_SS->Write();
          if (fSS.f)       fSS.f->Write();

          c->Write();
        }
      } // rbin

      // ---------------- Final P vs R canvas ----------------
      auto cP = std::make_unique<TCanvas>(
        Form("cP_vsR_%s_cent%dto%d", comb.Data(), centlow, centhigh),
        "P vs R", 1800, 1400
      );
      cP->SetLeftMargin(0.13);
      cP->SetBottomMargin(0.12);

      hP_LL_vsR->SetTitle(Form("Spin correlation P vs R (cent %d-%d, %s);R;P",
                               centlow, centhigh, comb.Data()));
      hP_LL_vsR->GetYaxis()->SetRangeUser(-0.2, 0.2);

      hP_LL_vsR->SetMarkerStyle(20);
      hP_ALAL_vsR->SetMarkerStyle(21);
      hP_SS_vsR->SetMarkerStyle(34);
      hP_US_vsR->SetMarkerStyle(22);

      hP_LL_vsR->Draw("E1");
      hP_ALAL_vsR->Draw("E1 SAME");
      hP_SS_vsR->Draw("E1 SAME");
      hP_US_vsR->Draw("E1 SAME");

      TLegend leg(0.50, 0.63, 0.89, 0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.AddEntry(hP_LL_vsR.get(),   "#Lambda#Lambda", "lp");
      leg.AddEntry(hP_ALAL_vsR.get(), "#bar{#Lambda}#bar{#Lambda}", "lp");
      leg.AddEntry(hP_SS_vsR.get(),   "Same-sign average", "lp");
      leg.AddEntry(hP_US_vsR.get(),   "Unlike-sign average", "lp");
      leg.Draw();

      const TString pngP = Form("%s/P_vsR_%s_cent%dto%d.png", pngDir.Data(), comb.Data(), centlow, centhigh);
      cP->SaveAs(pngP);

      auto cP_SB = std::make_unique<TCanvas>(
        Form("cP_SB_vsR_%s_cent%dto%d", comb.Data(), centlow, centhigh),
        "P vs R", 1800, 1400
      );
      cP_SB->SetLeftMargin(0.13);
      cP_SB->SetBottomMargin(0.12);

      hP_US_BB_vsR->SetTitle(Form("Spin correlation P vs R;R;P"
                               ));
      hP_US_BB_vsR->GetYaxis()->SetRangeUser(-0.2, 0.2);

      hP_US_BB_vsR->SetMarkerStyle(22);
      hP_US_BB_vsR->SetLineColor(kGreen+2);
      hP_US_BB_vsR->SetMarkerColor(kGreen+2);
      hP_US_BB_vsR->Draw("E1");


      if(do_YieldSS){
        hP_US_SB_vsR->SetMarkerStyle(20);
        hP_US_SB_vsR->SetLineColor(kRed);
        hP_US_SB_vsR->SetMarkerColor(kRed);
        
        hP_US_SB_vsR->Draw("E1 SAME");
        hP_SS_SB_vsR->SetMarkerStyle(34);
        hP_SS_SB_vsR->SetLineColor(kBlue);
        hP_SS_SB_vsR->SetMarkerColor(kBlue);
        hP_SS_SB_vsR->Draw("E1 SAME");
      }
     
      hP_SS_BB_vsR->SetMarkerStyle(33);
      hP_SS_BB_vsR->SetLineColor(kBlack);
      hP_SS_BB_vsR->SetMarkerColor(kBlack);
      hP_SS_BB_vsR->Draw("E1 SAME");

      TLegend leg1(0.50, 0.63, 0.89, 0.88);
      leg1.SetBorderSize(0);
      leg1.SetFillStyle(0);
      if(do_YieldSS){
        leg1.AddEntry(hP_US_SB_vsR.get(),   "2D fit;SB:Unlike-sign average", "lp");
        leg1.AddEntry(hP_SS_SB_vsR.get(),   "2D fit;SB:Same-sign average", "lp");
        leg1.AddEntry(hP_US_BB_vsR.get(),   "2D fit;BB:Unlike-sign average", "lp");
        leg1.AddEntry(hP_SS_BB_vsR.get(),   "2D fit;BB:Same-sign average", "lp");
      }else{
        leg1.AddEntry(hP_US_BB_vsR.get(),   "2D side band;BB:Unlike-sign average", "lp");
        leg1.AddEntry(hP_SS_BB_vsR.get(),   "2D side band;BB:Same-sign average", "lp");
      }
    
      leg1.Draw();
      //1.default options: ABCD method: signal:Double Gaussian ,background: Bern2,Bern3,Chebyshev2
      //2.fit method:combined fit method for signal and background:Double Gaussian ,background: Bern3,Bern4,Chebyshev2
      TString sg[] = {"Double Gaussian", "Bern3", "Chebyshev2"};
      TString bg[] = {"Bern3", "Bern4", "Chebyshev2"};

      TLatex label;
      label.SetNDC();
      label.SetTextSize(0.04);
      label.DrawLatex(0.15, 0.15, Form("Sig=%s, Bkg=%s", sg[sigstat].Data(), bg[bkgstat].Data()));

      const TString pngP_SB = Form("%s/P_vsR_SB_BB_%s_cent%dto%d.png", pngDir.Data(), comb.Data(), centlow, centhigh);
      cP_SB->SaveAs(pngP_SB);


      dComb->cd();
      hP_LL_vsR->Write("hP_LL_vsR");
      hP_ALAL_vsR->Write("hP_ALAL_vsR");
      hP_SS_vsR->Write("hP_SS_vsR");
      hP_US_vsR->Write("hP_US_vsR");
      hP_US_SB_vsR->Write("hP_US_SB_vsR");
      hP_US_BB_vsR->Write("hP_BCD_US_BB_vsR");
      hP_SS_BB_vsR->Write("hP_BCD_SS_BB_vsR");
      cP->Write("cP_vsR");
      cP_SB->Write("cP_vsR_SB");
      hP_SS_BPlusCMinusD_vsR->Write("hP_SS_BPlusCMinusD_vsR");
      hP_US_BPlusCMinusD_vsR->Write("hP_US_BPlusCMinusD_vsR");
      hP_US_BPlusC_vsR->Write("hP_US_BPlusC_vsR");
      hP_SS_BPlusC_vsR->Write("hP_SS_BPlusC_vsR");

      hP_LL_SPlusB_vsR->Write("hP_LL_SPlusB_vsR");
      hP_ALAL_SPlusB_vsR->Write("hP_ALAL_SPlusB_vsR");

      if(hP_SS_SPlusB_vsR)  hP_SS_SPlusB_vsR->Write();
      if(hP_US_SPlusB_vsR)  hP_US_SPlusB_vsR->Write();

      cP->Draw();
      cP_SB->Draw();
    }
  }

  fAll->Write();
  fAll->Close();

  std::cout << "\nDone.\n"
            << "  PNGs in: " << pngDir << "\n"
            << "  ROOT(all comb): " << outRootAll << "\n";
}

void UncertaintyOfSS(bool SS = true){

  TString SSandUS = SS ?"SS":"US";
  string histPathinput = SS ? "/hP_SS_vsR" : "/hP_US_vsR";


  vector<string> methods = {
    "YieldSS",
    "do_YieldSS_BCD"
  };

  vector<string> configs = {
    "sig0_bkg0",
    "sig0_bkg1",
    "sig0_bkg2"   // <-- default 
    //"sig2_bkg0",
    //"sig2_bkg1",
    //"sig2_bkg2"
  };
  //Open the result root file

  string prefix = "/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_";
  string suffix = "_Rdep_Output_pp_";
  string histPath = histPathinput;

  gSystem->mkdir(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/%s_Uncertainty/", SSandUS.Data()), kTRUE);

  string f_default_path = prefix + methods[1] + suffix + configs[2] + ".root";
  TFile* f_default = TFile::Open(f_default_path.c_str(), "READ");
  TH1D* h_default = (TH1D*)f_default->Get((configs[2] + histPath).c_str());
  //h_default->Draw();
  for(string method : methods){
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    

    for(string configs : configs){
      string defaultCfg = configs;
      string defaultFile = prefix + method + suffix + defaultCfg + ".root";

      TFile* f = TFile::Open(defaultFile.c_str(), "READ");
      
      TH1D* h = (TH1D*)f->Get((defaultCfg + histPath).c_str());
      h->Draw();
      c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/%s_Uncertainty/%s_%s.png", SSandUS.Data(), SSandUS.Data(), (method+defaultCfg).c_str()));

    }
     delete c;
  }
  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/%s_Uncertainty/%s_All_Differences.root", SSandUS.Data(), SSandUS.Data()), "RECREATE");
  // do uncertianty plot
  for(string method : methods){
    

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    
    for(string configs : configs){
      if(method == methods[0] && configs == configs[2]) continue;

      string defaultCfg = configs;
      string defaultFile = prefix + method + suffix + defaultCfg + ".root";

      TFile* f = TFile::Open(defaultFile.c_str(), "READ");
      TH1D* h = (TH1D*)f->Get((defaultCfg + histPath).c_str());
      TH1D* h_write = (TH1D*)h->Clone(Form("%s_%s_%s",SSandUS.Data(),method.c_str(), configs.c_str()));

      TH1D* h_diff = (TH1D*)h->Clone(Form("%s_diff_%s_%s",SSandUS.Data(),method.c_str(), configs.c_str()));
      h_diff->Add(h_default, -1);
      h_diff->GetYaxis()->SetTitle("Difference to default side band method+(sig0_bkg2)");
      h_diff->GetYaxis()->SetRangeUser(-0.01, 0.01);
      h_diff->SetLineWidth(2);
      outFile->cd();
      h_write->Write();
      h_diff->Write();
      h_diff->Draw("HIST");
      
      c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/%s_Uncertainty/%s_Diff_%s_%s.png",
                           SSandUS.Data(),SSandUS.Data(), method.c_str(), configs.c_str()));

    }
    delete c;
  }
  outFile->Write();
  outFile->Close();
}

void drawUncertaintyPlots(){

  //Draw the SS and US in sig0_bkg2
  TCanvas* c1 = new TCanvas("c1", "c1", 1800, 1200);
  TFile* f_US = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/US_Uncertainty/US_All_Differences.root", "READ");
  TFile* f_SS = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SS_Uncertainty/SS_All_Differences.root", "READ");

  TH1D* h_US_YieldSS_sig0_bkg2 = (TH1D*)f_US->Get("US_do_YieldSS_BCD_sig0_bkg2");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetRangeUser(-0.01, 0.01);

  h_US_YieldSS_sig0_bkg2->SetLineColor(kBlue);
  h_US_YieldSS_sig0_bkg2->SetLineWidth(2);
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitle("P_{#Lambda_{}#bar{#Lambda}}");
  TH1D* h_SS_YieldSS_sig0_bkg2 = (TH1D*)f_SS->Get("SS_do_YieldSS_BCD_sig0_bkg2");
  h_SS_YieldSS_sig0_bkg2->SetLineColor(kRed);
  h_SS_YieldSS_sig0_bkg2->SetLineWidth(2);
 
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetRangeUser(-0.2, 0.2);
  h_US_YieldSS_sig0_bkg2->Draw();
  h_SS_YieldSS_sig0_bkg2->Draw("SAME");  
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(h_US_YieldSS_sig0_bkg2, "Unlike-Sign:side band method;sig0+bkg2", "l");
  legend->AddEntry(h_SS_YieldSS_sig0_bkg2, "Like-Sign:side band method;sig0+bkg2", "l");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();
  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/SSandUS_fitmethod_sig0_bkg2.png");
  delete c1;
  // Draw the difference of US with two method
  TFile* outFile = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "RECREATE");

  TCanvas* c2 = new TCanvas("c2", "c2", 1800, 800);
  c2->Divide(2, 1);
  c2->cd(1);
  gPad->SetLeftMargin(0.15);
  TH1D* h_US_Diff_sig0_bkg0 = (TH1D*)f_US->Get("US_diff_YieldSS_sig0_bkg0");
  h_US_Diff_sig0_bkg0->GetYaxis()->SetRangeUser(-0.01, 0.01);
  h_US_Diff_sig0_bkg0->GetYaxis()->SetTitle("P-P_{default}");
  TH1D* h_US_Diff_sig0_bkg1 = (TH1D*)f_US->Get("US_diff_YieldSS_sig0_bkg1");
  TH1D* h_US_Diff_sig0_bkg2 = (TH1D*)f_US->Get("US_diff_YieldSS_sig0_bkg2");
  h_US_Diff_sig0_bkg0->SetLineColor(kBlue);
  h_US_Diff_sig0_bkg1->SetLineColor(kRed);
  h_US_Diff_sig0_bkg2->SetLineColor(kGreen+2);
  h_US_Diff_sig0_bkg0->Draw("HIST");
  h_US_Diff_sig0_bkg1->Draw("same HIST");
  h_US_Diff_sig0_bkg2->Draw("same HIST");
  

  TLegend* legend2 = new TLegend(0.2, 0.6, 0.9, 0.9);
  legend2->AddEntry(h_US_Diff_sig0_bkg0, "Unlike-Sign:fit method;Double Gaussian+Bern3", "l");
  legend2->AddEntry(h_US_Diff_sig0_bkg1, "Unlike-Sign:fit method;Double Gaussian+Bern4", "l");
  legend2->AddEntry(h_US_Diff_sig0_bkg2, "Unlike-Sign:fit method;Double Gaussian+Cheby2", "l");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->Draw();
  c2->cd(2);
  
  TH1D* h_US_Diff_BCD_sig0_bkg0 = (TH1D*)f_US->Get("US_diff_do_YieldSS_BCD_sig0_bkg0");
  TH1D* h_US_Diff_BCD_sig0_bkg1 = (TH1D*)f_US->Get("US_diff_do_YieldSS_BCD_sig0_bkg1");
  TH1D* h_US_Diff_BCD_sig0_bkg2 = (TH1D*)f_US->Get("US_diff_do_YieldSS_BCD_sig0_bkg2");
  h_US_Diff_BCD_sig0_bkg0->GetYaxis()->SetRangeUser(-0.01, 0.01);

  h_US_Diff_BCD_sig0_bkg0->SetLineColor(kBlue);
  h_US_Diff_BCD_sig0_bkg1->SetLineColor(kRed);
  h_US_Diff_BCD_sig0_bkg2->SetLineColor(kGreen+2);

  h_US_Diff_BCD_sig0_bkg0->SetLineStyle(1);
  h_US_Diff_BCD_sig0_bkg1->SetLineStyle(2);
  h_US_Diff_BCD_sig0_bkg2->SetLineStyle(2);

  h_US_Diff_BCD_sig0_bkg0->Draw("HIST");
  h_US_Diff_BCD_sig0_bkg1->Draw("same HIST");
  h_US_Diff_BCD_sig0_bkg2->Draw("same HIST");
  TLegend* legend3 = new TLegend(0.1, 0.6, 0.9, 0.9);
  legend3->AddEntry(h_US_Diff_BCD_sig0_bkg0, "Unlike-Sign:side band method;Double Gaussian+Bern3", "l");
  legend3->AddEntry(h_US_Diff_BCD_sig0_bkg1, "Unlike-Sign:side band method;Double Gaussian+Bern4", "l");
  legend3->AddEntry(h_US_Diff_BCD_sig0_bkg2, "Default:Unlike-Sign:side band method;Double Gaussian+Cheby2", "l");

  legend3->SetBorderSize(0);
  legend3->SetFillStyle(0);
  legend3->Draw();
  

  outFile->cd();
  outFile->Write();
  c2->Write("c_US_Differences");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/US_Diff.png");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/US_Diff.pdf");
  delete c2;
  outFile->Write();
  outFile->Close();

  // Draw the difference of US with two method

  TCanvas* c_like_sign = new TCanvas("c_like_sign", "c_like_sign", 1800, 800);
  c_like_sign->Divide(2, 1);
  c_like_sign->cd(1);
  gPad->SetLeftMargin(0.15);
  TH1D* h_SS_Diff_sig0_bkg0 = (TH1D*)f_SS->Get("SS_diff_YieldSS_sig0_bkg0");
  h_SS_Diff_sig0_bkg0->GetYaxis()->SetRangeUser(-0.01, 0.01);
  h_SS_Diff_sig0_bkg0->GetYaxis()->SetTitle("P-P_{default}");
  TH1D* h_SS_Diff_sig0_bkg1 = (TH1D*)f_SS->Get("SS_diff_YieldSS_sig0_bkg1");
  TH1D* h_SS_Diff_sig0_bkg2 = (TH1D*)f_SS->Get("SS_diff_YieldSS_sig0_bkg2");
  h_SS_Diff_sig0_bkg0->SetLineColor(kBlue);
  h_SS_Diff_sig0_bkg1->SetLineColor(kRed);
  h_SS_Diff_sig0_bkg2->SetLineColor(kGreen+2);
  h_SS_Diff_sig0_bkg0->Draw("HIST");
  h_SS_Diff_sig0_bkg1->Draw("same HIST");
  h_SS_Diff_sig0_bkg2->Draw("same HIST");
  

  TLegend* legend_like_sign = new TLegend(0.2, 0.6, 0.9, 0.9);
  legend_like_sign->AddEntry(h_SS_Diff_sig0_bkg0, "like-Sign:fit method;Double Gaussian+Bern3", "l");
  legend_like_sign->AddEntry(h_SS_Diff_sig0_bkg1, "like-Sign:fit method;Double Gaussian+Bern4", "l");
  legend_like_sign->AddEntry(h_SS_Diff_sig0_bkg2, "like-Sign:fit method;Double Gaussian+Cheby2", "l");
  legend_like_sign->SetBorderSize(0);
  legend_like_sign->SetFillStyle(0);
  legend_like_sign->Draw();
  
  

  c2->cd(2);
  
  TH1D* h_SS_Diff_BCD_sig0_bkg0 = (TH1D*)f_SS->Get("SS_diff_do_YieldSS_BCD_sig0_bkg0");
  TH1D* h_SS_Diff_BCD_sig0_bkg1 = (TH1D*)f_SS->Get("SS_diff_do_YieldSS_BCD_sig0_bkg1");
  TH1D* h_SS_Diff_BCD_sig0_bkg2 = (TH1D*)f_SS->Get("SS_diff_do_YieldSS_BCD_sig0_bkg2");
  h_SS_Diff_BCD_sig0_bkg0->GetYaxis()->SetRangeUser(-0.01, 0.01);

  h_SS_Diff_BCD_sig0_bkg0->SetLineColor(kBlue);
  h_SS_Diff_BCD_sig0_bkg1->SetLineColor(kRed);
  h_SS_Diff_BCD_sig0_bkg2->SetLineColor(kGreen+2);

  h_SS_Diff_BCD_sig0_bkg0->SetLineStyle(1);
  h_SS_Diff_BCD_sig0_bkg1->SetLineStyle(2);
  h_SS_Diff_BCD_sig0_bkg2->SetLineStyle(2);

  h_SS_Diff_BCD_sig0_bkg0->Draw("HIST");
  h_SS_Diff_BCD_sig0_bkg1->Draw("same HIST");
  h_SS_Diff_BCD_sig0_bkg2->Draw("same HIST");
  TLegend* legend_like_sign1 = new TLegend(0.1, 0.6, 0.9, 0.9);
  legend_like_sign1->AddEntry(h_SS_Diff_BCD_sig0_bkg0, "like-Sign:side band method;Double Gaussian+Bern3", "l");
  legend_like_sign1->AddEntry(h_SS_Diff_BCD_sig0_bkg1, "like-Sign:side band method;Double Gaussian+Bern4", "l");
  legend_like_sign1->AddEntry(h_SS_Diff_BCD_sig0_bkg2, "Default:like-Sign:side band method;Double Gaussian+Cheby2", "l");

  legend_like_sign1->SetBorderSize(0);
  legend_like_sign1->SetFillStyle(0);
  legend_like_sign1->Draw();
  c_like_sign->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/SS_Diff.png");
  c_like_sign->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/SS_Diff.pdf");
  
  //delete c_like_sign;
  //calculate the uncertainty from signal extraction

  TFile* outFile1 = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "UPDATE");

  TCanvas* c3 = new TCanvas("c3", "c3", 1200, 800);
  TH1D* h_RMS_US = (TH1D*)h_US_Diff_BCD_sig0_bkg0->Clone("h_RMS_US");h_RMS_US->Reset();
  h_RMS_US->GetYaxis()->SetTitle("RMS");
  gPad->SetLeftMargin(0.15);
  for(int bin = 1; bin <= h_US_YieldSS_sig0_bkg2->GetNbinsX(); ++bin){
    double diff1 = h_US_Diff_sig0_bkg0->GetBinContent(bin);
    double diff2 = h_US_Diff_sig0_bkg1->GetBinContent(bin);
    double diff3 = h_US_Diff_sig0_bkg2->GetBinContent(bin);
    double diff4 = h_US_Diff_BCD_sig0_bkg0->GetBinContent(bin);
    double diff5 = h_US_Diff_BCD_sig0_bkg1->GetBinContent(bin);
    double values[] = {diff1, diff2, diff3, diff4, diff5};
    double rms = TMath::RMS(5, values);
    cout << "Bin " << bin << ": RMS = " << rms << endl;
    h_RMS_US->SetBinContent(bin, rms);

  }
  //h_RMS_US->Draw();
  TH1D* h_SS_YieldSS_sig0_bkg0 = (TH1D*)f_SS->Get("SS_YieldSS_sig0_bkg0");
  TH1D* h_SS_YieldSS_sig0_bkg1 = (TH1D*)f_SS->Get("SS_YieldSS_sig0_bkg1");
  //TH1D* h_SS_YieldSS_sig0_bkg2 = (TH1D*)f_SS->Get("SS_YieldSS_sig0_bkg2");
  //TH1D* h_SS_Diff_BCD_sig0_bkg0 = (TH1D*)f_SS->Get("SS_diff_do_YieldSS_BCD_sig0_bkg0");
  //TH1D* h_SS_Diff_BCD_sig0_bkg1 = (TH1D*)f_SS->Get("SS_diff_do_YieldSS_BCD_sig0_bkg1");

  TH1D* h_RMS_SS = (TH1D*)h_SS_YieldSS_sig0_bkg2->Clone("h_RMS_SS");h_RMS_SS->Reset();
  h_RMS_SS->GetYaxis()->SetTitle("RMS");
  h_RMS_SS->GetYaxis()->SetRangeUser(-0.01, 0.01);
  gPad->SetLeftMargin(0.15);
  for(int bin = 1; bin <= h_US_YieldSS_sig0_bkg2->GetNbinsX(); ++bin){
    double diff1 = h_SS_YieldSS_sig0_bkg0->GetBinContent(bin);
    double diff2 = h_SS_YieldSS_sig0_bkg1->GetBinContent(bin);
    double diff3 = h_SS_YieldSS_sig0_bkg2->GetBinContent(bin);
    double diff4 = h_SS_Diff_BCD_sig0_bkg0->GetBinContent(bin);
    double diff5 = h_SS_Diff_BCD_sig0_bkg1->GetBinContent(bin);
    double values[] = {diff1, diff2, diff3, diff4, diff5};
    double rms = TMath::RMS(5, values);
    cout << "Bin " << bin << ": RMS = " << rms << endl;
    h_RMS_SS->SetBinContent(bin, rms);
  }


  outFile1->cd();
  h_RMS_US->Write("h_RMS_US");
  h_RMS_SS->Write("h_RMS_SS");
  outFile1->Write();
  outFile1->Close();
  delete c3;

  

}



//add RMS to the uncertainty of the default method
void drawPloarizationcorrPlots(){
  
  TFile* outFile = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "UPDATE");


  TCanvas* c1 = new TCanvas("c1", "c1", 1800, 1200);
  gPad->SetLeftMargin(0.15);
  TFile* f_US = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/US_Uncertainty/US_All_Differences.root", "READ");
  TFile* f_SS = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SS_Uncertainty/SS_All_Differences.root", "READ");

  TH1D* h_US_YieldSS_sig0_bkg2 = (TH1D*)f_US->Get("US_do_YieldSS_BCD_sig0_bkg2");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetRangeUser(-0.02, 0.05);

  h_US_YieldSS_sig0_bkg2->SetLineColor(kBlue);
  h_US_YieldSS_sig0_bkg2->SetLineWidth(2);
  h_US_YieldSS_sig0_bkg2->SetTitle("Spin correlation P_{#Lambda_{1},#Lambda_{2}} as a function of pair separation #DeltaR");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitle("P_{#Lambda_{1},#Lambda_{2}}");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitleSize(0.05);
  //h_US_YieldSS_sig0_bkg2->GetXaxis()->SetTitleSize(0.06);
  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}");
  TH1D* h_SS_YieldSS_sig0_bkg2 = (TH1D*)f_SS->Get("SS_do_YieldSS_BCD_sig0_bkg2");
  h_SS_YieldSS_sig0_bkg2->SetLineColor(kRed);
  h_SS_YieldSS_sig0_bkg2->SetLineWidth(2);
  
  
 

  h_US_YieldSS_sig0_bkg2->Draw("E1");
  h_US_YieldSS_sig0_bkg2->SetMarkerSize(3);
  h_US_YieldSS_sig0_bkg2->SetMarkerColor(kBlue);
  h_US_YieldSS_sig0_bkg2->GetXaxis()->CenterTitle();
  h_US_YieldSS_sig0_bkg2->GetYaxis()->CenterTitle();
  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetTitleOffset(1.30);
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitleOffset(1.0);

  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetLabelSize(0.05);  // X刻度变大
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetLabelSize(0.05);  // Y刻

  h_SS_YieldSS_sig0_bkg2->SetMarkerSize(3);
  h_SS_YieldSS_sig0_bkg2->SetMarkerColor(kRed);
  h_SS_YieldSS_sig0_bkg2->Draw("E1 SAME");  
 
  TFile* RMShistFile = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "READ");
  TH1D* h_RMS = (TH1D*)RMShistFile->Get("h_RMS_US");
  TH1D* h_RMS_SS = (TH1D*)RMShistFile->Get("h_RMS_SS");

  
    double box_width = 0.08;  
    double box_alpha = 0.3;   
    
    for(int bin = 1; bin <= h_US_YieldSS_sig0_bkg2->GetNbinsX(); ++bin){
        double bin_center = h_US_YieldSS_sig0_bkg2->GetBinCenter(bin);
        double bin_left = bin_center - box_width/2;
        double bin_right = bin_center + box_width/2;
        double pol_value = h_US_YieldSS_sig0_bkg2->GetBinContent(bin);
        double pol_value_SS = h_SS_YieldSS_sig0_bkg2->GetBinContent(bin);
        double sys_error = h_RMS->GetBinContent(bin);
        double sys_error_SS = h_RMS_SS->GetBinContent(bin);
        
        if (sys_error > 0) {
           
            double box_bottom = pol_value - sys_error;
            double box_top = pol_value + sys_error;
            TBox *box = new TBox(bin_left, box_bottom, bin_right, box_top);
            box->SetFillColor(kGray);
            box->SetFillStyle(0);  
            box->SetLineColor(kBlue);
            box->SetLineWidth(2);
            box->SetFillColorAlpha(kGray, box_alpha);  
            //box->Draw("SAME");

            double box_bottom_SS = pol_value_SS - sys_error_SS;
            double box_top_SS = pol_value_SS + sys_error_SS;
            TBox *box_SS = new TBox(bin_left, box_bottom_SS, bin_right, box_top_SS);
            box_SS->SetFillColor(kGray);
            box_SS->SetFillStyle(0);  
            box_SS->SetLineColor(kRed);
            box_SS->SetLineWidth(2);
            box_SS->SetFillColorAlpha(kGray, box_alpha);  
            //box_SS->Draw("SAME");
            

        }
    }
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.18, 0.85, "ALICE");
    tex->DrawLatex(0.18, 0.80, "p+p #sqrt{s} = 13.6 TeV");
    tex->DrawLatex(0.18, 0.75, "|y| < 0.5");
    tex->DrawLatex(0.18, 0.70, "0.8 GeV/c < p_{T} < 3 GeV/c");
    tex->DrawLatex(0.18, 0.65, "Helicity frame");
    tex->DrawLatex(0.18, 0.60, "Signal extraction:2D side band method");
    gPad->SetGridx();  
    gPad->SetGridy(); 
   

    TLine *line = new TLine(
    h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmin(),  // 左边界（自动匹配）
    0,                                              // y=0
    h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmax(),  // 右边界（自动匹配）
    0                                               // y=0
);

 
  line->SetLineStyle(10);
  
  line->SetLineColor(kBlack);
  // 设置宽度
  line->SetLineWidth(5);

  // 绘制线条
  line->Draw("same");

  double val1  = 0.015;
  double err1  = 0.002;
  double ymin1 = val1 - err1;
  double ymax1 = val1 + err1;
  double xmin1 = h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmin();
  double xmax1 = h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmax();

  TBox *box1 = new TBox(xmin1, ymin1, xmax1, ymax1);
  box1->SetFillColor(kRed);   // 灰色填充
  box1->SetFillStyle(3001);    // 透明灰色（不挡住曲线）
  box1->Draw();

  // 画中心线
  TLine *line1 = new TLine(xmin1, val1, xmax1, val1);
  line1->SetLineColor(kBlack);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);
  line1->Draw("same");

  TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.9);
  legend->AddEntry(h_US_YieldSS_sig0_bkg2, "Unlike-Sign : #Lambda#bar{#Lambda}", "lp");
  legend->AddEntry(h_SS_YieldSS_sig0_bkg2, "Like-Sign : #Lambda#Lambda", "lp");
  legend->AddEntry(box1, "STAR simulation:BJ:0.015 #pm 0.002", "f");  
  
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();


 



  TCanvas* c2= new TCanvas("c2", "c2", 1800, 1200);
  TFile* f1 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root", "READ");

  TH1D* h_LL = (TH1D*)f1->Get("sig0_bkg2/hP_LL_vsR;1");
  TH1D* h_ALAL = (TH1D*)f1->Get("sig0_bkg2/hP_ALAL_vsR;1");
  h_LL->GetYaxis()->SetRangeUser(-0.02, 0.05);
  h_LL->SetLineColor(kRed);
  h_LL->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}");

  h_LL->SetTitle("Spin correlation P_{#Lambda_{1},#Lambda_{2}} as a function of pair separation #DeltaR");
  h_LL->GetYaxis()->SetTitle("P_{#Lambda_{1},#Lambda_{2}}");
 
  
  


  h_LL->Draw("E1");
  h_LL->SetMarkerSize(3);
  h_LL->SetMarkerColor(kRed);
  h_LL->GetXaxis()->CenterTitle();
  h_LL->GetYaxis()->CenterTitle();
  h_LL->GetXaxis()->SetTitleOffset(1.30);
  h_LL->GetYaxis()->SetTitleOffset(1.0);
  h_LL->SetLineWidth(2);
  h_LL->GetYaxis()->SetTitleSize(0.05);
  h_LL->Draw();



  
  h_ALAL->SetLineColor(kBlue);
  h_ALAL->Draw("E1 SAME");
  h_ALAL->SetMarkerSize(3);
  h_ALAL->SetMarkerColor(kBlue);
  h_ALAL->SetLineWidth(2);
  h_ALAL->Draw("SAME");

  TLegend* legend1= new TLegend(0.6, 0.7, 0.88, 0.9);
  legend1->AddEntry(h_LL, "#Lambda#Lambda", "lp");
  legend1->AddEntry(h_ALAL, "#bar{#Lambda}#bar{#Lambda}", "lp");
  
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->Draw();
  gPad->SetGridx();
  gPad->SetGridy();

  TLatex *tex1 = new TLatex();
  tex1->SetNDC();
  tex1->SetTextSize(0.03);
  tex1->DrawLatex(0.18, 0.85, "ALICE");
  tex1->DrawLatex(0.18, 0.80, "p+p #sqrt{s} = 13.6 TeV");
  tex1->DrawLatex(0.18, 0.75, "|y| < 0.5");
  tex1->DrawLatex(0.18, 0.70, "0.8 GeV/c < p_{T} < 3 GeV/c");
  tex1->DrawLatex(0.18, 0.65, "Helicity frame");
  tex1->DrawLatex(0.18, 0.60, "Signal extraction:2D side band method");

  

   line->Draw("same");




 
  TFile* f2 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root", "READ");

  TCanvas* c3= new TCanvas("c3", "c3", 1800, 1200);
  TCanvas* h_BB = (TCanvas*)f2->Get("sig0_bkg2/cP_vsR_SB;1");
  h_BB->Draw();

  outFile->cd();
  c1->Write("c_P_with_sys");
  c2->Write("h_LL");
 // c3->Write("h_BB");
  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/P_.png");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_LL.png");
  c3->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_BB.png");

  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/P_.pdf");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_LL.pdf");
  c3->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_BB.pdf");
  outFile->Save();
  delete c1;
  delete c2;
 // delete c3;
}


void SEOverME(bool SS = true){
  //draw SE distribution for US
  TString SSandUS = SS ?"SS":"US";
  TString SSorUS  = SS ?"Like-sign(#Lambda#Lambda+#bar{#Lambda}#bar{#Lambda})":"Unlike-sign(#Lambda#bar{#Lambda}+#bar{#Lambda}#Lambda)";
  

  vector<string> methods = {
    "YieldSS",//0
    "do_YieldSS_BCD"//1
  };

  vector<string> configs = {
    "sig0_bkg0",//0
    "sig0_bkg1",//1
    "sig0_bkg2"   // <-- default 2
    //"sig2_bkg0",
    //"sig2_bkg1",
    //"sig2_bkg2"
  };
  
  gSystem->mkdir(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_distribution/"), kTRUE);

  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_distribution/%s_SEOverME.root",SSandUS.Data()),"recreate");

  for(int Rbin = 0; Rbin <= 5; ++Rbin){
    TCanvas* c = new TCanvas("c", Form("c_Rbin%d", Rbin), 1600, 1200);
    c->Divide(1, 2,0); 
    c->cd(1);
    gPad->SetBottomMargin(0.0); 
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.);   
    gPad->SetRightMargin(0.05);
    
     
    string prefix = "/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_";
    string suffix = "_Rdep_Output_pp_";
    string f_path = prefix + methods[1] + suffix + configs[2] + ".root";
    TFile* f_default = TFile::Open(f_path.c_str(), "READ");
    TH1D* hSE_US_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hSE_%s_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    TH1D* hME_US_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hME_%s_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    //hSE_US_sig0_bkg2->GetYaxis()->SetRangeUser(0,hSE_US_sig0_bkg2->GetMaximum()*1.5);

    hSE_US_sig0_bkg2->Scale(1.0/hSE_US_sig0_bkg2->Integral());
    hME_US_sig0_bkg2->Scale(1.0/hME_US_sig0_bkg2->Integral());

    TH1D* hRatio_SE_ME = (TH1D*)hSE_US_sig0_bkg2->Clone("hRatio_SE_ME");
    hRatio_SE_ME->Divide(hME_US_sig0_bkg2);

   

    hSE_US_sig0_bkg2->GetYaxis()->SetRangeUser(0.,0.25);
    hSE_US_sig0_bkg2->SetTitle("");
    hSE_US_sig0_bkg2->GetYaxis()->SetTitle("1/N dN/dcos#theta*");
    hSE_US_sig0_bkg2->SetLineColor(kRed);
    hSE_US_sig0_bkg2->SetMarkerColor(kRed);
    hSE_US_sig0_bkg2->SetLineWidth(2);
    hME_US_sig0_bkg2->SetLineWidth(2);
    
    
    

    //hSE_US_sig0_bkg2->SetMarkerStyle(20);
    //hME_US_sig0_bkg2->SetMarkerStyle(21);
    hSE_US_sig0_bkg2->Draw();
    hME_US_sig0_bkg2->Draw("SAME");

    const double REdges[] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};
    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextSize(0.05);
    title->SetTextAlign(22);
    title->DrawLatex(0.3, 0.95, Form("ALICE"));
    title->DrawLatex(0.3, 0.85, Form("0.8 GeV/c < p_{T} < 3 GeV/c"));
    title->DrawLatex(0.3, 0.75, Form("|y| < 0.5"));
    title->DrawLatex(0.3, 0.65, Form("Rbin %d: %.1f < R < %.1f", Rbin, REdges[Rbin], REdges[Rbin+1]));
    title->DrawLatex(0.3, 0.55, Form("Helicity frame"));

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hSE_US_sig0_bkg2, Form("SE_%s", SSorUS.Data()), "lp");
    legend->AddEntry(hME_US_sig0_bkg2, Form("ME_%s", SSorUS.Data()), "lp");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();

    
   
    /*
    c->cd(2);
   
    gPad->SetTopMargin(0.0);    
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.0);

    TH1D* hcosSS_SEerror = (TH1D*)hSE_US_sig0_bkg2->Clone("hRatio");hcosSS_SEerror->Reset();
    TH1D* hcosSS_MEerror = (TH1D*)hME_US_sig0_bkg2->Clone("hRatio");hcosSS_MEerror->Reset();
    hcosSS_SEerror->GetYaxis()->SetRangeUser(0, 0.004);
    hcosSS_SEerror->GetYaxis()->SetTitle("#sigma/N");
    hcosSS_SEerror->GetYaxis()->SetRangeUser(0, hcosSS_SEerror->GetMaximum()*1.5);
    for (int b=1; b<=hSE_US_sig0_bkg2->GetNbinsX(); ++b) {
      double se = hSE_US_sig0_bkg2->GetBinContent(b);
      double me = hME_US_sig0_bkg2->GetBinContent(b);
      double seErr = hSE_US_sig0_bkg2->GetBinError(b);
      double meErr = hME_US_sig0_bkg2->GetBinError(b);
      double se_error = seErr/se;
      double me_error = meErr/me;
      //cout<<b<<endl;
      hcosSS_SEerror->SetBinContent(b, se_error);
      hcosSS_MEerror->SetBinContent(b, me_error);
    }
    hcosSS_SEerror->Draw();
    hcosSS_MEerror->Draw("SAME");

    TLegend* legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(hcosSS_SEerror, Form("SE_%s", SSorUS.Data()), "lp");
    legend2->AddEntry(hcosSS_MEerror, Form("ME_%s", SSorUS.Data()), "lp");
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->Draw();
    */

    c->cd(2);
  
    gPad->SetTopMargin(0.0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(1);
  /*
    TH1D* hCorr_SS_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    hCorr_SS_sig0_bkg2->Draw();

    TF1* fit_SS_US_sig0_bkg2 = (TF1*)f_default->Get(Form("%s/rbin%d/fit_%s_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    //tex->DrawLatex(0.15, 0.85, Form("Fit: %f", fit_SS_US_sig0_bkg2->GetParameter(0)));
    */
    
    hRatio_SE_ME->GetYaxis()->SetRangeUser(0.8, 1.2);
    hRatio_SE_ME->GetYaxis()->SetTitle("SE/ME");
    hRatio_SE_ME->SetTitle("");
    hRatio_SE_ME->Draw();
    
    outFile->cd();
    outFile->Write();
    c->Write(Form("%s_SE_ME_and_Corr_Rbin%d", SSandUS.Data(), Rbin));
    
    c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_distribution/%s_SE_ME_and_Corr_Rbin%d.png", SSandUS.Data(), Rbin));
    c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_distribution/%s_SE_ME_and_Corr_Rbin%d.pdf", SSandUS.Data(), Rbin));

    delete c;

  }
  outFile->Close();
  

}

void x_P_Lambda_y(){
  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs.root"),"recreate");

  const int nMassBins1 = 6;
  
  double massBinCenter[nMassBins1+1] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};

    // ========== 1. 打开文件，读取直方图 ==========
    TFile* f_fit_fixed = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root","READ");
    if (!f_fit_fixed || f_fit_fixed->IsZombie()) {
        cout << "ERROR: 文件打开失败！" << endl;
        return;
    }
    TFile* f_fit_fixed3 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Sourav_analysis/Result_signalType0_bkgType2_fixedmusigma_1/output_corr_vsMass_cent00.root","READ");


    
    
    TH1D* hfit_unlike = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_US_vsR;1");   // ΛΛ̅ 异号
    TH1D* hfit_like_Lambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_LL_vsR;1");   // ΛΛ 同号
    TH1D* hfit_like_AntiLambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_ALAL_vsR;1");   // Λ̅ Λ̅  同号

    TH1D* hfit_unlike3 = (TH1D*)f_fit_fixed3->Get("hfit_Unlike_vsR;1");   // ΛΛ̅ 异号
    TH1D* hfit_like_Lambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle0_vsR;1");   // ΛΛ 同号
    TH1D* hfit_like_AntiLambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle3_vsR;1");   // Λ̅ Λ̅  同号
   
    for(int i = 0;i<=nMassBins1-1;i++){
    TCanvas *c = new TCanvas("c", Form("Spin-spin correlation: P vs #Lambda pairs_%d",i), 2000, 1200);
    c->Divide(2,1);
    
    c->cd(1);
    int nMassBins = hfit_unlike->GetNbinsX();

   
    

    TPad *pad1 = new TPad("pad1", "", 0.01, 0.01, 0.99, 0.99);
    pad1->SetLeftMargin(0.18);
    pad1->SetBottomMargin(0.15);
    pad1->SetTopMargin(0.1);
    pad1->SetRightMargin(0.05);
    pad1->SetGrid(0, 0);
    pad1->Draw();
    pad1->cd();

    
    TH2F *hframe1 = new TH2F("hframe1", "", 100, -0.25, 0.39, 4, 0.5, 4.5);
    hframe1->SetStats(0);
    hframe1->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
    hframe1->GetXaxis()->SetTitleOffset(1.4);
    hframe1->GetXaxis()->SetRangeUser(-0.245, 0.385);

   
    hframe1->GetYaxis()->SetBinLabel(1, "#Lambda#bar{#Lambda}");
    hframe1->GetYaxis()->SetBinLabel(2, "#Lambda#Lambda");
    hframe1->GetYaxis()->SetBinLabel(3, "#bar{#Lambda}#bar{#Lambda}");
    hframe1->GetYaxis()->SetBinLabel(4, "K^{0}_{S}K^{0}_{S}");
    hframe1->GetYaxis()->SetLabelSize(0.08);
    hframe1->GetYaxis()->SetLabelOffset(0.02);
    hframe1->GetXaxis()->CenterTitle();
    hframe1->GetYaxis()->SetTitle("");
    hframe1->GetXaxis()->SetTitleSize(0.05);
    
    hframe1->Draw();

    

   
    TLine *line0 = new TLine(0, 0.5, 0, 4.5);
    line0->SetLineStyle(2);
    line0->SetLineWidth(2);
    line0->SetLineColor(kBlack);
    line0->Draw("same");

    
   
      // 异号 ΛΛ̅
      double p_unlike    = hfit_unlike->GetBinContent(i+1);
      double err_unlike  = hfit_unlike->GetBinError(i+1);

      // 同号 ΛΛ
      double p_like     = hfit_like_Lambda->GetBinContent(i+1);
      double err_like   = hfit_like_Lambda->GetBinError(i+1);

      // 同号 Λ̅Λ̅
      double p_like_AntiLambda = hfit_like_AntiLambda->GetBinContent(i+1);
      double err_like_AntiLambda = hfit_like_AntiLambda->GetBinError(i+1);

      /////////////
      double p_unlike3    = hfit_unlike3->GetBinContent(i+1);
      double err_unlike3  = hfit_unlike3->GetBinError(i+1);

      // 同号 ΛΛ
      double p_like3     = hfit_like_Lambda3->GetBinContent(i+1);
      double err_like3  = hfit_like_Lambda3->GetBinError(i+1);

      // 同号 Λ̅Λ̅
      double p_like_AntiLambda3 = hfit_like_AntiLambda3->GetBinContent(i+1);
      double err_like_AntiLambda3 = hfit_like_AntiLambda3->GetBinError(i+1);

    
      const int n = 3;
      double x[n]  = {p_unlike,  p_like, p_like_AntiLambda };   // X = P 值
      double y[n]  = {1.0,       2.0, 3.0 };        // Y = 类别位置
      double ex[n] = {err_unlike, err_like, err_like_AntiLambda};// 误差棒
      double ey[n] = {0.0, 0.0,0.0,};

      double x3[n]  = {p_unlike3,  p_like3, p_like_AntiLambda3 };   // X = P 值
      double y3[n]  = {1.0,       2.0, 3.0 };        // Y = 类别位置
      double ex3[n] = {err_unlike3, err_like3, err_like_AntiLambda3};// 误差棒
      double ey3[n] = {0.0, 0.0,0.0};

      const int n_star1 = 4;
      double x_star1[n_star1]  = {0.181,  0.044, -0.067 ,0.035};   // X = P 值
      double y_star1[n_star1]  = {1.0,       2.0, 3.0 ,4.0};        // Y = 类别位置
      double ex_star1[n_star1] = {0.035, 0.109, 0.124,0.024};// 误差棒
      double ey_star1[n_star1] = {0.0, 0.0,0.0,0};


      const int n_star2 = 4;
      double x_star2[n_star1]  = {0.02,  -0.002, -0.022 ,-0.025};   // X = P 值
      double y_star2[n_star1]  = {1.0,       2.0, 3.0 ,4.0};        // Y = 类别位置
      double ex_star2[n_star1] = {0.023, 0.043, 0.052,0.014};// 误差棒
      double ey_star2[n_star1] = {0.0, 0.0,0.0,0};

      double sys_star1[n_star1] = {0.022, 0.022, 0.022, 0.018};
      double sys_star2[n_star1] = {0.022, 0.022, 0.022, 0.020};


      // ========== 8. 画带误差的点 ==========
     

      TGraphErrors *gr_star1 = new TGraphErrors(n_star1, x_star1, y_star1, ex_star1, ey_star1);
      gr_star1->SetMarkerStyle(20);
      gr_star1->SetMarkerSize(2.0);
      gr_star1->SetMarkerColor(kBlue);
      gr_star1->SetLineColor(kBlue);
      gr_star1->Draw("P same");


      TGraphErrors *gr_star2 = new TGraphErrors(n_star1, x_star2, y_star2, ex_star2, ey_star2);
      gr_star2->SetMarkerStyle(20);
      gr_star2->SetMarkerSize(2.0);
      gr_star2->SetMarkerColor(kBlue);
      gr_star2->SetLineColor(kBlue);
      //gr_star2->Draw("P same");

      

      TGraphErrors *gr_star1_sys = new TGraphErrors(n_star1, x_star1, y_star1, sys_star1, ey_star1);
      gr_star1_sys->SetLineColor(kBlue);       // 同色
      gr_star1_sys->SetLineStyle(2);           // 虚线（区分统计误差）
      gr_star1_sys->SetLineWidth(8);           // 细一点
      gr_star1_sys->SetMarkerSize(0);          // 不画点
      gr_star1_sys->Draw("EZ same");   

      TGraphErrors *gr_star2_sys = new TGraphErrors(n_star1, x_star2, y_star2, sys_star2, ey_star2);
      gr_star2_sys->SetLineColor(kBlue);       // 同色
      gr_star2_sys->SetLineStyle(2);           // 虚线（区分统计误差）
      gr_star2_sys->SetLineWidth(8);           // 细一点
      gr_star2_sys->SetMarkerSize(0);          // 不画点
      //gr_star2_sys->Draw("EZ same"); 

  


      TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(2.0);
      gr->SetMarkerColor(kRed);
      gr->SetLineColor(kRed);
      gr->Draw("P same");

    
      c->cd(2);
      // ========== 9. 文字标注 ==========
      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextSize(0.033);
      tex->SetTextColor(kRed);  
      tex->DrawLatex(0.10, 0.85, "ALICE");
      tex->DrawLatex(0.10, 0.80, "p+p #sqrt{s} = 13.6 TeV");
      tex->DrawLatex(0.10, 0.75, Form("#DeltaR = #sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}} : [%.2f-%.2f]", massBinCenter[i],massBinCenter[i+1]));
      tex->DrawLatex(0.10, 0.70, "|y| < 1.5");
      tex->DrawLatex(0.10, 0.65, "Helicity frame");

      
      TLatex *tex1 = new TLatex();
      tex1->SetNDC();
      tex1->SetTextSize(0.033);
      tex1->SetTextColor(kBlue);  
      tex1->DrawLatex(0.60, 0.85, "STAR");
      tex1->DrawLatex(0.60, 0.80, "p+p #sqrt{s} = 200 GeV");
      tex1->DrawLatex(0.60, 0.75, "|#Deltay| < 0.5,|#Delta#phi| < #pi/3");
      //tex1->DrawLatex(0.60, 0.75, "0.5 < |#Deltay| < 2.0 or |#Delta#phi| > #pi/3");
      tex1->DrawLatex(0.60, 0.70, "|y| < 1");
      tex1->DrawLatex(0.60, 0.65, "<p_{T,#Lambda}> = 1.35 GeV/c");
      tex1->DrawLatex(0.60, 0.60, "DCA_{#Lambda} < 1 cm");
      tex1->DrawLatex(0.60, 0.55, "STAR frame");




      TLegend *leg = new TLegend(0.10, 0.55, 0.82, 0.65);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.025);
      leg->AddEntry(gr, "Data", "lp");
      
      //leg->AddEntry(gr3, "2D fit method,Double Gaussian+Cheby2:fixed #mu,#sigma", "lp");
      leg->Draw("same");


      TLegend *leg1 = new TLegend(0.10, 0.45, 0.82, 0.55);
      leg1->SetBorderSize(0);
      leg1->SetFillStyle(0);
      leg1->SetTextSize(0.025);
      //leg1->AddEntry(gr, "Data", "lp");
      leg1->AddEntry(gr_star1, "Data", "lp");
      //leg->AddEntry(gr3, "2D fit method,Double Gaussian+Cheby2:fixed #mu,#sigma", "lp");
      leg1->Draw("same");
      // ========== 10. 保存图片 ==========
      outFile->cd();
      c->Write(Form("c_Spin_correlation_P_vs_Lambda_pairs_Rbin_%.2f_%.2f", massBinCenter[i],massBinCenter[i+1]));
      c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs_Rbin_%.2f_%.2f.png", massBinCenter[i],massBinCenter[i+1]));


    }
    outFile->Close();
    

    
}


void draw_STAR_x_P_Lambda_y(){
  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs.root"),"recreate");

  const int nMassBins1 = 6;
  
  double massBinCenter[nMassBins1+1] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};

    // ========== 1. 打开文件，读取直方图 ==========
    TFile* f_fit_fixed = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root","READ");
    if (!f_fit_fixed || f_fit_fixed->IsZombie()) {
        cout << "ERROR: 文件打开失败！" << endl;
        return;
    }
    TFile* f_fit_fixed3 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Sourav_analysis/Result_signalType0_bkgType2_fixedmusigma_1/output_corr_vsMass_cent00.root","READ");


    
    
    TH1D* hfit_unlike = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_US_vsR;1");   // ΛΛ̅ 异号
    TH1D* hfit_like_Lambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_LL_vsR;1");   // ΛΛ 同号
    TH1D* hfit_like_AntiLambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_ALAL_vsR;1");   // Λ̅ Λ̅  同号

    TH1D* hfit_unlike3 = (TH1D*)f_fit_fixed3->Get("hfit_Unlike_vsR;1");   // ΛΛ̅ 异号
    TH1D* hfit_like_Lambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle0_vsR;1");   // ΛΛ 同号
    TH1D* hfit_like_AntiLambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle3_vsR;1");   // Λ̅ Λ̅  同号
   
    for(int i = 0;i<=nMassBins1-1;i++){
      TCanvas *c = new TCanvas("c", Form("Spin-spin correlation: P vs #Lambda pairs_%d",i), 1800, 1200);

      
        
      int nMassBins = hfit_unlike->GetNbinsX();

    
      

      TPad *pad1 = new TPad("pad1", "", 0.01, 0.01, 0.99, 0.99);
      pad1->SetLeftMargin(0.18);
      pad1->SetBottomMargin(0.15);
      pad1->SetTopMargin(0.1);
      pad1->SetRightMargin(0.05);
      pad1->SetGrid(0, 0);
      pad1->Draw();
      pad1->cd();

      
      TH2F *hframe1 = new TH2F("hframe1", "", 100, -0.25, 0.39, 3, 0.5, 3.5);
      hframe1->SetStats(0);
      hframe1->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
      hframe1->GetXaxis()->SetTitleOffset(1.4);
      hframe1->GetXaxis()->SetRangeUser(-0.04, 0.04);

    
      hframe1->GetYaxis()->SetBinLabel(1, "#Lambda#bar{#Lambda}");
      hframe1->GetYaxis()->SetBinLabel(2, "#Lambda#Lambda");
      hframe1->GetYaxis()->SetBinLabel(3, "#bar{#Lambda}#bar{#Lambda}");
      hframe1->GetYaxis()->SetLabelSize(0.08);
      hframe1->GetYaxis()->SetLabelOffset(0.02);
      hframe1->GetYaxis()->SetTitle("");
      hframe1->Draw();

    
      TLine *line0 = new TLine(0, 0.5, 0, 3.5);
      line0->SetLineStyle(2);
      line0->SetLineWidth(4);  
      line0->SetLineColor(kBlack);
      line0->Draw("same");

    
  
      // 异号 ΛΛ̅
      double p_unlike    = hfit_unlike->GetBinContent(i+1);
      double err_unlike  = hfit_unlike->GetBinError(i+1);

      // 同号 ΛΛ
      double p_like     = hfit_like_Lambda->GetBinContent(i+1);
      double err_like   = hfit_like_Lambda->GetBinError(i+1);

      // 同号 Λ̅Λ̅
      double p_like_AntiLambda = hfit_like_AntiLambda->GetBinContent(i+1);
      double err_like_AntiLambda = hfit_like_AntiLambda->GetBinError(i+1);

      /////////////
      double p_unlike3    = hfit_unlike3->GetBinContent(i+1);
      double err_unlike3  = hfit_unlike3->GetBinError(i+1);

      // 同号 ΛΛ
      double p_like3     = hfit_like_Lambda3->GetBinContent(i+1);
      double err_like3  = hfit_like_Lambda3->GetBinError(i+1);

      // 同号 Λ̅Λ̅
      double p_like_AntiLambda3 = hfit_like_AntiLambda3->GetBinContent(i+1);
      double err_like_AntiLambda3 = hfit_like_AntiLambda3->GetBinError(i+1);

    
      const int n = 3;
      double x[n]  = {p_unlike,  p_like, p_like_AntiLambda };   // X = P 值
      double y[n]  = {1.0,       2.0, 3.0 };        // Y = 类别位置
      double ex[n] = {err_unlike, err_like, err_like_AntiLambda};// 误差棒
      double ey[n] = {0.0, 0.0,0.0};

      double x3[n]  = {p_unlike3,  p_like3, p_like_AntiLambda3 };   // X = P 值
      double y3[n]  = {1.0,       2.0, 3.0 };        // Y = 类别位置
      double ex3[n] = {err_unlike3, err_like3, err_like_AntiLambda3};// 误差棒
      double ey3[n] = {0.0, 0.0,0.0};



      // ========== 8. 画带误差的点 ==========
      TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(2.0);
      gr->SetMarkerColor(kRed);
      gr->SetLineColor(kRed);
      gr->Draw("P same");

      TGraphErrors *gr3 = new TGraphErrors(n, x3, y3, ex3, ey3);
      gr3->SetMarkerStyle(20);
      gr3->SetMarkerSize(2.0);
      gr3->SetMarkerColor(kBlue);
      gr3->SetLineColor(kBlue);
      //gr3->Draw("P same");

    

      // ========== 9. 文字标注 ==========
      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextSize(0.04);
      tex->DrawLatex(0.20, 0.85, "ALICE");
      tex->DrawLatex(0.20, 0.80, Form("#DeltaR =#sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}: [%.2f-%.2f]", massBinCenter[i],massBinCenter[i+1]));   
      tex->DrawLatex(0.20, 0.75, "p+p #sqrt{s} = 13.6 TeV");
      tex->DrawLatex(0.20, 0.70, "|y| < 1.5");
      tex->DrawLatex(0.20, 0.65, "Helicity frame");
      TLegend *leg = new TLegend(0.75, 0.80, 0.82, 0.95);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gr, "Data", "lp");
      //leg->AddEntry(gr3, "2D fit method,Double Gaussian+Cheby2:fixed #mu,#sigma", "lp");
      leg->Draw("same");
      // ========== 10. 保存图片 ==========
      outFile->cd();
      c->Write(Form("c_Spin_correlation_P_vs_Lambda_pairs_Rbin_%.2f_%.2f", massBinCenter[i],massBinCenter[i+1]));
      c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs_Rbin_%.2f_%.2f.png", massBinCenter[i],massBinCenter[i+1]));


    }
    outFile->Close();
    

    
}

void test_draw_STAR_x_P_Lambda_y(){
  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs.root"),"recreate");

  const int nMassBins1 = 6;
  double massBinCenter[nMassBins1+1] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};

  TFile* f_fit_fixed = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root","READ");
  if (!f_fit_fixed || f_fit_fixed->IsZombie()) {
    cout << "ERROR: 文件打开失败！" << endl;
    return;
  }
  TFile* f_fit_fixed3 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Sourav_analysis/Result_signalType0_bkgType2_fixedmusigma_1/output_corr_vsMass_cent00.root","READ");

  TH1D* hfit_unlike = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_US_vsR;1");
  TH1D* hfit_like_Lambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_LL_vsR;1");
  TH1D* hfit_like_AntiLambda  = (TH1D*)f_fit_fixed->Get("sig0_bkg2/hP_ALAL_vsR;1");

  TH1D* hfit_unlike3 = (TH1D*)f_fit_fixed3->Get("hfit_Unlike_vsR;1");
  TH1D* hfit_like_Lambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle0_vsR;1");
  TH1D* hfit_like_AntiLambda3  = (TH1D*)f_fit_fixed3->Get("hfit_particle3_vsR;1");

  TFile* RMShistFile = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "READ");
  TH1D* h_RMS_US = (TH1D*)RMShistFile->Get("h_RMS_US");
  TH1D* h_RMS_SS = (TH1D*)RMShistFile->Get("h_RMS_SS");


  TCanvas *c = new TCanvas("c", "Spin-spin correlation: P vs #Lambda pairs", 2000, 1400);
  c->Divide(3, 2, 0, 0);

  for(int i = 0; i < nMassBins1; i++){
    c->cd(i+1);

    TPad *pad = new TPad(Form("pad_%d",i), "", 0, 0, 1, 1);

    if(i == 0 || i == 3) pad->SetLeftMargin(0.18);
    else                 pad->SetLeftMargin(0.0);

    if(i == 2 || i == 5) pad->SetRightMargin(0.05);
    else                 pad->SetRightMargin(0.0);

    if(i < 3){
        pad->SetTopMargin(0.08);
        pad->SetBottomMargin(0.0);
    } else {
        pad->SetTopMargin(0.0);
        pad->SetBottomMargin(0.15);
    }

    pad->Draw();
    pad->cd();

    TH2F *hframe1 = new TH2F(Form("hframe1_%d",i), "", 100, -0.25, 0.39, 3, 0.5, 3.5);
    hframe1->SetStats(0);
    hframe1->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
    hframe1->GetXaxis()->SetTitleOffset(1.2);
    hframe1->GetXaxis()->SetRangeUser(-0.04, 0.04);
    hframe1->GetXaxis()->CenterTitle();
    hframe1->GetXaxis()->SetTitleSize(0.05);
    //gPad->UseCurrentStyle();
    gPad->SetFrameLineWidth(3); 

    if(i < 3){
        hframe1->GetXaxis()->SetLabelSize(0);
        //hframe1->GetXaxis()->SetTickLength(0);
    }

    if(i == 0 || i ==3){
      hframe1->GetYaxis()->SetBinLabel(1, "#Lambda#bar{#Lambda}");
      hframe1->GetYaxis()->SetBinLabel(2, "#Lambda#Lambda");
      hframe1->GetYaxis()->SetBinLabel(3, "#bar{#Lambda}#bar{#Lambda}");
      hframe1->GetYaxis()->SetLabelSize(0.075);
      hframe1->GetYaxis()->SetLabelOffset(0.02);
    } else {
      hframe1->GetYaxis()->SetBinLabel(1, "");
      hframe1->GetYaxis()->SetBinLabel(2, "");
      hframe1->GetYaxis()->SetBinLabel(3, "");
      hframe1->GetYaxis()->SetLabelSize(0);
    }

    hframe1->GetYaxis()->SetTitle("");
    hframe1->Draw();

    // ====================== 【正确：顶部居中大标题】 ======================
    if(i == 1){
      TLatex *title = new TLatex();
      title->SetNDC();
      title->SetTextSize(0.06);    // 字号可调：0.08 很大，0.07 中等
      title->SetTextAlign(22);     // 居中
      title->DrawLatex(0.5, 0.96, "Spin correlation P_{#Lambda_{1}#Lambda_{2}} vs #DeltaR");
    }

    TLine *line0 = new TLine(0, 0.5, 0, 3.5);
    line0->SetLineStyle(2);
    line0->SetLineWidth(3);
    line0->SetLineColor(kBlack);
    line0->Draw("same");

    double p_unlike    = hfit_unlike->GetBinContent(i+1);
    double err_unlike  = hfit_unlike->GetBinError(i+1);
    double p_like      = hfit_like_Lambda->GetBinContent(i+1);
    double err_like    = hfit_like_Lambda->GetBinError(i+1);
    double p_like_AntiLambda = hfit_like_AntiLambda->GetBinContent(i+1);
    double err_like_AntiLambda = hfit_like_AntiLambda->GetBinError(i+1);

    const int n = 3;
    double x[n]  = {p_unlike,  p_like, p_like_AntiLambda };
    double y[n]  = {1.0,       2.0, 3.0 };
    double ex[n] = {err_unlike, err_like, err_like_AntiLambda};
    double ey[n] = {0.0, 0.0, 0.0};

    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.8);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->Draw("P same");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.04);
   
    if(i == 0){
      tex->DrawLatex(0.2, 0.84, "#DeltaR=#sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}");
      tex->DrawLatex(0.66, 0.80, "ALICE");
      tex->DrawLatex(0.66, 0.72, "13.6 TeV pp");
      tex->DrawLatex(0.66, 0.68, "Helicity frame");

      TLegend *leg = new TLegend(0.70, 0.82, 0.90, 0.95);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(gr, "Data", "lp");
      leg->Draw("same");
    }
    tex->DrawLatex(0.66, 0.76, Form("#DeltaR: [%.2f-%.2f]", massBinCenter[i],massBinCenter[i+1]));
  }

  outFile->cd();
  c->Write("c_Spin_correlation_P_vs_Lambda_pairs_2x3");
  c->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Spin_correlation_P_vs_Lambda_pairs_2x3.png");
  outFile->Close();
}


void Data_Eventnumber(){
  TFile* f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/LHC23_pass4_thin_675285_AnalysisResult/merged_AnalysisResults_pass4_thin_675285.root");
  TH1D* hEventNum = (TH1D*)f->Get("lambdaspincorrelation/hEvtSelInfo");
  
  TCanvas* c = new TCanvas("c", "Event Number", 1600, 1200);
  c->SetLogy(); // 保持y轴对数刻度
  hEventNum->SetLineColor(kBlue);
  hEventNum->SetLineWidth(2);
  hEventNum->SetTitle("Event Selection Info;Selection Step;Number of Events");
  hEventNum->GetYaxis()->SetTitle("Number of Events");
  // 设置bin标签
  hEventNum->GetXaxis()->SetBinLabel(1, "|V_{z}|<10cm");
  hEventNum->GetXaxis()->SetBinLabel(2, "Other selections");
  hEventNum->GetXaxis()->SetBinLabel(3, "Number of V^{0} > 1");
  hEventNum->SetStats(0);
  hEventNum->Draw("HIST"); // 用HIST绘制确保bin边界清晰
  
  // ========== 核心新增：在每个bin上方显示content数值 ==========
  TLatex* latex = new TLatex();
  latex->SetTextFont(42);    // 字体（42=Helvetica，兼容ROOT默认）
  latex->SetTextSize(0.025); // 字体大小（适配1600x1200画布）
  latex->SetTextColor(kBlack); // 字体颜色（对比蓝色bin更清晰）
  latex->SetNDC(false);      // 使用像素坐标（而非归一化坐标）

  // 遍历3个bin，逐个绘制数值
  for (int i = 1; i <= 3; i++) {
    double binContent = hEventNum->GetBinContent(i); // 获取bin的content
    double binX = hEventNum->GetXaxis()->GetBinCenter(i); // bin的X中心坐标
    double binY = binContent * 1.2; // 数值显示在bin内容的1.2倍高度（避免重叠）
    
    // 格式化数值：如果是大数，用科学计数法；小数直接显示
    TString label;
    if (binContent > 1e4) {
      label.Form("%.8e", binContent); // 科学计数法（保留2位小数）
    } else {
      label.Form("%.0f", binContent); // 整数显示
    }
    
    // 在bin中心位置绘制数值
    latex->DrawLatex(binX-hEventNum->GetXaxis()->GetBinWidth(i)/2, binY, label);
  }

  c->SaveAs("Picture/EventNumber.png"); // 按需取消注释保存图片
}


void MC_Eventnumber(){
  TFile* f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysisresultfile/merged_AnalysisResults.root");
  TH1D* hEventNum = (TH1D*)f->Get("lambdaspincorrelation/hEvtSelInfo");
  
  TCanvas* c = new TCanvas("c", "Event Number", 1600, 1200);
  c->SetLogy(); // 保持y轴对数刻度
  hEventNum->SetLineColor(kBlue);
  hEventNum->SetLineWidth(2);
  hEventNum->SetTitle("Event Selection Info;Selection Step;Number of Events");
  hEventNum->GetYaxis()->SetTitle("Number of Events");
  // 设置bin标签
  hEventNum->GetXaxis()->SetBinLabel(1, "|V_{z}|<10cm");
  hEventNum->GetXaxis()->SetBinLabel(2, "Other selections");
  hEventNum->GetXaxis()->SetBinLabel(3, "Number of V^{0} > 1");
  hEventNum->SetStats(0);
  hEventNum->Draw("HIST"); // 用HIST绘制确保bin边界清晰
  
  // ========== 核心新增：在每个bin上方显示content数值 ==========
  TLatex* latex = new TLatex();
  latex->SetTextFont(42);    // 字体（42=Helvetica，兼容ROOT默认）
  latex->SetTextSize(0.025); // 字体大小（适配1600x1200画布）
  latex->SetTextColor(kBlack); // 字体颜色（对比蓝色bin更清晰）
  latex->SetNDC(false);      // 使用像素坐标（而非归一化坐标）

  // 遍历3个bin，逐个绘制数值
  for (int i = 1; i <= 3; i++) {
    double binContent = hEventNum->GetBinContent(i); // 获取bin的content
    double binX = hEventNum->GetXaxis()->GetBinCenter(i); // bin的X中心坐标
    double binY = binContent * 1.2; // 数值显示在bin内容的1.2倍高度（避免重叠）
    
    // 格式化数值：如果是大数，用科学计数法；小数直接显示
    TString label;
    if (binContent > 1e4) {
      label.Form("%.8e", binContent); // 科学计数法（保留2位小数）
    } else {
      label.Form("%.0f", binContent); // 整数显示
    }
    
    // 在bin中心位置绘制数值
    latex->DrawLatex(binX-hEventNum->GetXaxis()->GetBinWidth(i)/2, binY, label);
  }

  c->SaveAs("Picture/EventNumber.png"); // 按需取消注释保存图片
}

void plot_thn_mass2D()
{
    TCanvas *c = new TCanvas("c", "Mass vs Mass 2D", 2500, 1200);
    c->Divide(3,2);
    
    // 1. 打开文件
    TFile *f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_setting5.root");
    if (!f){cout<<"无法打开文件！"<<endl;} 

    // 2. 读取你的 THnSparseF
    THnSparseF *hn = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaAntiLambda;1");
    if (!hn) {
        cout << "找不到 hSparseLambdaAntiLambda" << endl;
        return;
    }

    // 3. 投影：轴0(质量) + 轴1(质量) → 二维质量图
    TH2D *h2d_mass = hn->Projection(0, 1);  // <-- 关键！

    // 4. 画图
   
    h2d_mass->SetTitle("2D #Lambda #bar{#Lambda} invariant mass;m_{#Lambda} (GeV);m_{#bar{#Lambda}} (GeV)");
    c->cd(1);
    h2d_mass->GetXaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");
    h2d_mass->GetYaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");

    h2d_mass->Draw("COLZ");  // 彩色热图
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.15, 0.85, "ALICE");
    tex->DrawLatex(0.15, 0.80, "0.5 < p_{T} < 3 GeV/c");
    tex->DrawLatex(0.55, 0.85, "p+p MB #sqrt{s} = 13.6 TeV");
    tex->DrawLatex(0.15, 0.75, "|y| < 0.5");
    //tex->DrawLatex(0.15, 0.70, "Helicity frame");

    // 5. 美化
    h2d_mass->SetStats(0);
    gPad->SetLogz(1); // 打开对数Z轴，看得更清楚
    c->cd(4);
    TH2D *h2d_mass_clone = (TH2D*)h2d_mass->Clone("h2d_mass_clone");
    gPad->SetLogz(1); // 打开对数Z轴，看得更清楚
    h2d_mass_clone->GetXaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");
    h2d_mass_clone->GetYaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");
    h2d_mass_clone->GetXaxis()->SetTitleOffset(2);
    h2d_mass_clone->GetYaxis()->SetTitleOffset(2);

    h2d_mass_clone->Draw("LEGO2");

    c->cd(2);

    THnSparseF *hnLL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaLambda;1");
    if (!hn) {
        cout << "找不到 hSparseLambdaAntiLambda" << endl;
        return;
    }

    // 3. 投影：轴0(质量) + 轴1(质量) → 二维质量图
    TH2D *h2d_mass_LL = hnLL->Projection(0, 1);  // <-- 关键！

    // 4. 画图
   
    h2d_mass_LL->SetTitle("2D #Lambda #Lambda invariant mass;m_{#Lambda} (GeV);m_{#Lambda} (GeV)");
   
    h2d_mass_LL->GetXaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");
    h2d_mass_LL->GetYaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");

    h2d_mass_LL->Draw("COLZ");  // 彩色热图
    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.04);
    tex2->DrawLatex(0.15, 0.85, "ALICE");
    tex2->DrawLatex(0.15, 0.80, "0.5 < p_{T} < 3 GeV/c");
    tex2->DrawLatex(0.55, 0.85, "p+p MB #sqrt{s} = 13.6 TeV");
    tex2->DrawLatex(0.15, 0.75, "|y| < 0.5");
    //tex2->DrawLatex(0.15, 0.70, "Helicity frame");
    gPad->SetLogz(1); 

    // 5. 美化
    h2d_mass_LL->SetStats(0);

    c->cd(5);

    TH2D *h2d_mass_clone_LL = (TH2D*)h2d_mass_LL->Clone("h2d_mass_clone");
    gPad->SetLogz(1); // 打开对数Z轴，看得更清楚
    h2d_mass_clone_LL->GetXaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");
    h2d_mass_clone_LL->GetYaxis()->SetTitle("M_{inv}(p#pi^{-}) GeV/c^{2}");
    h2d_mass_clone_LL->GetXaxis()->SetTitleOffset(2);
    h2d_mass_clone_LL->GetYaxis()->SetTitleOffset(2);

    h2d_mass_clone_LL->Draw("LEGO2");

    c->cd(3);
    THnSparseF *hnALAL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseAntiLambdaAntiLambda;1");
    if (!hn) {
        cout << "找不到 hSparseLambdaAntiLambda" << endl;
        return;
    } 
    // 3. 投影：轴0(质量) + 轴1(质量) → 二维质量图
    TH2D *h2d_mass_ALAL = hnALAL->Projection(0, 1);  // <-- 关键！
    // 4. 画图
    h2d_mass_ALAL->SetTitle("2D #bar{#Lambda} #bar{#Lambda} invariant mass;m_{#bar{#Lambda}} (GeV);m_{#bar{#Lambda}} (GeV)");
    h2d_mass_ALAL->GetXaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");
    h2d_mass_ALAL->GetYaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");
    h2d_mass_ALAL->Draw("COLZ");  // 彩色热图
    TLatex *tex3 = new TLatex();
    tex3->SetNDC();
    tex3->SetTextSize(0.04);
    tex3->DrawLatex(0.15, 0.85, "ALICE");
    tex3->DrawLatex(0.15, 0.80, "0.5 < p_{T} < 3 GeV/c");
    tex3->DrawLatex(0.55, 0.85, "p+p MB #sqrt{s} = 13.6 TeV");
    tex3->DrawLatex(0.15, 0.75, "|y| < 0.5");
    //tex3->DrawLatex(0.15, 0.70, "Helicity frame");
    gPad->SetLogz(1);
    // 5. 美化
    h2d_mass_ALAL->SetStats(0);
    c->cd(6);
    TH2D *h2d_mass_clone_ALAL = (TH2D*)h2d_mass_ALAL->Clone("h2d_mass_clone");
    gPad->SetLogz(1); // 打开对数Z轴，看得更清楚
    h2d_mass_clone_ALAL->GetXaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");
    h2d_mass_clone_ALAL->GetYaxis()->SetTitle("M_{inv}(#bar{p}#pi^{+}) GeV/c^{2}");
    h2d_mass_clone_ALAL->GetXaxis()->SetTitleOffset(2);
    h2d_mass_clone_ALAL->GetYaxis()->SetTitleOffset(2);
    h2d_mass_clone_ALAL->Draw("LEGO2");
    




    c->SaveAs("Picture/mass_2D_All_from_thnsparse.pdf");
}


void SEOverMEForBkgwitSandBandMethod(bool SS = true){
  

   //draw SE distribution for US
  TString SSandUS = SS ?"SS":"US";
  TString SSorUS  = SS ?"Like-sign(#Lambda#Lambda+#bar{#Lambda}#bar{#Lambda})":"Unlike-sign(#Lambda#bar{#Lambda}+#bar{#Lambda}#Lambda)";
  

  vector<string> methods = {
    "YieldSS",//0
    "do_YieldSS_BCD"//1
  };

  vector<string> configs = {
    "sig0_bkg0",//0
    "sig0_bkg1",//1
    "sig0_bkg2"   // <-- default 2
    //"sig2_bkg0",
    //"sig2_bkg1",
    //"sig2_bkg2"
  };
  
  gSystem->mkdir(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_Bkg_spin_with_SandBandMethod/"), kTRUE);
  TFile* outFile = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_Bkg_spin_with_SandBandMethod/%s_SEOverME.root",SSandUS.Data()),"recreate");

  TFile* Clone_f_default = nullptr;

  for(int Rbin = 0; Rbin <= 5; ++Rbin){
    TCanvas* c = new TCanvas("c", Form("c_Rbin%d", Rbin), 3400,3400);
    c->Divide(3,3,0);
      
    string prefix = "/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_";
    string suffix = "_Rdep_Output_pp_";
    string f_path = prefix + methods[1] + suffix + configs[2] + ".root";
    TFile* f_default = TFile::Open(f_path.c_str(), "READ");
    TH1D* hSE_B_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hSE_%s_B_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    TH1D* hME_B_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hME_%s_B_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    //hSE_US_sig0_bkg2->GetYaxis()->SetRangeUser(0,hSE_US_sig0_bkg2->GetMaximum()*1.5);
    Clone_f_default = f_default;
   



    hSE_B_sig0_bkg2->Scale(1.0/hSE_B_sig0_bkg2->Integral());
    hME_B_sig0_bkg2->Scale(1.0/hME_B_sig0_bkg2->Integral());

    


    TH1D* hRatio_SE_ME = (TH1D*)hSE_B_sig0_bkg2->Clone("hRatio_SE_ME");
    hRatio_SE_ME->Divide(hME_B_sig0_bkg2);

    c->cd(1);
    gPad->SetBottomMargin(0.0);
    gPad->SetRightMargin(0.00);
    
     

    hSE_B_sig0_bkg2->GetYaxis()->SetRangeUser(0.,0.25);
    hSE_B_sig0_bkg2->SetTitle("");
    hSE_B_sig0_bkg2->GetYaxis()->SetTitle("1/N dN/dcos#theta*");
    hSE_B_sig0_bkg2->SetLineColor(kRed);
    hSE_B_sig0_bkg2->SetMarkerColor(kRed);
    hSE_B_sig0_bkg2->SetLineWidth(2);
    hME_B_sig0_bkg2->SetLineWidth(2);
    hSE_B_sig0_bkg2->GetYaxis()->SetTitleOffset(1.4);
  
    //hSE_US_sig0_bkg2->SetMarkerStyle(20);
    //hME_US_sig0_bkg2->SetMarkerStyle(21);
    hSE_B_sig0_bkg2->Draw();
    hME_B_sig0_bkg2->Draw("SAME");
 
    const double REdges[] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};
    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextSize(0.03);
    title->SetTextAlign(22);
    title->DrawLatex(0.3, 0.95, Form("ALICE"));
    title->DrawLatex(0.3, 0.85, Form("0.8 GeV/c < p_{T} < 3 GeV/c"));
    title->DrawLatex(0.3, 0.75, Form("|y| < 0.5"));
    title->DrawLatex(0.3, 0.65, Form("Rbin %d: %.1f < R < %.1f", Rbin, REdges[Rbin], REdges[Rbin+1]));
    title->DrawLatex(0.3, 0.55, Form("Helicity frame"));

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(hSE_B_sig0_bkg2, Form("SE_B_%s", SSorUS.Data()), "lp");
    legend->AddEntry(hME_B_sig0_bkg2, Form("ME_B_%s", SSorUS.Data()), "lp");
    legend->SetTextSize(0.05);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();

    c->cd(2);

    gPad->SetBottomMargin(0.0);
    gPad->SetRightMargin(0.00);
    gPad->SetLeftMargin(0.0);
    TH1D* hSE_C_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hSE_%s_C_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    TH1D* hME_C_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hME_%s_C_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
   
    hSE_C_sig0_bkg2->Scale(1.0/hSE_C_sig0_bkg2->Integral());
    hME_C_sig0_bkg2->Scale(1.0/hME_C_sig0_bkg2->Integral());

    hSE_C_sig0_bkg2->GetYaxis()->SetRangeUser(0.,0.25);
    hSE_C_sig0_bkg2->SetTitle("");
    hSE_C_sig0_bkg2->GetYaxis()->SetTitle("1/N dN/dcos#theta*");
    hSE_C_sig0_bkg2->SetLineColor(kRed);
    hSE_C_sig0_bkg2->SetMarkerColor(kRed);
    hSE_C_sig0_bkg2->SetLineWidth(2);
    hME_C_sig0_bkg2->SetLineWidth(2);

    hSE_C_sig0_bkg2->Draw();
    hME_C_sig0_bkg2->Draw("SAME");

    TLegend* legend2 = new TLegend(0.4, 0.7, 0.8, 0.9);
    legend2->AddEntry(hSE_C_sig0_bkg2, Form("SE_C_%s", SSorUS.Data()), "lp");
    legend2->AddEntry(hME_C_sig0_bkg2, Form("ME_C_%s", SSorUS.Data()), "lp");
    legend2->SetTextSize(0.05);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->Draw();

    c->cd(3);

    gPad->SetBottomMargin(0.0);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.0);
    TH1D* hSE_D_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hSE_%s_D_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    TH1D* hME_D_sig0_bkg2 = (TH1D*)f_default->Get(Form("%s/rbin%d/hME_%s_D_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
   
    hSE_D_sig0_bkg2->Scale(1.0/hSE_D_sig0_bkg2->Integral());
    hME_D_sig0_bkg2->Scale(1.0/hME_D_sig0_bkg2->Integral());

    hSE_D_sig0_bkg2->GetYaxis()->SetRangeUser(0.,0.25);
    hSE_D_sig0_bkg2->SetTitle("");
    hSE_D_sig0_bkg2->GetYaxis()->SetTitle("1/N dN/dcos#theta*");
    hSE_D_sig0_bkg2->SetLineColor(kRed);
    hSE_D_sig0_bkg2->SetMarkerColor(kRed);
    hSE_D_sig0_bkg2->SetLineWidth(2);
    hME_D_sig0_bkg2->SetLineWidth(2);

    

    hSE_D_sig0_bkg2->Draw();
    hME_D_sig0_bkg2->Draw("SAME");

    TLegend* legend3 = new TLegend(0.4, 0.7, 0.8, 0.9);
    legend3->AddEntry(hSE_D_sig0_bkg2, Form("SE_D_%s", SSorUS.Data()), "lp");
    legend3->AddEntry(hME_D_sig0_bkg2, Form("ME_D_%s", SSorUS.Data()), "lp");
    legend3->SetTextSize(0.05);
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->Draw();


    c->cd(4);

    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.0);
     gPad->SetBottomMargin(0.1);
    TH1D* hRatio_SE_ME_B =  (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_B_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    
    //hRatio_SE_ME_B->GetYaxis()->SetRangeUser(0.95,1.05);
    hRatio_SE_ME_B->GetYaxis()->SetTitle("SE/ME");
    hRatio_SE_ME_B->SetTitle("");

    hRatio_SE_ME_B->GetYaxis()->SetTitleOffset(1.4);
    hRatio_SE_ME_B->Draw();

    c->cd(5);

    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.0);
    gPad->SetLeftMargin(0.0);
     gPad->SetBottomMargin(0.1);
    TH1D* hRatio_SE_ME_C =  (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_C_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    
    //hRatio_SE_ME_C->GetYaxis()->SetRangeUser(0.95,1.05);
    hRatio_SE_ME_C->SetTitle("");

    hRatio_SE_ME_C->Draw();

    c->cd(6);

    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.0);
    gPad->SetBottomMargin(0.1);
    TH1D* hRatio_SE_ME_D =  (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_D_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    
    //hRatio_SE_ME_D->GetYaxis()->SetRangeUser(0.95,1.05);
    hRatio_SE_ME_D->SetTitle("");

    hRatio_SE_ME_D->Draw();

    c->cd(7);
    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.01);
    //gPad->SetLeftMargin(0.0);
    TH1D* hCorr_US_BPlusC = (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_BPlusC_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    hCorr_US_BPlusC->GetYaxis()->SetTitle("SE/ME");
    hCorr_US_BPlusC->SetTitle("B+C");
    //hCorr_US_BPlusC->GetYaxis()->SetRangeUser(1.9,2.1);
    hCorr_US_BPlusC->Draw();

    c->cd(8);
    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.1);
    TH1D* hCorr_US_BPlusCMinusD = (TH1D*)f_default->Get(Form("%s/rbin%d/hCorr_%s_BPlusCMinusD_sig0_bkg2_rbin%d",configs[2].c_str(), Rbin, SSandUS.Data(), Rbin));
    hCorr_US_BPlusCMinusD->GetYaxis()->SetTitle("SE/ME");
    hCorr_US_BPlusCMinusD->SetTitle("B+C-D");
    //hCorr_US_BPlusCMinusD->GetYaxis()->SetRangeUser(0.9,1.1);
    hCorr_US_BPlusCMinusD->Draw();

    c->cd(9);

    gPad->SetTopMargin(0.0);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.1);

    // Divide the 9th pad into two vertically
    TPad *pad9_top = new TPad("pad9_top", "top", 0.0, 0.5, 1.0, 1.0);
    TPad *pad9_bottom = new TPad("pad9_bottom", "bottom", 0.0, 0.0, 1.0, 0.5);
    pad9_top->SetBottomMargin(0.0);
    pad9_top->SetTopMargin(0.1);
    pad9_bottom->SetTopMargin(0.0);
    pad9_bottom->SetBottomMargin(0.15);

    pad9_top->Draw();
    pad9_bottom->Draw();
    
    pad9_top->cd();

    TH1D* hP_US_BPlusC = (TH1D*)f_default->Get(Form("%s/hP_%s_BPlusC_vsR",configs[2].c_str(), SSandUS.Data()));
    hP_US_BPlusC->Draw();
    hP_US_BPlusC->SetTitle("spin correction from B+C ");
    //hP_US_BPlusC->GetYaxis()->SetRangeUser(-0.02,0.02);
    
    hP_US_BPlusC->GetYaxis()->SetTitleSize(0.05);
    hP_US_BPlusC->GetYaxis()->SetLabelSize(0.05);
    
    

    pad9_bottom->cd();

    TH1D* hP_US_BPlusCMinusD = (TH1D*)f_default->Get(Form("%s/hP_%s_BPlusCMinusD_vsR",configs[2].c_str(), SSandUS.Data()));
    hP_US_BPlusCMinusD->SetTitle("spin correction from B+C-D");
    hP_US_BPlusCMinusD->GetYaxis()->SetTitleSize(0.05);
    hP_US_BPlusCMinusD->GetXaxis()->SetTitleSize(0.05);
   
    hP_US_BPlusCMinusD->GetYaxis()->SetLabelSize(0.05);
    //hP_US_BPlusCMinusD->GetYaxis()->SetRangeUser(-0.02,0.02);
    hP_US_BPlusCMinusD->Draw();


    
   

    
    
    outFile->cd();
    outFile->Write();
    c->Write(Form("%s_SE_ME_B_and_Corr_Rbin%d", SSandUS.Data(), Rbin));
    
    c->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_Bkg_spin_with_SandBandMethod/%s_SE_ME_B_and_Corr_Rbin%d.png", SSandUS.Data(), Rbin));
    delete c;

  }
 

  TCanvas* c1 = new TCanvas("c1", "c1", 1600, 1200);
  TH1D* hP_US_BPlusC = (TH1D*)Clone_f_default->Get(Form("%s/hP_%s_BPlusC_vsR",configs[2].c_str(), SSandUS.Data()));
  TH1D* hP_US_BPlusCMinusD = (TH1D*)Clone_f_default->Get(Form("%s/hP_%s_BPlusCMinusD_vsR",configs[2].c_str(), SSandUS.Data()));

  TH1D* hP_BCD_US_BB = (TH1D*)Clone_f_default->Get(Form("%s/hP_BCD_%s_BB_vsR;1",configs[2].c_str(), SSandUS.Data()));



  hP_US_BPlusC->SetLineColor(kBlue);
  hP_US_BPlusC->SetLineWidth(4);
  hP_US_BPlusCMinusD->SetLineColor(kRed);
  hP_US_BPlusCMinusD->SetLineWidth(4);
  hP_BCD_US_BB->SetLineColor(kGreen);
  hP_BCD_US_BB->SetLineWidth(4);

  hP_US_BPlusC->GetYaxis()->SetRangeUser(-0.02,0.02);
  hP_US_BPlusC->GetYaxis()->SetTitleOffset(0.7);

  hP_US_BPlusC->Draw();
  hP_US_BPlusCMinusD->Draw("same");
  hP_BCD_US_BB->Draw("same");

  TLegend* leg = new TLegend(0.4, 0.7, 0.9, 0.9);
  leg->AddEntry(hP_US_BPlusC, Form("%s_B+C:#frac{SE(B)}{ME(B)}+#frac{SE(C)}{ME(C)}", SSandUS.Data()), "l");
  leg->AddEntry(hP_US_BPlusCMinusD, Form("%s_B+C-D:#frac{SE(B)}{ME(B)}+#frac{SE(C)}{ME(C)}-#frac{SE(D)}{ME(D)}", SSandUS.Data()), "l");
  leg->AddEntry(hP_BCD_US_BB, Form("%s_B+C-D:#frac{SE(B+C-D)}{ME(B+C-D)}", SSandUS.Data()), "l");
  leg->SetTextSize(0.02);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c1->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SEOverME_Bkg_spin_with_SandBandMethod/%s_Spin_Corr_BPlusC_vs_R.png", SSandUS.Data()));

  outFile->cd();
  outFile->Write();
  c1->Write(Form("hP_%s_vsR", SSandUS.Data()));

  delete c1;

  outFile->Close();

}

void CalculationOfSpinCorrelation(){

  TFile* outFile = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AddsignalFrac/Spincore_withFrac.root", "RECREATE");

  // ============================================
  // 1. 加载数据
  // ============================================
  
  // US (Unlike-Sign: Lambda + AntiLambda) 的 fs
  TFile* fUS = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Output_pp_sig0_bkg2_hSparseLambdaAntiLambda/fs_avg.root");
  if (!fUS || fUS->IsZombie()) {
    std::cerr << "ERROR: Cannot open US fs_avg.root" << std::endl;
    return;
  }
  TH1D* hFs_US = (TH1D*)fUS->Get("hFsAvg_vs_R");
  if (!hFs_US) {
    std::cerr << "ERROR: Cannot find hFsAvg_vs_R in US file" << std::endl;
    return;
  }

  // SS (Like-Sign: LambdaLambda + AntiLambdaAntiLambda) 的 fs
  TFile* fSS = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Output_pp_sig0_bkg2_hSparseLambdaLambda/fs_avg.root");
  if (!fSS || fSS->IsZombie()) {
    std::cerr << "ERROR: Cannot open SS fs_avg.root" << std::endl;
    return;
  }
  TH1D* hFs_SS = (TH1D*)fSS->Get("hFsAvg_vs_R");
  if (!hFs_SS) {
    std::cerr << "ERROR: Cannot find hFsAvg_vs_R in SS file" << std::endl;
    return;
  }

  // 结果文件
  TFile* fResults = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root");
  if (!fResults || fResults->IsZombie()) {
    std::cerr << "ERROR: Cannot open Results file" << std::endl;
    return;
  }

  // US: 信号+背景 和 纯背景
  TH1D* hP_US_SPlusB_vsR = (TH1D*)fResults->Get("sig0_bkg2/hP_US_vsR");
  TH1D* hP_US_B_vsR = (TH1D*)fResults->Get("sig0_bkg2/hP_BCD_US_BB_vsR");
  
  // SS: 信号+背景 和 纯背景
  TH1D* hP_SS_SPlusB_vsR = (TH1D*)fResults->Get("sig0_bkg2/hP_SS_vsR");
  TH1D* hP_SS_B_vsR = (TH1D*)fResults->Get("sig0_bkg2/hP_BCD_SS_BB_vsR");

  if (!hP_US_SPlusB_vsR || !hP_US_B_vsR || !hP_SS_SPlusB_vsR || !hP_SS_B_vsR) {
    std::cerr << "ERROR: Cannot find required histograms" << std::endl;
    std::cerr << "  hP_US_vsR: " << (hP_US_SPlusB_vsR ? "OK" : "MISSING") << std::endl;
    std::cerr << "  hP_US_BPlusCMinusD_vsR: " << (hP_US_B_vsR ? "OK" : "MISSING") << std::endl;
    std::cerr << "  hP_SS_vsR: " << (hP_SS_SPlusB_vsR ? "OK" : "MISSING") << std::endl;
    std::cerr << "  hP_SS_BPlusCMinusD_vsR: " << (hP_SS_B_vsR ? "OK" : "MISSING") << std::endl;
    return;
  }

  int nRbins = hFs_US->GetNbinsX();
  
  // 检查 bin 数一致性
  if (hFs_SS->GetNbinsX() != nRbins ||
      hP_US_SPlusB_vsR->GetNbinsX() != nRbins ||
      hP_US_B_vsR->GetNbinsX() != nRbins ||
      hP_SS_SPlusB_vsR->GetNbinsX() != nRbins ||
      hP_SS_B_vsR->GetNbinsX() != nRbins) {
    std::cerr << "ERROR: Number of bins mismatch!" << std::endl;
    return;
  }

  // ============================================
  // 2. 计算纯信号的自旋关联
  // ============================================
  
  // 创建输出直方图
  TH1D* hP_US_S = (TH1D*)hP_US_SPlusB_vsR->Clone("hP_US_S_vsR");
  hP_US_S->Reset();
  hP_US_S->SetTitle("Pure Signal Spin Correlation (Unlike-Sign);#DeltaR;P_{S}");
  hP_US_S->Sumw2();

  TH1D* hP_SS_S = (TH1D*)hP_SS_SPlusB_vsR->Clone("hP_SS_S_vsR");
  hP_SS_S->Reset();
  hP_SS_S->SetTitle("Pure Signal Spin Correlation (Same-Sign);#DeltaR;P_{S}");
  hP_SS_S->Sumw2();

  // 存储中间结果
  TH1D* hFs_US_copy = (TH1D*)hFs_US->Clone("hFs_US_vs_R");
  TH1D* hFs_SS_copy = (TH1D*)hFs_SS->Clone("hFs_SS_vs_R");
  TH1D* hP_US_SPlusB_copy = (TH1D*)hP_US_SPlusB_vsR->Clone("hP_US_SPlusB_vs_R");
  TH1D* hP_SS_SPlusB_copy = (TH1D*)hP_SS_SPlusB_vsR->Clone("hP_SS_SPlusB_vs_R");
  TH1D* hP_US_B_copy = (TH1D*)hP_US_B_vsR->Clone("hP_US_B_vs_R");
  TH1D* hP_SS_B_copy = (TH1D*)hP_SS_B_vsR->Clone("hP_SS_B_vs_R");

  std::cout << "\n=== Calculation of Pure Signal Spin Correlation ===" << std::endl;
  std::cout << "\n--- Unlike-Sign (Lambda+AntiLambda) ---" << std::endl;
  std::cout << "Bin\tf_R\t\tP_{S+B}\t\tP_B\t\tP_S" << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;

  for (int i = 1; i <= nRbins; i++) {
    
    // ===== US (Unlike-Sign) =====
    double P_US_SPlusB = hP_US_SPlusB_vsR->GetBinContent(i);
    double err_P_US_SPlusB = hP_US_SPlusB_vsR->GetBinError(i);
    
    double P_US_B = hP_US_B_vsR->GetBinContent(i);
    double err_P_US_B = hP_US_B_vsR->GetBinError(i);
    
    double f_US = hFs_US->GetBinContent(i);
    double err_f_US = hFs_US->GetBinError(i);

    // 检查 f_R 的有效性
    if (f_US <= 1e-6 || f_US > 1.0) {
      std::cout << "WARNING: US Bin " << i << " has invalid f_R = " << f_US << std::endl;
      hP_US_S->SetBinContent(i, 0.0);
      hP_US_S->SetBinError(i, 0.0);
    } else {
      // 公式: P_S = (P_{S+B} - (1-f_R)*P_B) / f_R
      double P_US_S = (P_US_SPlusB - (1.0 - f_US) * P_US_B) / f_US;
      
      // 误差传递
      double dP_dPSB = 1.0 / f_US;
      double dP_dPB = -(1.0 - f_US) / f_US;
      double dP_df = (P_US_B - P_US_SPlusB) / (f_US * f_US);
      
      double err_P_US_S = sqrt(
        pow(dP_dPSB * err_P_US_SPlusB, 2) +
        pow(dP_dPB * err_P_US_B, 2) +
        pow(dP_df * err_f_US, 2)
      );
      
      hP_US_S->SetBinContent(i, P_US_S);
      hP_US_S->SetBinError(i, err_P_US_S);
      
      printf("%d\t%.4f ± %.4f\t%.4f ± %.4f\t%.4f ± %.4f\t%.4f ± %.4f\n",
             i, f_US, err_f_US, P_US_SPlusB, err_P_US_SPlusB, 
             P_US_B, err_P_US_B, P_US_S, err_P_US_S);
    }

    // ===== SS (Same-Sign) =====
    double P_SS_SPlusB = hP_SS_SPlusB_vsR->GetBinContent(i);
    double err_P_SS_SPlusB = hP_SS_SPlusB_vsR->GetBinError(i);
    
    double P_SS_B = hP_SS_B_vsR->GetBinContent(i);
    double err_P_SS_B = hP_SS_B_vsR->GetBinError(i);
    
    double f_SS = hFs_SS->GetBinContent(i);
    double err_f_SS = hFs_SS->GetBinError(i);

    if (f_SS <= 1e-6 || f_SS > 1.0) {
      std::cout << "WARNING: SS Bin " << i << " has invalid f_R = " << f_SS << std::endl;
      hP_SS_S->SetBinContent(i, 0.0);
      hP_SS_S->SetBinError(i, 0.0);
    } else {
      double P_SS_S = (P_SS_SPlusB - (1.0 - f_SS) * P_SS_B) / f_SS;
      
      double dP_dPSB = 1.0 / f_SS;
      double dP_dPB = -(1.0 - f_SS) / f_SS;
      double dP_df = (P_SS_B - P_SS_SPlusB) / (f_SS * f_SS);
      
      double err_P_SS_S = sqrt(
        pow(dP_dPSB * err_P_SS_SPlusB, 2) +
        pow(dP_dPB * err_P_SS_B, 2) +
        pow(dP_df * err_f_SS, 2)
      );
      
      hP_SS_S->SetBinContent(i, P_SS_S);
      hP_SS_S->SetBinError(i, err_P_SS_S);
    }
  }

  std::cout << "\n--- Same-Sign (LambdaLambda + AntiLambdaAntiLambda) ---" << std::endl;
  std::cout << "Bin\tf_R\t\tP_{S+B}\t\tP_B\t\tP_S" << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;
  
  for (int i = 1; i <= nRbins; i++) {
    double f_SS = hFs_SS->GetBinContent(i);
    double err_f_SS = hFs_SS->GetBinError(i);
    double P_SS_SPlusB = hP_SS_SPlusB_vsR->GetBinContent(i);
    double err_P_SS_SPlusB = hP_SS_SPlusB_vsR->GetBinError(i);
    double P_SS_B = hP_SS_B_vsR->GetBinContent(i);
    double err_P_SS_B = hP_SS_B_vsR->GetBinError(i);
    double P_SS_S = hP_SS_S->GetBinContent(i);
    double err_P_SS_S = hP_SS_S->GetBinError(i);
    
    printf("%d\t%.4f ± %.4f\t%.4f ± %.4f\t%.4f ± %.4f\t%.4f ± %.4f\n",
           i, f_SS, err_f_SS, P_SS_SPlusB, err_P_SS_SPlusB, 
           P_SS_B, err_P_SS_B, P_SS_S, err_P_SS_S);
  }
  
  std::cout << "========================================================" << std::endl;

  // ============================================
  // 3. 保存直方图
  // ============================================
  outFile->cd();
  
  // 纯信号结果
  hP_US_S->Write();
  hP_SS_S->Write();
  
  // 信号纯度
  hFs_US_copy->Write();
  hFs_SS_copy->Write();
  
  // 信号+背景
  hP_US_SPlusB_copy->Write();
  hP_SS_SPlusB_copy->Write();
  
  // 纯背景
  hP_US_B_copy->Write();
  hP_SS_B_copy->Write();

  // ============================================
  // 4. 绘制对比图
  // ============================================
  
  // 图1: US 对比
  TCanvas* cUS = new TCanvas("cUS", "US Spin Correlation Comparison", 1200, 800);
  cUS->SetLeftMargin(0.13);
  cUS->SetBottomMargin(0.12);
  
  hP_US_S->SetMarkerStyle(20);
  hP_US_S->SetMarkerColor(kRed);
  hP_US_S->SetLineColor(kRed);
  hP_US_S->GetYaxis()->SetRangeUser(-0.05, 0.05);
  hP_US_S->GetXaxis()->SetTitle("#DeltaR");
  hP_US_S->GetYaxis()->SetTitle("P");
  hP_US_S->SetTitle("Spin Correlation Comparison (Unlike-Sign: #Lambda#bar{#Lambda}); #DeltaR; P");
  
  hP_US_S->Draw("E1");
  hP_US_SPlusB_copy->SetMarkerStyle(21);
  hP_US_SPlusB_copy->SetMarkerColor(kBlue);
  hP_US_SPlusB_copy->SetLineColor(kBlue);
  hP_US_SPlusB_copy->Draw("E1 SAME");
  hP_US_B_copy->SetMarkerStyle(22);
  hP_US_B_copy->SetMarkerColor(kGreen+2);
  hP_US_B_copy->SetLineColor(kGreen+2);
  hP_US_B_copy->Draw("E1 SAME");
  
  TLegend* legUS = new TLegend(0.65, 0.70, 0.88, 0.88);
  legUS->SetBorderSize(0);
  legUS->SetFillStyle(0);
  legUS->AddEntry(hP_US_S, "Pure Signal P_{S}", "lp");
  legUS->AddEntry(hP_US_SPlusB_copy, "Signal+Background P_{S+B}", "lp");
  legUS->AddEntry(hP_US_B_copy, "Pure Background P_{B}", "lp");
  legUS->Draw();
  
  cUS->Write();
  cUS->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AddsignalFrac/SpinCorrelation_Comparison_US.png");
  delete cUS;

  // 图2: SS 对比
  TCanvas* cSS = new TCanvas("cSS", "SS Spin Correlation Comparison", 1200, 800);
  cSS->SetLeftMargin(0.13);
  cSS->SetBottomMargin(0.12);
  
  hP_SS_S->SetMarkerStyle(20);
  hP_SS_S->SetMarkerColor(kRed);
  hP_SS_S->SetLineColor(kRed);
  hP_SS_S->GetYaxis()->SetRangeUser(-0.05, 0.05);
  hP_SS_S->GetXaxis()->SetTitle("#DeltaR");
  hP_SS_S->GetYaxis()->SetTitle("P");
  hP_SS_S->SetTitle("Spin Correlation Comparison (Same-Sign: #Lambda#Lambda + #bar{#Lambda}#bar{#Lambda}); #DeltaR; P");
 
  
  hP_SS_S->Draw("E1");
  hP_SS_SPlusB_copy->SetMarkerStyle(21);
  hP_SS_SPlusB_copy->SetMarkerColor(kBlue);
  hP_SS_SPlusB_copy->SetLineColor(kBlue);
  hP_SS_SPlusB_copy->Draw("E1 SAME");
  hP_SS_B_copy->SetMarkerStyle(22);
  hP_SS_B_copy->SetMarkerColor(kGreen+2);
  hP_SS_B_copy->SetLineColor(kGreen+2);
  hP_SS_B_copy->Draw("E1 SAME");
  
  TLegend* legSS = new TLegend(0.65, 0.70, 0.88, 0.88);
  legSS->SetBorderSize(0);
  legSS->SetFillStyle(0);
  legSS->AddEntry(hP_SS_S, "Pure Signal P_{S}", "lp");
  legSS->AddEntry(hP_SS_SPlusB_copy, "Signal+Background P_{S+B}", "lp");
  legSS->AddEntry(hP_SS_B_copy, "Pure Background P_{B}", "lp");
  legSS->Draw();
  
  cSS->Write();
  cSS->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AddsignalFrac/SpinCorrelation_Comparison_SS.png");
  delete cSS;

  // 图3: US vs SS 纯信号对比
  TCanvas* cCompare = new TCanvas("cCompare", "Pure Signal Comparison", 1200, 800);
  cCompare->SetLeftMargin(0.13);
  cCompare->SetBottomMargin(0.12);
  
  hP_US_S->SetMarkerStyle(20);
  hP_US_S->SetMarkerColor(kBlue);
  hP_US_S->SetLineColor(kBlue);
  hP_US_S->GetYaxis()->SetRangeUser(-0.05, 0.05);
  hP_US_S->GetXaxis()->SetTitle("#DeltaR");
  hP_US_S->GetYaxis()->SetTitle("P_{S}");
  hP_US_S->SetTitle("Pure Signal Spin Correlation Comparison; #DeltaR; P_{S}");
  
  hP_US_S->Draw("E1");
  hP_SS_S->SetMarkerStyle(21);
  hP_SS_S->SetMarkerColor(kRed);
  hP_SS_S->SetLineColor(kRed);
  hP_SS_S->Draw("E1 SAME");
  
  TLegend* legComp = new TLegend(0.65, 0.70, 0.88, 0.88);
  legComp->SetBorderSize(0);
  legComp->SetFillStyle(0);
  legComp->AddEntry(hP_US_S, "Unlike-Sign (#Lambda#bar{#Lambda})", "lp");
  legComp->AddEntry(hP_SS_S, "Same-Sign (#Lambda#Lambda + #bar{#Lambda}#bar{#Lambda})", "lp");
  legComp->Draw();
  
  cCompare->Write();
  cCompare->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AddsignalFrac/PureSignal_Comparison_US_vs_SS.png");
  delete cCompare;


  TCanvas* sigFrac = new TCanvas("sigFrac", "Pure Signal Comparison", 1200, 600);
  sigFrac->Divide(2,1);
  sigFrac->cd(1);
  gPad->SetLeftMargin(0.15);
  hFs_SS->Draw();

  
  
  sigFrac->cd(2);
  gPad->SetLeftMargin(0.15);
  hFs_US->Draw();
  
  sigFrac->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AddsignalFrac/frac_US_SS.png");
  delete sigFrac;


  // ============================================
  // 5. 清理
  // ============================================
  outFile->Close();
  
  fUS->Close();
  fSS->Close();
  fResults->Close();
  
  delete fUS;
  delete fSS;
  delete fResults;
  delete outFile;

  std::cout << "\nDone! Output saved to: Spincore_withFrac.root" << std::endl;
  std::cout << "  - hP_US_S_vsR: Pure signal P for Unlike-Sign" << std::endl;
  std::cout << "  - hP_SS_S_vsR: Pure signal P for Same-Sign" << std::endl;
  std::cout << "Comparison plots saved to:" << std::endl;
  std::cout << "  - SpinCorrelation_Comparison_US.png" << std::endl;
  std::cout << "  - SpinCorrelation_Comparison_SS.png" << std::endl;
  std::cout << "  - PureSignal_Comparison_US_vs_SS.png" << std::endl;
}


//draw the data and simulation
//add RMS to the uncertainty of the default method
//add RMS to the uncertainty of the default method

void drawPloarizationcorrPlots_And_Simulation(){
  
  TFile* outFile = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "UPDATE");


  TCanvas* c1 = new TCanvas("c1", "c1", 2400, 1200);
  gPad->SetLeftMargin(0.15);
  TFile* f_US = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/US_Uncertainty/US_All_Differences.root", "READ");
  TFile* f_SS = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/SS_Uncertainty/SS_All_Differences.root", "READ");

  TH1D* h_US_YieldSS_sig0_bkg2 = (TH1D*)f_US->Get("US_do_YieldSS_BCD_sig0_bkg2");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetRangeUser(-0.1, 0.1);

  h_US_YieldSS_sig0_bkg2->SetLineColor(kBlue);
  h_US_YieldSS_sig0_bkg2->SetLineWidth(2);
  h_US_YieldSS_sig0_bkg2->SetTitle("Spin correlation P_{#Lambda_{1},#Lambda_{2}} as a function of pair separation #DeltaR");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitle("P_{#Lambda_{1},#Lambda_{2}}");
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitleSize(0.05);
  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}");
  
  TH1D* h_SS_YieldSS_sig0_bkg2 = (TH1D*)f_SS->Get("SS_do_YieldSS_BCD_sig0_bkg2");
  h_SS_YieldSS_sig0_bkg2->SetLineColor(kRed);
  h_SS_YieldSS_sig0_bkg2->SetLineWidth(2);
  
  h_US_YieldSS_sig0_bkg2->Draw("E1");
  h_US_YieldSS_sig0_bkg2->SetMarkerSize(3);
  h_US_YieldSS_sig0_bkg2->SetMarkerColor(kBlue);
  h_US_YieldSS_sig0_bkg2->GetXaxis()->CenterTitle();
  h_US_YieldSS_sig0_bkg2->GetYaxis()->CenterTitle();
  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetTitleOffset(1.30);
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetTitleOffset(1.0);
  h_US_YieldSS_sig0_bkg2->GetXaxis()->SetLabelSize(0.05);
  h_US_YieldSS_sig0_bkg2->GetYaxis()->SetLabelSize(0.05);

  h_SS_YieldSS_sig0_bkg2->SetMarkerSize(3);
  h_SS_YieldSS_sig0_bkg2->SetMarkerColor(kRed);
  h_SS_YieldSS_sig0_bkg2->Draw("E1 SAME");  
 
  TFile* RMShistFile = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/All_Plots_Diff.root", "READ");
  TH1D* h_RMS = (TH1D*)RMShistFile->Get("h_RMS_US");
  TH1D* h_RMS_SS = (TH1D*)RMShistFile->Get("h_RMS_SS");

  double box_width = 0.08;  
  double box_alpha = 0.3;   
  
  for(int bin = 1; bin <= h_US_YieldSS_sig0_bkg2->GetNbinsX(); ++bin){
    double bin_center = h_US_YieldSS_sig0_bkg2->GetBinCenter(bin);
    double bin_left = bin_center - box_width/2;
    double bin_right = bin_center + box_width/2;
    double pol_value = h_US_YieldSS_sig0_bkg2->GetBinContent(bin);
    double pol_value_SS = h_SS_YieldSS_sig0_bkg2->GetBinContent(bin);
    double sys_error = h_RMS->GetBinContent(bin);
    double sys_error_SS = h_RMS_SS->GetBinContent(bin);
    
    if (sys_error > 0) {
      double box_bottom = pol_value - sys_error;
      double box_top = pol_value + sys_error;
      TBox *box = new TBox(bin_left, box_bottom, bin_right, box_top);
      box->SetFillColor(kGray);
      box->SetFillStyle(0);  
      box->SetLineColor(kBlue);
      box->SetLineWidth(2);
      box->SetFillColorAlpha(kGray, box_alpha);  
      box->Draw("SAME");

      double box_bottom_SS = pol_value_SS - sys_error_SS;
      double box_top_SS = pol_value_SS + sys_error_SS;
      TBox *box_SS = new TBox(bin_left, box_bottom_SS, bin_right, box_top_SS);
      box_SS->SetFillColor(kGray);
      box_SS->SetFillStyle(0);  
      box_SS->SetLineColor(kRed);
      box_SS->SetLineWidth(2);
      box_SS->SetFillColorAlpha(kGray, box_alpha);  
      box_SS->Draw("SAME");
    }
  }
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->DrawLatex(0.18, 0.40, "ALICE");
  tex->DrawLatex(0.18, 0.35, "p+p #sqrt{s} = 13.6 TeV");
  tex->DrawLatex(0.18, 0.30, "|y| < 0.5");
  tex->DrawLatex(0.18, 0.25, "0.8 GeV/c < p_{T} < 3 GeV/c");
  tex->DrawLatex(0.18, 0.20, "Helicity frame");
  tex->DrawLatex(0.18, 0.15, "Signal extraction:2D side band method");
  gPad->SetGridx();  
  gPad->SetGridy(); 

  TLine *line = new TLine(
    h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmin(),
    0,
    h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmax(),
    0
  );
  line->SetLineStyle(10);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1);
  line->Draw("same");

  // ==================== 读取模拟数据 ====================
  TFile *f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Pythia_code/merged_output_CR1_v1.root");
 
    TProfile *prof_R_SU6_nominal   = (TProfile*)f->Get("prof_R_SU6_nominal");
    TProfile *prof_R_BJ_nominal    = (TProfile*)f->Get("prof_R_BJ_nominal");
    TProfile *prof_R_SU6_pert      = (TProfile*)f->Get("prof_R_SU6_pert");
    TProfile *prof_R_BJ_pert       = (TProfile*)f->Get("prof_R_BJ_pert");
    TProfile *prof_R_SU6_nonpert   = (TProfile*)f->Get("prof_R_SU6_nonpert");
    TProfile *prof_R_BJ_nonpert    = (TProfile*)f->Get("prof_R_BJ_nonpert");

    prof_R_SU6_nonpert->SetStats(0);

    // Non-pert SU6
    prof_R_SU6_nonpert->SetLineColor(kRed);
    prof_R_SU6_nonpert->SetLineWidth(3);
    prof_R_SU6_nonpert->SetMarkerColor(kRed);
    prof_R_SU6_nonpert->SetMarkerStyle(20);
    prof_R_SU6_nonpert->SetMarkerSize(0);
    prof_R_SU6_nonpert->SetFillColor(kRed);
    prof_R_SU6_nonpert->Draw("E3 SAME");
    prof_R_SU6_nonpert->GetYaxis()->SetRangeUser(-0.02, 0.1);

    // Non-pert BJ
    prof_R_BJ_nonpert->SetLineColor(kBlue);
    prof_R_BJ_nonpert->SetLineWidth(3);
    prof_R_BJ_nonpert->SetMarkerColor(kBlue);
    prof_R_BJ_nonpert->SetMarkerStyle(21);
    prof_R_BJ_nonpert->SetMarkerSize(0);
    prof_R_BJ_nonpert->SetFillColor(kBlue);
    prof_R_BJ_nonpert->Draw("SAME E3");

    // Pert SU6
    prof_R_SU6_pert->SetLineColor(kGreen);
    prof_R_SU6_pert->SetLineWidth(3);
    prof_R_SU6_pert->SetMarkerColor(kGreen);
    prof_R_SU6_pert->SetMarkerStyle(22);
    prof_R_SU6_pert->SetMarkerSize(0);
    prof_R_SU6_pert->SetFillColor(kGreen);
    prof_R_SU6_pert->Draw("SAME E3");

    // Pert BJ
    prof_R_BJ_pert->SetLineColor(kMagenta);
    prof_R_BJ_pert->SetLineWidth(3);
    prof_R_BJ_pert->SetMarkerColor(kMagenta);
    prof_R_BJ_pert->SetMarkerStyle(23);
    prof_R_BJ_pert->SetMarkerSize(0);
    prof_R_BJ_pert->SetFillColor(kMagenta);
    prof_R_BJ_pert->Draw("SAME E3");

    // Total SU6
    prof_R_SU6_nominal->SetLineColor(kBlack);
    prof_R_SU6_nominal->SetLineWidth(3);
    prof_R_SU6_nominal->SetMarkerColor(kBlack);
    prof_R_SU6_nominal->SetMarkerStyle(33);
    prof_R_SU6_nominal->SetMarkerSize(0);
    prof_R_SU6_nominal->SetFillColor(kBlack);
    prof_R_SU6_nominal->Draw("SAME E3");

    // Total BJ
    prof_R_BJ_nominal->SetLineColor(kOrange);
    prof_R_BJ_nominal->SetLineWidth(3);
    prof_R_BJ_nominal->SetMarkerColor(kOrange);
    prof_R_BJ_nominal->SetMarkerStyle(34);
    prof_R_BJ_nominal->SetMarkerSize(0);
    prof_R_BJ_nominal->SetFillColor(kOrange);
    prof_R_BJ_nominal->Draw("SAME E3");

    // SU6 理论预测条带 (使用不同的变量名)
    double val_SU6 = 0.096;
    double err_SU6 = 0.004;
    double ymin_SU6 = val_SU6 - err_SU6;
    double ymax_SU6 = val_SU6 + err_SU6;
    double xmin_sim = prof_R_SU6_nonpert->GetXaxis()->GetXmin();
    double xmax_sim = prof_R_SU6_nonpert->GetXaxis()->GetXmax();

    TBox *box_SU6 = new TBox(xmin_sim, ymin_SU6, xmax_sim, ymax_SU6);
    box_SU6->SetFillColor(kGray);
    box_SU6->SetFillStyle(3001);
    box_SU6->Draw();

    TLine *line_SU6 = new TLine(xmin_sim, val_SU6, xmax_sim, val_SU6);
    line_SU6->SetLineColor(kBlack);
    line_SU6->SetLineWidth(2);
    line_SU6->SetLineStyle(2);
    line_SU6->Draw();

    // BJ 理论预测条带 (使用不同的变量名)
    double val_BJ = 0.015;
    double err_BJ = 0.002;
    double ymin_BJ = val_BJ - err_BJ;
    double ymax_BJ = val_BJ + err_BJ;

    TBox *box_BJ = new TBox(xmin_sim, ymin_BJ, xmax_sim, ymax_BJ);
    box_BJ->SetFillColor(kRed);
    box_BJ->SetFillStyle(3001);
    box_BJ->Draw();

    TLine *line_BJ = new TLine(xmin_sim, val_BJ, xmax_sim, val_BJ);
    line_BJ->SetLineColor(kBlack);
    line_BJ->SetLineWidth(2);
    line_BJ->SetLineStyle(2);
    line_BJ->Draw();

    // 图例
    TLegend *leg_sim = new TLegend(0.7, 0.1, 0.95, 0.35);
    leg_sim->AddEntry(prof_R_SU6_nominal, "SU(6) nominal", "lf");
    leg_sim->AddEntry(prof_R_BJ_nominal, "BJ nominal", "lf");
    leg_sim->AddEntry(prof_R_SU6_nonpert, "SU(6) non-pert", "lf");
    leg_sim->AddEntry(prof_R_BJ_nonpert, "BJ non-pert", "lf");
    leg_sim->AddEntry(prof_R_SU6_pert, "SU(6) pert", "lf");
    leg_sim->AddEntry(prof_R_BJ_pert, "BJ pert", "lf");
    leg_sim->SetBorderSize(0);
    leg_sim->SetFillStyle(0);
    leg_sim->Draw();

   

    TLatex latex_sim;
    latex_sim.SetNDC();
    latex_sim.SetTextSize(0.03);
    latex_sim.DrawLatex(0.48, 0.45, "pp #sqrt{s} = 13.6 TeV, PYTHIA8 Monash, CR model 1");
    latex_sim.DrawLatex(0.48, 0.40, "#Lambda#bar{#Lambda} pairs:|y| < 0.5, 0.8 < p_{T} < 3.0 GeV/c");
  

  // ==================== 第一个理论预测条带 (BJ) ====================
  double val_BJ1 = 0.015;
  double err_BJ1 = 0.002;
  double ymin_BJ1 = val_BJ1 - err_BJ1;
  double ymax_BJ1 = val_BJ1 + err_BJ1;
  double xmin_data = h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmin();
  double xmax_data = h_US_YieldSS_sig0_bkg2->GetXaxis()->GetXmax();

  TBox *box_BJ_data = new TBox(xmin_data, ymin_BJ1, xmax_data, ymax_BJ1);
  box_BJ_data->SetFillColor(kRed);
  box_BJ_data->SetFillStyle(3001);
  box_BJ_data->Draw();

  TLine *line_BJ_data = new TLine(xmin_data, val_BJ1, xmax_data, val_BJ1);
  line_BJ_data->SetLineColor(kBlack);
  line_BJ_data->SetLineWidth(2);
  line_BJ_data->SetLineStyle(2);
  line_BJ_data->Draw("same");

  // 主图例
  

  TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.9);
  legend->AddEntry(h_US_YieldSS_sig0_bkg2, "Unlike-Sign : #Lambda#bar{#Lambda}", "lp");
  legend->AddEntry(h_SS_YieldSS_sig0_bkg2, "Like-Sign : #Lambda#Lambda", "lp");
  legend->AddEntry(box_BJ_data, "STAR simulation: BJ: 0.015 #pm 0.002", "f");  
  legend->AddEntry(box_SU6, "STAR simulation: SU(6): 0.096 #pm 0.004", "f");  


  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();

  // ==================== c2: LL vs ALAL ====================
  TCanvas* c2 = new TCanvas("c2", "c2", 1800, 1200);
  TFile* f1 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root", "READ");

  TH1D* h_LL = (TH1D*)f1->Get("sig0_bkg2/hP_LL_vsR;1");
  TH1D* h_ALAL = (TH1D*)f1->Get("sig0_bkg2/hP_ALAL_vsR;1");
  h_LL->GetYaxis()->SetRangeUser(-0.02, 0.05);
  h_LL->SetLineColor(kRed);
  h_LL->GetXaxis()->SetTitle("#DeltaR = #sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}}");
  h_LL->SetTitle("Spin correlation P_{#Lambda_{1},#Lambda_{2}} as a function of pair separation #DeltaR");
  h_LL->GetYaxis()->SetTitle("P_{#Lambda_{1},#Lambda_{2}}");
  
  h_LL->Draw("E1");
  h_LL->SetMarkerSize(3);
  h_LL->SetMarkerColor(kRed);
  h_LL->GetXaxis()->CenterTitle();
  h_LL->GetYaxis()->CenterTitle();
  h_LL->GetXaxis()->SetTitleOffset(1.30);
  h_LL->GetYaxis()->SetTitleOffset(1.0);
  h_LL->SetLineWidth(2);
  h_LL->GetYaxis()->SetTitleSize(0.05);
  h_LL->Draw();

  h_ALAL->SetLineColor(kBlue);
  h_ALAL->Draw("E1 SAME");
  h_ALAL->SetMarkerSize(3);
  h_ALAL->SetMarkerColor(kBlue);
  h_ALAL->SetLineWidth(2);
  h_ALAL->Draw("SAME");

  TLegend* legend1 = new TLegend(0.6, 0.7, 0.88, 0.9);
  legend1->AddEntry(h_LL, "#Lambda#Lambda", "lp");
  legend1->AddEntry(h_ALAL, "#bar{#Lambda}#bar{#Lambda}", "lp");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->Draw();
  gPad->SetGridx();
  gPad->SetGridy();

  TLatex *tex1 = new TLatex();
  tex1->SetNDC();
  tex1->SetTextSize(0.03);
  tex1->DrawLatex(0.18, 0.85, "ALICE");
  tex1->DrawLatex(0.18, 0.80, "p+p #sqrt{s} = 13.6 TeV");
  tex1->DrawLatex(0.18, 0.75, "|y| < 0.5");
  tex1->DrawLatex(0.18, 0.70, "0.8 GeV/c < p_{T} < 3 GeV/c");
  tex1->DrawLatex(0.18, 0.65, "Helicity frame");
  tex1->DrawLatex(0.18, 0.60, "Signal extraction:2D side band method");

  TLine *line_zero = new TLine(
    h_LL->GetXaxis()->GetXmin(),
    0,
    h_LL->GetXaxis()->GetXmax(),
    0
  );
  line_zero->SetLineStyle(10);
  line_zero->SetLineColor(kBlack);
  line_zero->SetLineWidth(5);
  line_zero->Draw("same");

  // ==================== c3: BB ====================
  TCanvas* c3 = new TCanvas("c3", "c3", 1800, 1200);
  TFile* f2 = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Results_do_YieldSS_BCD_Rdep_Output_pp_sig0_bkg2.root", "READ");
  TCanvas* h_BB = (TCanvas*)f2->Get("sig0_bkg2/cP_vsR_SB;1");
  h_BB->Draw();

  // ==================== 保存 ====================
  outFile->cd();
  c1->Write("c_P_with_sys");
  c2->Write("h_LL");
  
  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/P_with_sys_sim.png");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_LL.png");
  c3->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_BB.png");

  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/P_with_sys_sim.pdf");
  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_LL.pdf");
  c3->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/h_BB.pdf");
  
  outFile->Save();
  
  delete c1;
  delete c2;
  delete c3;
  
  f_US->Close();
  f_SS->Close();
  f1->Close();
  f2->Close();
  if (f) f->Close();
  RMShistFile->Close();
  outFile->Close();
}

void readTheEffSigma(){
  const double REdges[] = {0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1};
  for(int i = 0;i<=5;i++){
    TCanvas* c1 = new TCanvas("c1", "c1", 2000, 800);
    c1->Divide(2,1);
    // US (Unlike-Sign: Lambda + AntiLambda) 的 fs
    TFile* fUS = TFile::Open(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Output_pp_sig0_bkg2_hSparseLambdaAntiLambda/fit_Rbin%d.root",i));
  
    TH1D* hSigma_X = (TH1D*)fUS->Get("hSigmaEffX_vsCos");
    
    TH1D* hSigma_Y = (TH1D*)fUS->Get("hSigmaEffY_vsCos");
    c1->cd(1);
    gPad->SetLeftMargin(0.15); 
    hSigma_X->Draw();
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.16, 0.86, Form("R bin:[%2.2f, %.2f]", REdges[i], REdges[i+1]));
    c1->cd(2);
    gPad->SetLeftMargin(0.15); 
    hSigma_Y->Draw("same");
  
    c1->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/EffSigmafigure/Eff_Sigma_Rbin%d.png",i));
    delete c1;

    TCanvas* c2= new TCanvas("c1", "c1", 2000, 800);
    c2->Divide(2,1);
    c2->cd(1);
    gPad->SetLeftMargin(0.15); 
    TH1D* hWinXlow_SE = (TH1D*)fUS->Get("hWinXlow_SE");
    TH1D* hWinXhi_SE = (TH1D*)fUS->Get("hWinXhi_SE");
    hWinXlow_SE->GetYaxis()->SetRangeUser(1.1,1.13);
    hWinXlow_SE->GetYaxis()->SetTitle("GeV/c^{2}");
    hWinXlow_SE->Draw();
    hWinXhi_SE->Draw("same");
    TLatex latx;
    latx.SetNDC(true);
    latx.SetTextSize(0.040);
    latx.DrawLatex(0.16, 0.86, Form("R bin:[%2.2f, %.2f]", REdges[i], REdges[i+1]));
   

    c2->cd(2);
    gPad->SetLeftMargin(0.15); 
    TH1D* hWinXlow_ME = (TH1D*)fUS->Get("hWinXlow_used_ME");
    TH1D* hWinXhi_ME = (TH1D*)fUS->Get("hWinXhi_used_ME");
    hWinXlow_ME->GetYaxis()->SetRangeUser(1.1,1.13);
    hWinXlow_ME->GetYaxis()->SetTitle("GeV/c^{2}");
    hWinXlow_ME->Draw();
    hWinXhi_ME->Draw("same");

    c2->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/EffSigmafigure/SignalWindowX_Rbin%d.png",i));

    delete c2;

    TCanvas* c3= new TCanvas("c1", "c1", 2000, 800);
    c3->Divide(2,1);
    c3->cd(1);
   
    gPad->SetLeftMargin(0.15); 
    TH1D* hWinYlow_SE = (TH1D*)fUS->Get("hWinYlow_SE");
    TH1D* hWinYhi_SE = (TH1D*)fUS->Get("hWinYhi_SE");
    hWinYlow_SE->GetYaxis()->SetRangeUser(1.1,1.13);
    hWinYlow_SE->GetYaxis()->SetTitle("GeV/c^{2}");
   
    hWinYlow_SE->Draw();
     TLatex lat1;
    lat1.SetNDC(true);
    lat1.SetTextSize(0.040);
    lat1.DrawLatex(0.16, 0.86, Form("R bin:[%2.2f, %.2f]", REdges[i], REdges[i+1]));
    hWinYhi_SE->Draw("same");
   
   

    c3->cd(2);
    gPad->SetLeftMargin(0.15); 
    TH1D* hWinYlow_ME = (TH1D*)fUS->Get("hWinYlow_used_ME");
    TH1D* hWinYhi_ME = (TH1D*)fUS->Get("hWinYhi_used_ME");
    hWinYlow_ME->GetYaxis()->SetRangeUser(1.1,1.13);
    hWinYlow_ME->GetYaxis()->SetTitle("GeV/c^{2}");
    hWinYlow_ME->Draw();
    hWinYhi_ME->Draw("same");
    

    c3->SaveAs(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/EffSigmafigure/SignalWindowY_Rbin%d.png",i));
    delete c3;


  }
 
}

void drawRdistribution(){
  // 1. 打开文件
  TCanvas* c1 = new TCanvas("c1", "c1", 2000, 1400);
  c1->Divide(2,2);
  ///home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/MC_file/MC_data/AnalysisResults_MC_unweight_691539.root
  //home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_setting5.root
  TFile *f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_675285_MEV3.root");
  if (!f){cout<<"无法打开文件！"<<endl; return;} 

  // ========== 设置 R 范围 ==========
  double R_low = 0.0; // 根据需要修改
  double R_high = 3.1;  // 根据需要修改

  // ========== Unlike-sign (Lambda + AntiLambda) ==========
  THnSparseF *hnSE = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaAntiLambda;1");
  THnSparseF *hnME = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaAntiLambdaMixed;1");
  
  if (!hnSE || !hnME) { cout << "找不到稀疏直方图" << endl; return; }

  // 设置 R 轴范围（索引3）
  hnSE->GetAxis(3)->SetRange(hnSE->GetAxis(3)->FindBin(R_low + 1e-6), 
                             hnSE->GetAxis(3)->FindBin(R_high - 1e-6));
  hnME->GetAxis(3)->SetRange(hnME->GetAxis(3)->FindBin(R_low + 1e-6), 
                             hnME->GetAxis(3)->FindBin(R_high - 1e-6));

  // 投影 cosθ（索引2）
  TH1D *hRSE = hnSE->Projection(3); 
  TH1D *hRME = hnME->Projection(3); 

  hRSE->SetLineColor(kRed);
  hRSE->SetLineWidth(4);
  hRSE->Scale(1.0/hRSE->Integral());

  hRME->SetLineColor(kBlue);
  hRME->SetLineWidth(4);
  hRME->Scale(1.0/hRME->Integral());
  
  c1->cd(1);
  hRSE->SetTitle(Form("Unlike-sign R distribution (R in [%.1f, %.1f]);R;1/N dN/dR", R_low, R_high));
  hRSE->SetStats(0);
  hRSE->Draw();
  hRME->Draw("same");
  
  TLegend *leg1 = new TLegend(0.7, 0.7, 0.88, 0.88);
  leg1->AddEntry(hRSE, "SE", "l");
  leg1->AddEntry(hRME, "ME", "l");
  leg1->Draw();

  c1->cd(2);
  TH1D *hRatio = (TH1D*)hRSE->Clone("hRatio");
  
  hRatio->Divide(hRME);
  hRatio->GetYaxis()->SetRangeUser(0.98, 1.02);
  hRatio->SetTitle(Form("SE/ME ratio (R in [%.1f, %.1f]);R;SE/ME", R_low, R_high));
  hRatio->Draw();

  // ========== Like-sign (LL + ALAL) ==========
  c1->cd(3);
  
  // LambdaLambda
  THnSparseF *hnSE_LL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaLambda;1");
  THnSparseF *hnME_LL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseLambdaLambdaMixed;1");
  
  // AntiLambdaAntiLambda
  THnSparseF *hnSE_ALAL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseAntiLambdaAntiLambda;1");
  THnSparseF *hnME_ALAL = (THnSparseF*)f->Get("lambdaspincorrderived/hSparseAntiLambdaAntiLambdaMixed;1");

  // 设置 R 轴范围
  hnSE_LL->GetAxis(3)->SetRange(hnSE_LL->GetAxis(3)->FindBin(R_low + 1e-6), 
                                hnSE_LL->GetAxis(3)->FindBin(R_high - 1e-6));
  hnME_LL->GetAxis(3)->SetRange(hnME_LL->GetAxis(3)->FindBin(R_low + 1e-6), 
                                hnME_LL->GetAxis(3)->FindBin(R_high - 1e-6));
  hnSE_ALAL->GetAxis(3)->SetRange(hnSE_ALAL->GetAxis(3)->FindBin(R_low + 1e-6), 
                                  hnSE_ALAL->GetAxis(3)->FindBin(R_high - 1e-6));
  hnME_ALAL->GetAxis(3)->SetRange(hnME_ALAL->GetAxis(3)->FindBin(R_low + 1e-6), 
                                  hnME_ALAL->GetAxis(3)->FindBin(R_high - 1e-6));

  // 投影 cosθ
  TH1D *hRSE_LL = hnSE_LL->Projection(3); 
  TH1D *hRME_LL = hnME_LL->Projection(3); 
  TH1D *hRSE_ALAL = hnSE_ALAL->Projection(3); 
  TH1D *hRME_ALAL = hnME_ALAL->Projection(3); 

  // 合并 LL + ALAL
  hRSE_LL->Add(hRSE_ALAL);
  hRME_LL->Add(hRME_ALAL);

  hRSE_LL->SetLineColor(kRed);
  hRSE_LL->SetLineWidth(4);
  hRSE_LL->Scale(1.0/hRSE_LL->Integral());

  hRME_LL->SetLineColor(kBlue);
  hRME_LL->SetLineWidth(4);
  hRME_LL->Scale(1.0/hRME_LL->Integral());

  hRSE_LL->SetTitle(Form("Like-sign R distribution (R in [%.1f, %.1f]);R;1/N dN/dR", R_low, R_high));
  hRSE_LL->SetStats(0);
  hRSE_LL->Draw();
  hRME_LL->Draw("same");
  
  TLegend *leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
  leg2->AddEntry(hRSE_LL, "SE", "l");
  leg2->AddEntry(hRME_LL, "ME", "l");
  leg2->Draw();

  c1->cd(4);
  TH1D *hRatioLL = (TH1D*)hRSE_LL->Clone("hRatioLL");
  hRatioLL->Divide(hRME_LL);
  hRatioLL->GetYaxis()->SetRangeUser(0.98, 1.02);
  hRatioLL->SetTitle(Form("SE/ME ratio (R in [%.1f, %.1f]);R;SE/ME", R_low, R_high));
  hRatioLL->Draw();

  // 保存
  c1->SaveAs(Form("Picture/cosTheta_distribution_R_%.1f_%.1f.png", R_low, R_high));
  delete c1;
}

void CheckPtEtaYSEAndMEdistribution(){
  TCanvas* c1 = new  TCanvas("c1","c1", 1800, 1400);
  c1->Divide(2,2);
  //home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_setting5.root
 
   TFile *f = TFile::Open("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_675285_MEV3.root");
  if (!f){cout<<"无法打开文件！"<<endl; return;} 

 
  TH2D *hPtYSE = (TH2D*)f->Get("lambdaspincorrderived/hPtYSame;1");
  TH2D *hPtYME = (TH2D*)f->Get("lambdaspincorrderived/hPtYMix;1");

  TH1D *hPtSE = hPtYSE->ProjectionX("hPtSE");
  TH1D *hPtME = hPtYME->ProjectionX("hPtME");
 
 
 
  hPtSE->Scale(1.0/hPtSE->Integral());
  hPtME->Scale(1.0/hPtME->Integral());
  hPtSE->SetLineColor(kRed);
  hPtME->SetLineColor(kBlue);
  hPtSE->SetLineWidth(4);
  hPtME->SetLineWidth(4);
  hPtSE->SetStats(0);
  hPtSE->SetTitle("p_{T} distribution; p_{T} (GeV/c); 1/N dN/dp_{T}");
  c1->cd(1);
  hPtSE->Draw();
  hPtME->Draw("same");

  TLegend* legend1 = new TLegend(0.8, 0.8, 0.9, 0.9);
  legend1->AddEntry(hPtSE, "SE", "lp");
  legend1->AddEntry(hPtME, "ME", "lp");
  legend1->Draw();

  c1->cd(3);
  TH1D* hPtRatio = (TH1D*)hPtSE->Clone("hPtRatio");
  hPtRatio->Divide(hPtME);
  hPtRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hPtRatio->SetTitle("SE/ME ratio; p_{T} (GeV/c); SE/ME");
  hPtRatio->Draw();

  c1->cd(2);
  TH1D *hYSE = hPtYSE->ProjectionY("hYSE");
  TH1D *hYME = hPtYME->ProjectionY("hYME");
  hYSE->Scale(1.0/hYSE->Integral());
  hYME->Scale(1.0/hYME->Integral());
  hYSE->SetLineColor(kRed);
  hYME->SetLineColor(kBlue);
  hYSE->SetLineWidth(4);
  hYME->SetLineWidth(4);
  hYSE->SetStats(0);
  hYSE->SetTitle("Y distribution; Y; 1/N dN/dY");
  hYSE->GetYaxis()->SetTitleOffset(1.0);
  hYSE->Draw();
  hYME->Draw("same");

  TLegend* legend2 = new TLegend(0.8, 0.8, 0.9, 0.9);
  legend2->AddEntry(hYSE, "SE", "lp");
  legend2->AddEntry(hYME, "ME", "lp");
  legend2->Draw();

  c1->cd(4);
  TH1D* hYRatio = (TH1D*)hYSE->Clone("hYRatio");
  hYRatio->Divide(hYME);
  hYRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hYRatio->SetTitle("SE/ME ratio; Y; SE/ME");
  hYRatio->Draw();



  c1->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Pt_Y_Distribution_SE_ME.png");

  TCanvas* c2 = new  TCanvas("c2","c2", 1800, 1400);
  c2->Divide(2,2);

  c2->cd(1);
  TH2D *hPhiEtaSE = (TH2D*)f->Get("lambdaspincorrderived/hPhiEtaSame;1");
  TH2D *hPhiEtaME = (TH2D*)f->Get("lambdaspincorrderived/hPhiEtaMix;1");

  TH1D *hPhiSE = hPhiEtaSE->ProjectionX("hPhiEtaSE_proj");
  TH1D *hPhiME = hPhiEtaME->ProjectionX("hPhiEtaME_proj");
  hPhiSE->Scale(1.0/hPhiSE->Integral());
  hPhiME->Scale(1.0/hPhiME->Integral());
  hPhiSE->SetLineColor(kRed);
  hPhiME->SetLineColor(kBlue);
  hPhiSE->SetLineWidth(4);
  hPhiME->SetLineWidth(4);
  hPhiSE->SetStats(0);
  hPhiSE->SetTitle("#phi distribution; #phi (rad); 1/N dN/d#phi");
  hPhiSE->GetXaxis()->SetTitle("#phi");
  hPhiSE->GetYaxis()->SetTitleOffset(1.0);
  hPhiSE->Draw();
  hPhiME->Draw("same");

  TLegend* legend3 = new TLegend(0.8, 0.8, 0.9, 0.9);
  legend3->AddEntry(hPhiSE, "SE", "lp");
  legend3->AddEntry(hPhiME, "ME", "lp");
  legend3->Draw();

 

  c2->cd(3);
  TH1D* hPhiRatio = (TH1D*)hPhiSE->Clone("hPhiRatio");
  hPhiRatio->Divide(hPhiME);
  hPhiRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hPhiRatio->SetTitle("SE/ME ratio; #phi (rad); SE/ME");
  hPhiRatio->Draw();  
  
  c2->cd(2);
  TH1D *hEtaSE = hPhiEtaSE->ProjectionY("hEtaSE_proj");
  TH1D *hEtaME = hPhiEtaME->ProjectionY("hEtaME_proj");
  hEtaSE->Scale(1.0/hEtaSE->Integral());
  hEtaME->Scale(1.0/hEtaME->Integral());
  hEtaSE->SetLineColor(kRed);
  hEtaME->SetLineColor(kBlue);
  hEtaSE->SetLineWidth(4);
  hEtaME->SetLineWidth(4);
  hEtaSE->SetStats(0);
  hEtaSE->SetTitle("#eta distribution; #eta; 1/N dN/d#eta");
  hEtaSE->GetYaxis()->SetTitleOffset(1.0);
  hEtaSE->Draw();
  hEtaME->Draw("same");

  TLegend* legend4 = new TLegend(0.8 ,0.8, 0.9, 0.9);
  legend4->AddEntry(hEtaSE, "SE", "lp");
  legend4->AddEntry(hEtaME, "ME", "lp");
  legend4->Draw();


  

  c2->cd(4);
  TH1D* hEtaRatio = (TH1D*)hEtaSE->Clone("hEtaRatio");
  hEtaRatio->Divide(hEtaME);
  hEtaRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hEtaRatio->SetTitle("SE/ME ratio; #eta; SE/ME");
  hEtaRatio->Draw();  

  c2->SaveAs("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Picture/Phi_Eta_Distribution_SE_ME.png");
  
  delete c1;
  delete c2;

}



int main(){
  //void Result2(int centlow = 0, int centhigh = 100, int nRbins = 6, TString inputfold = "Output_pp",
  //const double* rEdges = nullptr, int sigstat = 2, int bkgstat = 0, bool do_YieldSS = true, bool do_YieldSS_BCD = false)
  
  //do 2D fit
  Result2(0, 100, 6, "Output_pp", nullptr, 0, 2, false, true);
  //Result2(0, 100, 6, "Output_pp", nullptr, 0, 2, true, false);
  //do ABCD
  bool do_fit = false;
  if(do_fit){
    for(int i=0;i<=2;i++){
      Result2(0, 100, 6, "Output_pp", nullptr, 0, i, false, true);
      Result2(0, 100, 6, "Output_pp", nullptr, 0, i, true, false);
    }
  }
  
  bool do_Uncertianty = false;
  if(do_Uncertianty){
    UncertaintyOfSS(false); //US
    UncertaintyOfSS(true); // SS
  }
 
  bool do_PolarizationCorr = false;
  if(do_PolarizationCorr){
    drawUncertaintyPlots();
    drawPloarizationcorrPlots();
  }
  
 
 //SEOverME(true);
  //SEOverME(false);

  //x_P_Lambda_y();
  //draw_STAR_x_P_Lambda_y();
  //test_draw_STAR_x_P_Lambda_y();
  //Data_Eventnumber();
  //plot_thn_mass2D();

  bool do_SEOverME_Bkg_SandBandMethod = true;
  if(do_SEOverME_Bkg_SandBandMethod){
    SEOverMEForBkgwitSandBandMethod(false);
    SEOverMEForBkgwitSandBandMethod(true);
  }

  CalculationOfSpinCorrelation();

 //drawPloarizationcorrPlots_And_Simulation();

 //readTheEffSigma();
  
  return 0; 
}