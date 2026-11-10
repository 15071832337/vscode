// BuildUploadWeightAndQA_ULcombined.C
// One-step ROOT macro for stage-1 replacement residual weights.
// UL is built as LAL+ALL by default, matching the O2 task that applies one combined UL map.
//
// One-step ROOT macro for stage-1 replacement residual weights:
//   1) read SE/TGT and ME/REP TH3 maps from AnalysisResults.root
//   2) build the shape-only residual map
//        epsilon(phi,y,pT) = (ME/intME)/(SE/intSE)
//                           = (REP/intREP)/(TGT/intTGT)
//   3) write six CCDB-ready ROOT files with TH3D named "ccdb_object"
//   4) optionally upload them to CCDB
//   5) make shape-normalized 2D and 1D QA plots before and after weighting
//
// Weight modes:
//   useRaw3D=true:
//     epsilon(phi,y,pT) = [(ME/intME)/(SE/intSE)]_3D
//   useRaw3D=false:
//     epsilon(phi,y,pT) = epsilon_phi(phi) * epsilon_ypt(y,pT)
//     where both factors are separately shape-normalized.
//
// Measurement acceptance:
//   The histogram axes may be wider than the analysis acceptance. By default,
//   only phi=[0,2pi), y/eta=[-0.5,0.5), and pT=[0.5,3.0) enter the weight,
//   normalization, diagnostics, and QA. The CCDB epsilon is exactly 1 outside.
//
// Histogram convention assumed in the task:
//   x = phi, y = y/eta, z = pT
//
// Categories:
//   LL   = Lambda-Lambda
//   UL   = unlike, built as LAL+ALL by default
//   ALAL = AntiLambda-AntiLambda
//
// Default is a safe dry run: doUpload=false.
// To upload, pass doUpload=true.
//
// Example, fixed leg1 / replace leg2, factorized weights, upload false:
//   root -l -b -q 'BuildUploadWeightAndQA_shapeNormalized.C+("AnalysisResults.root","lambdaspincorrderived","Stage1Weights_factorized_ULcombined",1,"Users/s/skundu/2026_ppspincorrMC/stage1residual_factorized_ULcombined","https://alice-ccdb.cern.ch",false,false,true,0.01,10.0,2.0)'
//
// Same, with upload:
//   root -l -b -q 'BuildUploadWeightAndQA_shapeNormalized.C+("AnalysisResults.root","lambdaspincorrderived","Stage1Weights_factorized_ULcombined",1,"Users/s/skundu/2026_ppspincorrMC/stage1residual_factorized_ULcombined","https://alice-ccdb.cern.ch",true,false,true,0.01,10.0,2.0)'

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <iomanip>
#include <memory>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <initializer_list>

#include "CCDB/CcdbApi.h"

namespace weightqa {

static constexpr long kDefaultStartValidityMS = 1L;
static constexpr long kDefaultEndValidityMS   = 2524604400000L; // 2050-01-01

static bool good(double x) { return std::isfinite(x); }

struct Spec {
  TString tag;  // output/CCDB tag, e.g. REP_LL_leg2
  TString cat;  // LL, UL, ALAL
  int leg;      // 1 or 2
};

static const std::vector<Spec> kSpecs = {
  {"REP_LL_leg1",   "LL",   1},
  {"REP_UL_leg1",   "UL",   1},
  {"REP_ALAL_leg1", "ALAL", 1},
  {"REP_LL_leg2",   "LL",   2},
  {"REP_UL_leg2",   "UL",   2},
  {"REP_ALAL_leg2", "ALAL", 2}
};

static TObject* findObjectRecursive(TDirectory* dir, const TString& name)
{
  if (!dir) return nullptr;
  if (auto* obj = dir->Get(name)) return obj;

  TIter next(dir->GetListOfKeys());
  while (auto* keyObj = next()) {
    auto* key = dynamic_cast<TKey*>(keyObj);
    if (!key) continue;
    TObject* obj = key->ReadObj();
    auto* sub = dynamic_cast<TDirectory*>(obj);
    if (!sub) continue;
    TObject* found = findObjectRecursive(sub, name);
    if (found) return found;
  }
  return nullptr;
}

static TH3* getH3(TDirectory* taskDir, const TString& name, bool warn = true)
{
  auto* obj = findObjectRecursive(taskDir, name);
  auto* h = dynamic_cast<TH3*>(obj);
  if (!h && warn) std::cerr << "[WARN] missing TH3: " << name << "\n";
  return h;
}

static std::vector<double> getAxisEdges(const TAxis* axis)
{
  std::vector<double> edges(axis->GetNbins() + 1);
  for (int ib = 1; ib <= axis->GetNbins(); ++ib) {
    edges[ib - 1] = axis->GetBinLowEdge(ib);
  }
  edges[axis->GetNbins()] = axis->GetBinUpEdge(axis->GetNbins());
  return edges;
}

static TH3D* cloneH3D(const TH3* h, const TString& name)
{
  if (!h) return nullptr;

  // Always copy the physical bin edges. Reconstructing only from xmin/xmax
  // silently changes variable-width axes into uniform axes.
  const auto xEdges = getAxisEdges(h->GetXaxis());
  const auto yEdges = getAxisEdges(h->GetYaxis());
  const auto zEdges = getAxisEdges(h->GetZaxis());
  auto* out = new TH3D(name, h->GetTitle(),
                       h->GetNbinsX(), xEdges.data(),
                       h->GetNbinsY(), yEdges.data(),
                       h->GetNbinsZ(), zEdges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();

  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
        out->SetBinContent(ix, iy, iz, h->GetBinContent(ix, iy, iz));
        out->SetBinError(ix, iy, iz, h->GetBinError(ix, iy, iz));
      }
    }
  }
  return out;
}

static double integral3D(const TH3* h)
{
  if (!h) return 0.0;
  const double s = h->Integral(1, h->GetNbinsX(), 1, h->GetNbinsY(), 1, h->GetNbinsZ());
  return good(s) ? s : 0.0;
}

struct Acceptance {
  double phiMin;
  double phiMax;
  double yMin;
  double yMax;
  double ptMin;
  double ptMax;

  Acceptance(double inPhiMin = 0.0,
             double inPhiMax = 6.283185307179586,
             double inYMin = -0.5,
             double inYMax = 0.5,
             double inPtMin = 0.5,
             double inPtMax = 3.0)
    : phiMin(inPhiMin),
      phiMax(inPhiMax),
      yMin(inYMin),
      yMax(inYMax),
      ptMin(inPtMin),
      ptMax(inPtMax)
  {
  }
};

static bool insideAcceptance(const TH3* h, int ix, int iy, int iz,
                             const Acceptance& acc)
{
  if (!h) return false;
  const double phi = h->GetXaxis()->GetBinCenter(ix);
  const double y = h->GetYaxis()->GetBinCenter(iy);
  const double pt = h->GetZaxis()->GetBinCenter(iz);
  return phi >= acc.phiMin && phi < acc.phiMax &&
         y >= acc.yMin && y < acc.yMax &&
         pt >= acc.ptMin && pt < acc.ptMax;
}

static double integral3DInAcceptance(const TH3* h, const Acceptance& acc)
{
  if (!h) return 0.0;
  double sum = 0.0;
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
        if (!insideAcceptance(h, ix, iy, iz, acc)) continue;
        const double value = h->GetBinContent(ix, iy, iz);
        if (good(value)) sum += value;
      }
    }
  }
  return sum;
}

static TH3D* maskOutsideAcceptance(const TH3* h, const TString& name,
                                   const Acceptance& acc)
{
  TH3D* out = cloneH3D(h, name);
  if (!out) return nullptr;
  for (int ix = 1; ix <= out->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= out->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= out->GetNbinsZ(); ++iz) {
        if (insideAcceptance(out, ix, iy, iz, acc)) continue;
        out->SetBinContent(ix, iy, iz, 0.0);
        out->SetBinError(ix, iy, iz, 0.0);
      }
    }
  }
  return out;
}

static bool sameAxis(const TAxis* a, const TAxis* b)
{
  if (!a || !b) return false;
  if (a->GetNbins() != b->GetNbins()) return false;
  for (int ib = 1; ib <= a->GetNbins(); ++ib) {
    if (std::abs(a->GetBinLowEdge(ib) - b->GetBinLowEdge(ib)) > 1e-9) return false;
  }
  if (std::abs(a->GetBinUpEdge(a->GetNbins()) -
               b->GetBinUpEdge(b->GetNbins())) > 1e-9) return false;
  return true;
}

static bool sameBinning(const TH3* a, const TH3* b)
{
  return a && b &&
         sameAxis(a->GetXaxis(), b->GetXaxis()) &&
         sameAxis(a->GetYaxis(), b->GetYaxis()) &&
         sameAxis(a->GetZaxis(), b->GetZaxis());
}

static void fillUnity(TH3D* h)
{
  if (!h) return;
  h->Reset("ICES");
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
        h->SetBinContent(ix, iy, iz, 1.0);
        h->SetBinError(ix, iy, iz, 0.0);
      }
    }
  }
}

static TH3D* makeTemplate(TDirectory* taskDir)
{
  const std::vector<TString> tryNames = {
    "TGT_LL_leg2", "REP_LL_leg2", "TGT_UL_leg2", "REP_UL_leg2",
    "TGT_LAL_leg2", "REP_LAL_leg2", "TGT_ALL_leg2", "REP_ALL_leg2",
    "TGT_ALAL_leg2", "REP_ALAL_leg2",
    "TGT_LL_leg1", "REP_LL_leg1", "TGT_UL_leg1", "REP_UL_leg1",
    "TGT_LAL_leg1", "REP_LAL_leg1", "TGT_ALL_leg1", "REP_ALL_leg1",
    "TGT_ALAL_leg1", "REP_ALAL_leg1"
  };

  for (const auto& n : tryNames) {
    if (auto* h = getH3(taskDir, n, false)) {
      TH3D* out = cloneH3D(h, "unity_template");
      fillUnity(out);
      return out;
    }
  }
  return nullptr;
}

static TH3D* getDirect(TDirectory* taskDir, const TString& kind, const TString& cat, int leg)
{
  const TString name = Form("%s_%s_leg%d", kind.Data(), cat.Data(), leg);
  return cloneH3D(getH3(taskDir, name, false), name + "_clone");
}

static TH3D* sumLALALL(TDirectory* taskDir, const TString& kind, int leg)
{
  std::unique_ptr<TH3D> lal(getDirect(taskDir, kind, "LAL", leg));
  std::unique_ptr<TH3D> all(getDirect(taskDir, kind, "ALL", leg));
  if (!lal && !all) return nullptr;

  TH3D* out = lal ? dynamic_cast<TH3D*>(lal->Clone(Form("%s_UL_leg%d_LALplusALL", kind.Data(), leg)))
                  : dynamic_cast<TH3D*>(all->Clone(Form("%s_UL_leg%d_LALplusALL", kind.Data(), leg)));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  out->Sumw2();
  if (lal && all) out->Add(all.get());
  out->SetTitle(Form("%s UL leg%d = LAL + ALL", kind.Data(), leg));
  return out;
}

static TH3D* getCat(TDirectory* taskDir, const TString& kind, const TString& cat, int leg,
                    bool forceULFromLALALL = true)
{
  // The O2 task has one combined UL CCDB map. Therefore, by default, build UL as LAL+ALL.
  // Use direct TGT_UL/REP_UL only as a fallback, or if forceULFromLALALL=false.
  if (cat == "UL" && forceULFromLALALL) {
    if (auto* ul = sumLALALL(taskDir, kind, leg)) return ul;
    std::cerr << "[WARN] " << kind << "_LAL/ALL_leg" << leg
              << " not found. Falling back to direct " << kind << "_UL_leg" << leg << " if present.\n";
  }

  if (auto* direct = getDirect(taskDir, kind, cat, leg)) return direct;

  if (cat == "UL") {
    if (auto* ul = sumLALALL(taskDir, kind, leg)) return ul;
  }

  return nullptr;
}

static double capEpsilon(double e, double epsMin, double epsMax, int& nLow, int& nHigh)
{
  if (!good(e) || e <= 0.0) return 1.0;
  if (e < epsMin) { ++nLow; return epsMin; }
  if (e > epsMax) { ++nHigh; return epsMax; }
  return e;
}

static TH3D* buildEpsilonMEoverSE(const TH3* hSE,
                                  const TH3* hME,
                                  const TH3* unityTemplate,
                                  const TString& tag,
                                  bool useRaw3D,
                                  double epsMin,
                                  double epsMax,
                                  double minSEContent,
                                  const Acceptance& acc,
                                  FILE* summary)
{
  TH3D* eps = nullptr;
  if (hSE) eps = cloneH3D(hSE, "ccdb_object");
  else if (unityTemplate) eps = cloneH3D(unityTemplate, "ccdb_object");

  if (!eps) {
    std::cerr << "[ERROR] no template for " << tag << "\n";
    return nullptr;
  }

  eps->SetName("ccdb_object");
  eps->SetTitle(tag + ": shape epsilon=(ME/intME)/(SE/intSE);#phi;y/#eta;p_{T}");
  fillUnity(eps);

  const double intSE = integral3DInAcceptance(hSE, acc);
  const double intME = integral3DInAcceptance(hME, acc);
  const double global = (intSE > 0.0) ? intME / intSE : 1.0;

  if (!hSE || !hME || !sameBinning(hSE, hME) || !(intSE > 0.0) || !(intME > 0.0)) {
    std::cerr << "[WARN] " << tag
              << ": missing/empty/mismatched SE/ME inside the acceptance. Writing unity.\n";
    if (summary) fprintf(summary, "%s,unity,%.17g,%.17g,1,1,1,1,0,0,0,0,0,0,0\n", tag.Data(), intSE, intME);
    return eps;
  }

  int nDefault = 0, nCapLow = 0, nCapHigh = 0;
  int nZeroTGTWithREP = 0, nZeroREPWithTGT = 0;
  int nInsideBins = 0;
  for (int ix = 1; ix <= eps->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= eps->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= eps->GetNbinsZ(); ++iz) {
        if (!insideAcceptance(eps, ix, iy, iz, acc)) continue;
        ++nInsideBins;
        const double se = hSE->GetBinContent(ix, iy, iz);
        const double me = hME->GetBinContent(ix, iy, iz);
        if (se <= 0.0 && me > 0.0) ++nZeroTGTWithREP;
        if (se > 0.0 && me <= 0.0) ++nZeroREPWithTGT;
      }
    }
  }
  const int nAllBins = eps->GetNbinsX() * eps->GetNbinsY() * eps->GetNbinsZ();
  const int nOutsideBins = nAllBins - nInsideBins;
  double minFinal = 1e30, maxFinal = -1e30;

  if (useRaw3D) {
    for (int ix = 1; ix <= eps->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= eps->GetNbinsY(); ++iy) {
        for (int iz = 1; iz <= eps->GetNbinsZ(); ++iz) {
          if (!insideAcceptance(eps, ix, iy, iz, acc)) {
            eps->SetBinContent(ix, iy, iz, 1.0);
            eps->SetBinError(ix, iy, iz, 0.0);
            continue;
          }
          const double se = hSE->GetBinContent(ix, iy, iz);
          const double me = hME->GetBinContent(ix, iy, iz);
          double e = 1.0;
          if (good(se) && good(me) && se > 0.0 &&
              se >= minSEContent && me > 0.0 && global > 0.0) {
            // Shape-only map: remove the global REP/TGT normalization.
            e = (me / se) / global;
          } else {
            ++nDefault;
          }
          e = capEpsilon(e, epsMin, epsMax, nCapLow, nCapHigh);
          eps->SetBinContent(ix, iy, iz, e);
          eps->SetBinError(ix, iy, iz, 0.0);
          minFinal = std::min(minFinal, e);
          maxFinal = std::max(maxFinal, e);
        }
      }
    }
  } else {
    // Factorized shape-only fallback:
    // epsilon(phi,y,pT) = epsilon_phi_shape(phi) * epsilon_ypt_shape(y,pT).
    const int nx = eps->GetNbinsX();
    const int ny = eps->GetNbinsY();
    const int nz = eps->GetNbinsZ();
    std::vector<double> sePhi(nx + 1, 0.0), mePhi(nx + 1, 0.0);
    std::vector<double> seYpt((ny + 1) * (nz + 1), 0.0);
    std::vector<double> meYpt((ny + 1) * (nz + 1), 0.0);
    auto yptIndex = [nz](int iy, int iz) { return iy * (nz + 1) + iz; };

    for (int ix = 1; ix <= nx; ++ix) {
      for (int iy = 1; iy <= ny; ++iy) {
        for (int iz = 1; iz <= nz; ++iz) {
          if (!insideAcceptance(eps, ix, iy, iz, acc)) continue;
          const double se = hSE->GetBinContent(ix, iy, iz);
          const double me = hME->GetBinContent(ix, iy, iz);
          sePhi[ix] += se;
          mePhi[ix] += me;
          seYpt[yptIndex(iy, iz)] += se;
          meYpt[yptIndex(iy, iz)] += me;
        }
      }
    }

    for (int ix = 1; ix <= nx; ++ix) {
      const double seP = sePhi[ix];
      const double meP = mePhi[ix];
      const double ePhi = (seP > 0.0 && seP >= minSEContent &&
                           meP > 0.0 && global > 0.0)
                            ? (meP / seP) / global
                            : 1.0;

      for (int iy = 1; iy <= ny; ++iy) {
        for (int iz = 1; iz <= nz; ++iz) {
          if (!insideAcceptance(eps, ix, iy, iz, acc)) {
            eps->SetBinContent(ix, iy, iz, 1.0);
            eps->SetBinError(ix, iy, iz, 0.0);
            continue;
          }
          const double seYP = seYpt[yptIndex(iy, iz)];
          const double meYP = meYpt[yptIndex(iy, iz)];
          const double eYpt = (seYP > 0.0 && seYP >= minSEContent &&
                               meYP > 0.0 && global > 0.0)
                                ? (meYP / seYP) / global : 1.0;
          double e = ePhi * eYpt;
          if (hSE->GetBinContent(ix, iy, iz) < minSEContent ||
              hME->GetBinContent(ix, iy, iz) <= 0.0) {
            ++nDefault;
          }
          e = capEpsilon(e, epsMin, epsMax, nCapLow, nCapHigh);
          eps->SetBinContent(ix, iy, iz, e);
          eps->SetBinError(ix, iy, iz, 0.0);
          minFinal = std::min(minFinal, e);
          maxFinal = std::max(maxFinal, e);
        }
      }
    }
  }

  // Enforce a genuinely shape-only correction. In particular, the product of
  // two separately normalized factorized components need not itself preserve
  // the total yield when phi and (y,pT) are correlated. After this rescaling,
  // sum[ME/epsilon] = sum[ME].
  double correctedIntegral = 0.0;
  for (int ix = 1; ix <= eps->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= eps->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= eps->GetNbinsZ(); ++iz) {
        if (!insideAcceptance(eps, ix, iy, iz, acc)) continue;
        const double e = eps->GetBinContent(ix, iy, iz);
        const double me = hME->GetBinContent(ix, iy, iz);
        if (good(e) && e > 0.0 && good(me)) correctedIntegral += me / e;
      }
    }
  }
  const double shapeScale = (correctedIntegral > 0.0 && intME > 0.0)
                              ? correctedIntegral / intME : 1.0;
  minFinal = 1e30;
  maxFinal = -1e30;
  for (int ix = 1; ix <= eps->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= eps->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= eps->GetNbinsZ(); ++iz) {
        if (!insideAcceptance(eps, ix, iy, iz, acc)) {
          // Outside the measurement acceptance the CCDB map must be neutral.
          eps->SetBinContent(ix, iy, iz, 1.0);
          eps->SetBinError(ix, iy, iz, 0.0);
          continue;
        }
        const double e = eps->GetBinContent(ix, iy, iz) * shapeScale;
        eps->SetBinContent(ix, iy, iz, e);
        minFinal = std::min(minFinal, e);
        maxFinal = std::max(maxFinal, e);
      }
    }
  }

  if (summary) {
    fprintf(summary, "%s,built,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%d,%d\n",
            tag.Data(), intSE, intME, global, shapeScale, minFinal, maxFinal,
            nInsideBins, nOutsideBins, nDefault, nCapLow, nCapHigh,
            nZeroTGTWithREP, nZeroREPWithTGT);
  }

  std::cout << "[OK] " << tag
            << " SE=" << intSE
            << " ME=" << intME
            << " global ME/SE=" << global
            << " shape scale=" << shapeScale
            << " epsilon=" << minFinal << "--" << maxFinal
            << " inside/outside bins=" << nInsideBins << "/" << nOutsideBins
            << " defaults=" << nDefault
            << " caps=" << nCapLow << "/" << nCapHigh << "\n";
  if (nZeroTGTWithREP > 0 || nZeroREPWithTGT > 0) {
    std::cout << "     zero-support cells: TGT=0,REP>0: " << nZeroTGTWithREP
              << "; REP=0,TGT>0 (not recoverable by weighting): "
              << nZeroREPWithTGT << "\n";
  }

  return eps;
}

static TH3D* makeInverseWeight(const TH3* eps, const TString& name, double epsMin, double epsMax)
{
  if (!eps) return nullptr;
  TH3D* out = cloneH3D(eps, name);
  if (!out) return nullptr;
  out->Reset("ICES");
  const double wMin = (epsMax > 0.0) ? 1.0 / epsMax : 0.0;
  const double wMax = (epsMin > 0.0) ? 1.0 / epsMin : 1e12;
  for (int ix = 1; ix <= eps->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= eps->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= eps->GetNbinsZ(); ++iz) {
        const double e = eps->GetBinContent(ix, iy, iz);
        double w = (good(e) && e > 0.0) ? 1.0 / e : 1.0;
        if (w < wMin) w = wMin;
        if (w > wMax) w = wMax;
        out->SetBinContent(ix, iy, iz, w);
        out->SetBinError(ix, iy, iz, 0.0);
      }
    }
  }
  out->SetTitle("shape residual weight = 1/epsilon;#phi;y/#eta;p_{T}");
  return out;
}

static TH3D* makeCorrectedME(const TH3* hME, const TH3* eps, const TString& name)
{
  if (!hME || !eps || !sameBinning(hME, eps)) return nullptr;
  TH3D* out = cloneH3D(hME, name);
  if (!out) return nullptr;
  for (int ix = 1; ix <= out->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= out->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= out->GetNbinsZ(); ++iz) {
        const double e = eps->GetBinContent(ix, iy, iz);
        const double v = out->GetBinContent(ix, iy, iz);
        const double er = out->GetBinError(ix, iy, iz);
        if (good(e) && e > 0.0) {
          out->SetBinContent(ix, iy, iz, v / e);
          out->SetBinError(ix, iy, iz, er / e);
        } else {
          out->SetBinContent(ix, iy, iz, 0.0);
          out->SetBinError(ix, iy, iz, 0.0);
        }
      }
    }
  }
  return out;
}

static TH2D* projectPhiY(const TH3* h, const TString& name)
{
  if (!h) return nullptr;
  const TAxis* ax = h->GetXaxis(); // phi
  const TAxis* ay = h->GetYaxis(); // y/eta
  const auto xEdges = getAxisEdges(ax);
  const auto yEdges = getAxisEdges(ay);
  auto* out = new TH2D(name, ";#phi;y/#eta",
                       ax->GetNbins(), xEdges.data(),
                       ay->GetNbins(), yEdges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      double sum = 0.0, err2 = 0.0;
      for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
        const double v = h->GetBinContent(ix, iy, iz);
        const double e = h->GetBinError(ix, iy, iz);
        sum += v;
        err2 += e * e;
      }
      out->SetBinContent(ix, iy, sum);
      out->SetBinError(ix, iy, std::sqrt(err2));
    }
  }
  return out;
}

static TH2D* projectPtY(const TH3* h, const TString& name)
{
  if (!h) return nullptr;
  const TAxis* ay = h->GetYaxis(); // y/eta
  const TAxis* az = h->GetZaxis(); // pT
  const auto yEdges = getAxisEdges(ay);
  const auto zEdges = getAxisEdges(az);
  auto* out = new TH2D(name, ";p_{T} (GeV/#it{c});y/#eta",
                       az->GetNbins(), zEdges.data(),
                       ay->GetNbins(), yEdges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();
  for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      double sum = 0.0, err2 = 0.0;
      for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        const double v = h->GetBinContent(ix, iy, iz);
        const double e = h->GetBinError(ix, iy, iz);
        sum += v;
        err2 += e * e;
      }
      out->SetBinContent(iz, iy, sum);
      out->SetBinError(iz, iy, std::sqrt(err2));
    }
  }
  return out;
}

static TH1D* project1D(const TH3* h, int axis, const TString& name)
{
  if (!h || axis < 0 || axis > 2) return nullptr;

  const TAxis* sourceAxis = (axis == 0) ? h->GetXaxis() :
                            (axis == 1) ? h->GetYaxis() : h->GetZaxis();
  const auto edges = getAxisEdges(sourceAxis);

  const TString axisTitle = (axis == 0) ? "#phi" :
                            (axis == 1) ? "y/#eta" : "p_{T} (GeV/#it{c})";
  auto* out = new TH1D(name, Form(";%s;counts", axisTitle.Data()),
                       sourceAxis->GetNbins(), edges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();

  std::vector<double> error2(sourceAxis->GetNbins() + 1, 0.0);
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= h->GetNbinsZ(); ++iz) {
        const int ib = (axis == 0) ? ix : (axis == 1) ? iy : iz;
        const double value = h->GetBinContent(ix, iy, iz);
        const double error = h->GetBinError(ix, iy, iz);
        out->AddBinContent(ib, value);
        error2[ib] += error * error;
      }
    }
  }
  for (int ib = 1; ib <= out->GetNbinsX(); ++ib) {
    out->SetBinError(ib, std::sqrt(error2[ib]));
  }
  return out;
}

static TH2D* makeAbsShapeRatioMinus1(const TH2* num, const TH2* den,
                                     const TString& name, double minDen = 1e-12)
{
  if (!num || !den) return nullptr;
  if (num->GetNbinsX() != den->GetNbinsX() || num->GetNbinsY() != den->GetNbinsY()) {
    std::cerr << "[ERROR] 2D bin mismatch for " << name << "\n";
    return nullptr;
  }
  auto* out = dynamic_cast<TH2D*>(num->Clone(name));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  out->Reset("ICES");
  out->Sumw2();

  const double intNum = num->Integral(1, num->GetNbinsX(), 1, num->GetNbinsY());
  const double intDen = den->Integral(1, den->GetNbinsX(), 1, den->GetNbinsY());
  const double scale = (intNum > 0.0 && intDen > 0.0) ? intDen / intNum : 1.0;

  for (int ix = 1; ix <= out->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= out->GetNbinsY(); ++iy) {
      const double a = scale * num->GetBinContent(ix, iy);
      const double b = den->GetBinContent(ix, iy);
      double v = 0.0;
      if (good(a) && good(b) && b > minDen) v = std::abs(a / b - 1.0);
      else if (a > 0.0) v = 1.0;
      out->SetBinContent(ix, iy, good(v) ? v : 0.0);
      out->SetBinError(ix, iy, 0.0);
    }
  }
  out->SetMinimum(0.0);
  out->SetMaximum(1.0);
  out->GetZaxis()->SetTitle("|shape ratio - 1|");
  return out;
}

static TH1D* makeShapeRatio(const TH1* num, const TH1* den, const TString& name,
                            double minDen = 1e-12)
{
  if (!num || !den || num->GetNbinsX() != den->GetNbinsX()) return nullptr;
  auto* out = dynamic_cast<TH1D*>(num->Clone(name));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  out->Reset("ICES");
  out->Sumw2();

  const double intNum = num->Integral(1, num->GetNbinsX());
  const double intDen = den->Integral(1, den->GetNbinsX());
  const double scale = (intNum > 0.0 && intDen > 0.0) ? intDen / intNum : 1.0;

  for (int ib = 1; ib <= out->GetNbinsX(); ++ib) {
    const double a = scale * num->GetBinContent(ib);
    const double ea = scale * num->GetBinError(ib);
    const double b = den->GetBinContent(ib);
    const double eb = den->GetBinError(ib);

    if (good(a) && good(b) && b > minDen) {
      const double ratio = a / b;
      double error2 = (ea / b) * (ea / b);
      if (eb > 0.0) error2 += (a * eb / (b * b)) * (a * eb / (b * b));
      out->SetBinContent(ib, ratio);
      out->SetBinError(ib, std::sqrt(std::max(0.0, error2)));
    } else {
      out->SetBinContent(ib, 0.0);
      out->SetBinError(ib, 0.0);
    }
  }
  out->GetYaxis()->SetTitle("shape ratio");
  return out;
}

struct Stats {
  double mean = 0.0;
  double rms = 0.0;
  double max = 0.0;
  int n = 0;
};

static Stats calcStats(const TH2* h, double xMin, double xMax,
                       double yMin, double yMax)
{
  Stats s;
  if (!h) return s;
  double sum = 0.0, sum2 = 0.0;
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      const double x = h->GetXaxis()->GetBinCenter(ix);
      const double y = h->GetYaxis()->GetBinCenter(iy);
      if (x < xMin || x >= xMax || y < yMin || y >= yMax) continue;
      const double v = h->GetBinContent(ix, iy);
      if (!good(v)) continue;
      sum += v;
      sum2 += v * v;
      s.max = std::max(s.max, v);
      ++s.n;
    }
  }
  if (s.n > 0) {
    s.mean = sum / s.n;
    s.rms = std::sqrt(std::max(0.0, sum2 / s.n - s.mean * s.mean));
  }
  return s;
}

static Stats calcStats1D(const TH1* h, double xMin, double xMax)
{
  Stats s;
  if (!h) return s;
  double sum = 0.0, sum2 = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double x = h->GetXaxis()->GetBinCenter(ib);
    if (x < xMin || x >= xMax) continue;
    const double ratio = h->GetBinContent(ib);
    if (!good(ratio) || ratio <= 0.0) continue;
    const double v = std::abs(ratio - 1.0);
    sum += v;
    sum2 += v * v;
    s.max = std::max(s.max, v);
    ++s.n;
  }
  if (s.n > 0) {
    s.mean = sum / s.n;
    s.rms = std::sqrt(std::max(0.0, sum2 / s.n - s.mean * s.mean));
  }
  return s;
}

static void style2D(TH2* h, const TString& title)
{
  if (!h) return;
  h->SetTitle(title);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetZaxis()->SetTitleOffset(1.25);
}

static void drawLabel(const TString& txt)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.12, 0.93, txt);
}

static void makeQAPlot(const TString& outPng2D,
                       const TString& outPng1D,
                       const TString& tag,
                       const TH3* hSE,
                       const TH3* hME,
                       const TH3* hEps,
                       const Acceptance& acc,
                       FILE* qacsv,
                       TDirectory* outDir = nullptr)
{
  if (!hSE || !hME || !hEps) return;


  

  // QA uses only the declared measurement acceptance. Histogram bins outside
  // it are intentionally ignored rather than treated as failed closure bins.
  std::unique_ptr<TH3D> hSEAcc(maskOutsideAcceptance(hSE, tag + "_SE_acceptance", acc));
  std::unique_ptr<TH3D> hMEAcc(maskOutsideAcceptance(hME, tag + "_ME_acceptance", acc));
  std::unique_ptr<TH3D> hMEw(makeCorrectedME(hMEAcc.get(), hEps, tag + "_ME_weighted_3D"));

  std::unique_ptr<TH2D> sePhiY(projectPhiY(hSEAcc.get(), tag + "_SE_phiY"));
  std::unique_ptr<TH2D> mePhiY(projectPhiY(hMEAcc.get(), tag + "_ME_phiY"));
  std::unique_ptr<TH2D> mewPhiY(projectPhiY(hMEw.get(), tag + "_MEw_phiY"));

  std::unique_ptr<TH2D> sePtY(projectPtY(hSEAcc.get(), tag + "_SE_ptY"));
  std::unique_ptr<TH2D> mePtY(projectPtY(hMEAcc.get(), tag + "_ME_ptY"));
  std::unique_ptr<TH2D> mewPtY(projectPtY(hMEw.get(), tag + "_MEw_ptY"));

  std::unique_ptr<TH2D> rawPhiY(makeAbsShapeRatioMinus1(mePhiY.get(), sePhiY.get(), tag + "_raw_abs_phiY"));
  std::unique_ptr<TH2D> rawPtY(makeAbsShapeRatioMinus1(mePtY.get(), sePtY.get(), tag + "_raw_abs_ptY"));
  std::unique_ptr<TH2D> corrPhiY(makeAbsShapeRatioMinus1(mewPhiY.get(), sePhiY.get(), tag + "_corr_abs_phiY"));
  std::unique_ptr<TH2D> corrPtY(makeAbsShapeRatioMinus1(mewPtY.get(), sePtY.get(), tag + "_corr_abs_ptY"));

  style2D(rawPhiY.get(),  tag + ": before, |normalized ME/SE - 1|;#phi;y/#eta");
  style2D(rawPtY.get(),   tag + ": before, |normalized ME/SE - 1|;p_{T};y/#eta");
  style2D(corrPhiY.get(), tag + ": after, |normalized ME_{w}/SE - 1|;#phi;y/#eta");
  style2D(corrPtY.get(),  tag + ": after, |normalized ME_{w}/SE - 1|;p_{T};y/#eta");
  for (TH2* h : {static_cast<TH2*>(rawPhiY.get()), static_cast<TH2*>(corrPhiY.get())}) {
    if (h) {
      h->GetXaxis()->SetRangeUser(acc.phiMin, acc.phiMax);
      h->GetYaxis()->SetRangeUser(acc.yMin, acc.yMax);
    }
  }
  for (TH2* h : {static_cast<TH2*>(rawPtY.get()), static_cast<TH2*>(corrPtY.get())}) {
    if (h) {
      h->GetXaxis()->SetRangeUser(acc.ptMin, acc.ptMax);
      h->GetYaxis()->SetRangeUser(acc.yMin, acc.yMax);
    }
  }

  TCanvas c2D("cQA2D", "cQA2D", 1500, 1200);
  c2D.Divide(2, 2);
  TH2* arr[4] = {rawPhiY.get(), rawPtY.get(), corrPhiY.get(), corrPtY.get()};
  const char* lab[4] = {"before: #phi vs y", "before: p_{T} vs y", "after: #phi vs y", "after: p_{T} vs y"};
  for (int i = 0; i < 4; ++i) {
    c2D.cd(i + 1);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    if (arr[i]) arr[i]->Draw("COLZ");
    drawLabel(lab[i]);
  }
  c2D.SaveAs(outPng2D);

  std::unique_ptr<TH1D> sePhi(project1D(hSEAcc.get(), 0, tag + "_SE_phi"));
  std::unique_ptr<TH1D> mePhi(project1D(hMEAcc.get(), 0, tag + "_ME_phi"));
  std::unique_ptr<TH1D> mewPhi(project1D(hMEw.get(), 0, tag + "_MEw_phi"));
  std::unique_ptr<TH1D> seY(project1D(hSEAcc.get(), 1, tag + "_SE_y"));
  std::unique_ptr<TH1D> meY(project1D(hMEAcc.get(), 1, tag + "_ME_y"));
  std::unique_ptr<TH1D> mewY(project1D(hMEw.get(), 1, tag + "_MEw_y"));
  std::unique_ptr<TH1D> sePt(project1D(hSEAcc.get(), 2, tag + "_SE_pt"));
  std::unique_ptr<TH1D> mePt(project1D(hMEAcc.get(), 2, tag + "_ME_pt"));
  std::unique_ptr<TH1D> mewPt(project1D(hMEw.get(), 2, tag + "_MEw_pt"));

  std::unique_ptr<TH1D> rawPhi(makeShapeRatio(mePhi.get(), sePhi.get(), tag + "_raw_shape_phi"));
  std::unique_ptr<TH1D> rawY(makeShapeRatio(meY.get(), seY.get(), tag + "_raw_shape_y"));
  std::unique_ptr<TH1D> rawPt(makeShapeRatio(mePt.get(), sePt.get(), tag + "_raw_shape_pt"));
  std::unique_ptr<TH1D> corrPhi(makeShapeRatio(mewPhi.get(), sePhi.get(), tag + "_corr_shape_phi"));
  std::unique_ptr<TH1D> corrY(makeShapeRatio(mewY.get(), seY.get(), tag + "_corr_shape_y"));
  std::unique_ptr<TH1D> corrPt(makeShapeRatio(mewPt.get(), sePt.get(), tag + "_corr_shape_pt"));



  // ========== 创建 Ratio 直方图 (修正前: ME/SE) ==========
// 注意：先归一化再计算 ratio

  // 1. 复制并归一化 mePhi 和 sePhi
  std::unique_ptr<TH1D> mePhi_norm((TH1D*)mePhi->Clone(Form("%s_ME_phi_norm", tag.Data())));
  std::unique_ptr<TH1D> sePhi_norm((TH1D*)sePhi->Clone(Form("%s_SE_phi_norm", tag.Data())));
  if (mePhi_norm && sePhi_norm) {
      if (mePhi_norm->Integral() > 0) mePhi_norm->Scale(1.0 / mePhi_norm->Integral());
      if (sePhi_norm->Integral() > 0) sePhi_norm->Scale(1.0 / sePhi_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioPhi_before((TH1D*)mePhi_norm->Clone(Form("%s_ratio_phi_before", tag.Data())));
  if (ratioPhi_before) {
      ratioPhi_before->Divide(sePhi_norm.get());
      ratioPhi_before->SetTitle(Form("%s: ME/SE ratio vs #phi (before weight);#phi;ME/SE", tag.Data()));
      ratioPhi_before->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  // 2. Y轴 (y/eta)
  std::unique_ptr<TH1D> meY_norm((TH1D*)meY->Clone(Form("%s_ME_y_norm", tag.Data())));
  std::unique_ptr<TH1D> seY_norm((TH1D*)seY->Clone(Form("%s_SE_y_norm", tag.Data())));
  if (meY_norm && seY_norm) {
      if (meY_norm->Integral() > 0) meY_norm->Scale(1.0 / meY_norm->Integral());
      if (seY_norm->Integral() > 0) seY_norm->Scale(1.0 / seY_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioY_before((TH1D*)meY_norm->Clone(Form("%s_ratio_y_before", tag.Data())));
  if (ratioY_before) {
      ratioY_before->Divide(seY_norm.get());
      ratioY_before->SetTitle(Form("%s: ME/SE ratio vs y (before weight);y;ME/SE", tag.Data()));
      ratioY_before->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  // 3. pT
  std::unique_ptr<TH1D> mePt_norm((TH1D*)mePt->Clone(Form("%s_ME_pt_norm", tag.Data())));
  std::unique_ptr<TH1D> sePt_norm((TH1D*)sePt->Clone(Form("%s_SE_pt_norm", tag.Data())));
  if (mePt_norm && sePt_norm) {
      if (mePt_norm->Integral() > 0) mePt_norm->Scale(1.0 / mePt_norm->Integral());
      if (sePt_norm->Integral() > 0) sePt_norm->Scale(1.0 / sePt_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioPt_before((TH1D*)mePt_norm->Clone(Form("%s_ratio_pt_before", tag.Data())));
  if (ratioPt_before) {
      ratioPt_before->Divide(sePt_norm.get());
      ratioPt_before->SetTitle(Form("%s: ME/SE ratio vs p_{T} (before weight);p_{T};ME/SE", tag.Data()));
      ratioPt_before->GetYaxis()->SetRangeUser(0.8, 1.2);
  }



  // ========== 创建 Ratio 直方图 (修正后: MEw/SE) ==========
  // 同样需要归一化

  // 1. phi
  std::unique_ptr<TH1D> mewPhi_norm((TH1D*)mewPhi->Clone(Form("%s_MEw_phi_norm", tag.Data())));
  if (mewPhi_norm) {
      if (mewPhi_norm->Integral() > 0) mewPhi_norm->Scale(1.0 / mewPhi_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioPhi_after((TH1D*)mewPhi_norm->Clone(Form("%s_ratio_phi_after", tag.Data())));
  if (ratioPhi_after) {
      ratioPhi_after->Divide(sePhi_norm.get());  // 使用已归一化的 sePhi_norm
      ratioPhi_after->SetTitle(Form("%s: MEw/SE ratio vs #phi (after weight);#phi;MEw/SE", tag.Data()));
      ratioPhi_after->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  // 2. y
  std::unique_ptr<TH1D> mewY_norm((TH1D*)mewY->Clone(Form("%s_MEw_y_norm", tag.Data())));
  if (mewY_norm) {
      if (mewY_norm->Integral() > 0) mewY_norm->Scale(1.0 / mewY_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioY_after((TH1D*)mewY_norm->Clone(Form("%s_ratio_y_after", tag.Data())));
  if (ratioY_after) {
      ratioY_after->Divide(seY_norm.get());  // 使用已归一化的 seY_norm
      ratioY_after->SetTitle(Form("%s: MEw/SE ratio vs y (after weight);y;MEw/SE", tag.Data()));
      ratioY_after->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  // 3. pT
  std::unique_ptr<TH1D> mewPt_norm((TH1D*)mewPt->Clone(Form("%s_MEw_pt_norm", tag.Data())));
  if (mewPt_norm) {
      if (mewPt_norm->Integral() > 0) mewPt_norm->Scale(1.0 / mewPt_norm->Integral());
  }

  std::unique_ptr<TH1D> ratioPt_after((TH1D*)mewPt_norm->Clone(Form("%s_ratio_pt_after", tag.Data())));
  if (ratioPt_after) {
      ratioPt_after->Divide(sePt_norm.get());  // 使用已归一化的 sePt_norm
      ratioPt_after->SetTitle(Form("%s: MEw/SE ratio vs p_{T} (after weight);p_{T};MEw/SE", tag.Data()));
      ratioPt_after->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  if (outDir) {
    outDir->cd();  // 切换到指定的目录
   
    if (sePhi) sePhi->Write();
    if (mePhi) mePhi->Write();
    if (mewPhi) mewPhi->Write();
    if (seY) seY->Write();
    if (meY) meY->Write();
    if (mewY) mewY->Write();
    if (sePt) sePt->Write();      // 这就是您要保存的一维pT分布
    if (mePt) mePt->Write();
    if (mewPt) mewPt->Write();


    // ========== 保存 Ratio 直方图 (修正前和修正后) ==========
    if (ratioPhi_before) ratioPhi_before->Write();
    if (ratioY_before) ratioY_before->Write();
    if (ratioPt_before) ratioPt_before->Write();
    
    if (ratioPhi_after) ratioPhi_after->Write();
    if (ratioY_after) ratioY_after->Write();
    if (ratioPt_after) ratioPt_after->Write();
    

  }

  TCanvas c1D("cQA1D", "cQA1D", 1800, 1000);
  c1D.Divide(3, 2);
  TH1* arr1D[6] = {rawPhi.get(), rawY.get(), rawPt.get(),
                   corrPhi.get(), corrY.get(), corrPt.get()};
  const char* lab1D[6] = {"before: #phi", "before: y/#eta", "before: p_{T}",
                          "after: #phi", "after: y/#eta", "after: p_{T}"};
  std::vector<std::unique_ptr<TLine>> unityLines;
  for (int i = 0; i < 6; ++i) {
    c1D.cd(i + 1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    if (arr1D[i]) {
      if (i == 0 || i == 3) arr1D[i]->GetXaxis()->SetRangeUser(acc.phiMin, acc.phiMax);
      if (i == 1 || i == 4) arr1D[i]->GetXaxis()->SetRangeUser(acc.yMin, acc.yMax);
      if (i == 2 || i == 5) arr1D[i]->GetXaxis()->SetRangeUser(acc.ptMin, acc.ptMax);
      arr1D[i]->SetTitle(Form("%s: %s", tag.Data(), lab1D[i]));
      arr1D[i]->SetMarkerStyle(20);
      arr1D[i]->SetMarkerSize(0.7);
      arr1D[i]->SetMinimum(0.0);
      arr1D[i]->SetMaximum(2.0);
      arr1D[i]->Draw("E1");
      unityLines.emplace_back(std::make_unique<TLine>(
        arr1D[i]->GetXaxis()->GetXmin(), 1.0,
        arr1D[i]->GetXaxis()->GetXmax(), 1.0));
      unityLines.back()->SetLineColor(kRed + 1);
      unityLines.back()->SetLineStyle(2);
      unityLines.back()->Draw();
    }
  }
  c1D.SaveAs(outPng1D);

  if (qacsv) {
    const Stats s0 = calcStats(rawPhiY.get(), acc.phiMin, acc.phiMax, acc.yMin, acc.yMax);
    const Stats s1 = calcStats(rawPtY.get(), acc.ptMin, acc.ptMax, acc.yMin, acc.yMax);
    const Stats s2 = calcStats(corrPhiY.get(), acc.phiMin, acc.phiMax, acc.yMin, acc.yMax);
    const Stats s3 = calcStats(corrPtY.get(), acc.ptMin, acc.ptMax, acc.yMin, acc.yMax);
    fprintf(qacsv, "%s,before,phiY,%.17g,%.17g,%.17g,%d\n", tag.Data(), s0.mean, s0.rms, s0.max, s0.n);
    fprintf(qacsv, "%s,before,ptY,%.17g,%.17g,%.17g,%d\n", tag.Data(), s1.mean, s1.rms, s1.max, s1.n);
    fprintf(qacsv, "%s,after,phiY,%.17g,%.17g,%.17g,%d\n", tag.Data(), s2.mean, s2.rms, s2.max, s2.n);
    fprintf(qacsv, "%s,after,ptY,%.17g,%.17g,%.17g,%d\n", tag.Data(), s3.mean, s3.rms, s3.max, s3.n);

    const Stats s4 = calcStats1D(rawPhi.get(), acc.phiMin, acc.phiMax);
    const Stats s5 = calcStats1D(rawY.get(), acc.yMin, acc.yMax);
    const Stats s6 = calcStats1D(rawPt.get(), acc.ptMin, acc.ptMax);
    const Stats s7 = calcStats1D(corrPhi.get(), acc.phiMin, acc.phiMax);
    const Stats s8 = calcStats1D(corrY.get(), acc.yMin, acc.yMax);
    const Stats s9 = calcStats1D(corrPt.get(), acc.ptMin, acc.ptMax);
    fprintf(qacsv, "%s,before,phi1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s4.mean, s4.rms, s4.max, s4.n);
    fprintf(qacsv, "%s,before,y1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s5.mean, s5.rms, s5.max, s5.n);
    fprintf(qacsv, "%s,before,pt1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s6.mean, s6.rms, s6.max, s6.n);
    fprintf(qacsv, "%s,after,phi1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s7.mean, s7.rms, s7.max, s7.n);
    fprintf(qacsv, "%s,after,y1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s8.mean, s8.rms, s8.max, s8.n);
    fprintf(qacsv, "%s,after,pt1D,%.17g,%.17g,%.17g,%d\n", tag.Data(), s9.mean, s9.rms, s9.max, s9.n);
  }
}

static bool uploadOne(const TString& inFile,
                      const TString& ccdbPath,
                      const TString& ccdbUrl,
                      long startValidity,
                      long endValidity,
                      const Acceptance& acc,
                      bool doUpload)
{
  std::unique_ptr<TFile> f(TFile::Open(inFile, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] cannot open weight file " << inFile << "\n";
    return false;
  }
  auto* h0 = dynamic_cast<TH3D*>(f->Get("ccdb_object"));
  if (!h0) {
    std::cerr << "[ERROR] missing TH3D ccdb_object in " << inFile << "\n";
    return false;
  }
  std::unique_ptr<TH3D> h(dynamic_cast<TH3D*>(h0->Clone("ccdb_object")));
  if (!h) return false;
  h->SetDirectory(nullptr);

  std::cout << "Upload object: " << inFile << " -> " << ccdbPath << "\n";
  std::cout << "  bins = " << h->GetNbinsX() << " x " << h->GetNbinsY() << " x " << h->GetNbinsZ()
            << ", integral = " << std::setprecision(12) << integral3D(h.get()) << "\n";

  if (!doUpload) {
    std::cout << "  [DRY RUN] not uploaded\n";
    return true;
  }

  o2::ccdb::CcdbApi api;
  api.init(ccdbUrl.Data());

  std::map<std::string, std::string> meta;
  meta["Description"] = "Stage-1 shape-only residual epsilon map; task divides mixed-pair weight by epsilon";
  meta["ObjectName"] = "ccdb_object";
  meta["partName"] = "blob";
  meta["Convention"] = "TH3 axes: phi, y_or_eta, pT";
  meta["EpsilonDefinition"] = "(ME/intME)/(SE/intSE) = (REP/intREP)/(TGT/intTGT)";
  meta["ResidualWeight"] = "shape residual weight = 1/epsilon";
  meta["MeasurementAcceptance"] =
    Form("phi:[%.17g,%.17g),y_or_eta:[%.17g,%.17g),pT:[%.17g,%.17g); epsilon=1 outside",
         acc.phiMin, acc.phiMax, acc.yMin, acc.yMax, acc.ptMin, acc.ptMax);

  try {
    const int rc = api.storeAsTFileAny(h.get(), ccdbPath.Data(), meta, startValidity, endValidity);
    if (rc == 0) {
      std::cout << "  [OK] uploaded\n";
      return true;
    }
    std::cerr << "  [ERROR] storeAsTFileAny returned rc=" << rc << "\n";
  } catch (const std::exception& e) {
    std::cerr << "  [FATAL] upload failed: " << e.what() << "\n";
  }
  return false;
}

} // namespace weightqa


//Stage1Weights_shapeNormalized3D_ULcombined
void BuildUploadWeightAndQA_shapeNormalized(
  TString inputFile = "/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_woTPC.root",
  TString taskDirName = "lambdaspincorrderived_noTPCPID",
  //TString taskDirName = "lambdaspincorrderived",
  TString outDir = "/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/Stage1Weights_factorized_ULcombined",
  int legMode = 1,
  TString baseCcdbPath = "Users/y/yosu/Nov_5_noTPCPID",
  TString ccdbUrl = "https://alice-ccdb.cern.ch",
  bool doUpload =1,
  bool useRaw3D = 0,
  bool forceULFromLALALL = true,
  double epsMin = 0.001,
  double epsMax = 100.0,
  double minSEContent = 5.0,
  long startValidity = weightqa::kDefaultStartValidityMS,
  long endValidity = weightqa::kDefaultEndValidityMS,
  double acceptancePhiMin = 0.0,
  double acceptancePhiMax = 6.29,
  double acceptanceYMin = -0.5,
  double acceptanceYMax = 0.5,
  double acceptancePtMin = 0.6,
  double acceptancePtMax = 3.0)
{
  using namespace weightqa;
  const Acceptance acc{acceptancePhiMin, acceptancePhiMax,
                       acceptanceYMin, acceptanceYMax,
                       acceptancePtMin, acceptancePtMax};
  if (!(acc.phiMin < acc.phiMax) ||
      !(acc.yMin < acc.yMax) ||
      !(acc.ptMin < acc.ptMax)) {
    std::cerr << "[ERROR] invalid measurement-acceptance limits\n";
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);
  gSystem->mkdir(outDir + "/qa", kTRUE);

  std::unique_ptr<TFile> fin(TFile::Open(inputFile, "READ"));
  if (!fin || fin->IsZombie()) {
    std::cerr << "[ERROR] cannot open input file " << inputFile << "\n";
    return;
  }

  auto* taskDir = dynamic_cast<TDirectory*>(fin->Get(taskDirName));
  if (!taskDir) {
    std::cerr << "[ERROR] cannot find task directory " << taskDirName << "\n";
    return;
  }

  std::unique_ptr<TH3D> unityTemplate(makeTemplate(taskDir));
  if (!unityTemplate) {
    std::cerr << "[ERROR] cannot find any TGT/REP template histogram in " << taskDirName << "\n";
    return;
  }

  std::cout << "\nBuildUploadWeightAndQA_shapeNormalized\n"
            << "  inputFile:      " << inputFile << "\n"
            << "  taskDir:        " << taskDirName << "\n"
            << "  outDir:         " << outDir << "\n"
            << "  legMode:        " << legMode << "  (0=leg1, 1=leg2, 2=both)\n"
            << "  weight mode:    "
            << (useRaw3D ? "direct shape-normalized 3D"
                         : "shape decomposition: phi * 2D(y/eta,pT)") << "\n"
            << "  UL map:         " << (forceULFromLALALL ? "LAL+ALL combined" : "direct UL if present") << "\n"
            << "  epsilon:        (ME/intME)/(SE/intSE)\n"
            << "  acceptance:     phi=[" << acc.phiMin << "," << acc.phiMax
            << "), y/eta=[" << acc.yMin << "," << acc.yMax
            << "), pT=[" << acc.ptMin << "," << acc.ptMax << ")\n"
            << "  outside acc.:   epsilon = 1; excluded from normalization and QA\n"
            << "  fixed leg:      unity maps are written/uploaded for task compatibility\n"
            << "  task usage:     final mixed weight = mixpairweight / epsilon\n"
            << "  upload:         " << (doUpload ? "true" : "false") << "\n"
            << "  CCDB base path: " << baseCcdbPath << "\n\n";

  TString summaryName = outDir + "/weight_summary.csv";
  FILE* summary = fopen(summaryName.Data(), "w");
  if (summary) {
    fprintf(summary, "tag,status,integral_SE_TGT,integral_ME_REP,global_ME_over_SE,shapeScale,epsMinFinal,epsMaxFinal,nInsideBins,nOutsideBins,nDefault,nCapLow,nCapHigh,nZeroTGTWithREP,nZeroREPWithTGT\n");
  }

  TString qaCsvName = outDir + "/qa/qa_summary.csv";
  FILE* qacsv = fopen(qaCsvName.Data(), "w");
  if (qacsv) {
    fprintf(qacsv, "tag,state,projection,mean_abs_rminus1,rms_abs_rminus1,max_abs_rminus1,nBins\n");
  }

  int nBuilt = 0;
  int nUploaded = 0;
  int nRequested = 0;

  TString qaRootFile = outDir + "/qa/qa_histograms.root";
  TFile* qaFile = new TFile(qaRootFile, "RECREATE");



  for (const auto& sp : kSpecs) {
    const bool active = (legMode == 2) ||
                        (legMode == 0 && sp.leg == 1) ||
                        (legMode == 1 && sp.leg == 2);
    ++nRequested;

    // The analysis task retrieves all six CCDB paths during initialization,
    // even when only one leg is replaced. Supply neutral maps for the fixed
    // leg so the task can initialize without changing its retrieval logic.
    if (!active) {
      std::unique_ptr<TH3D> eps(cloneH3D(unityTemplate.get(), "ccdb_object"));
      if (!eps) {
        std::cerr << "[ERROR] cannot create inactive-leg unity map for "
                  << sp.tag << "\n";
        continue;
      }
      fillUnity(eps.get());
      eps->SetTitle(sp.tag + ": fixed/inactive leg; unity epsilon");

      std::unique_ptr<TH3D> residual(
        makeInverseWeight(eps.get(), "residual_weight_SE_over_ME", epsMin, epsMax));

      TString outFile = outDir + "/" + sp.tag + ".root";
      TFile fout(outFile, "RECREATE");
      eps->Write("ccdb_object");
      eps->Write("epsilon_ME_over_SE");
      if (residual) residual->Write("residual_weight_SE_over_ME");
      fout.Close();
      ++nBuilt;

      if (summary) {
        fprintf(summary,
                "%s,inactive_leg_unity,0,0,1,1,1,1,0,0,0,0,0,0,0\n",
                sp.tag.Data());
      }

      std::cout << "[OK] " << sp.tag
                << ": fixed/inactive leg; wrote unity CCDB map\n";

      TString ccdbPath = baseCcdbPath + "/" + sp.tag;
      if (uploadOne(outFile, ccdbPath, ccdbUrl,
                    startValidity, endValidity, acc, doUpload)) {
        if (doUpload) ++nUploaded;
      }
      continue;
    }

    std::unique_ptr<TH3D> hSE(getCat(taskDir, "TGT", sp.cat, sp.leg, forceULFromLALALL));
    std::unique_ptr<TH3D> hME(getCat(taskDir, "REP", sp.cat, sp.leg, forceULFromLALALL));
    if (!hSE || !hME) {
      std::cerr << "[ERROR] " << sp.tag
                << ": missing TGT or REP source. Skipping creation and upload.\n";
      if (summary) {
        fprintf(summary, "%s,missing_source,%.17g,%.17g,1,1,1,1,0,0,0,0,0,0,0\n",
                sp.tag.Data(), integral3D(hSE.get()), integral3D(hME.get()));
      }
      continue;
    }

    std::unique_ptr<TH3D> eps(buildEpsilonMEoverSE(hSE.get(), hME.get(), unityTemplate.get(), sp.tag,
                                                  useRaw3D, epsMin, epsMax, minSEContent,
                                                  acc, summary));
    if (!eps) continue;

    std::unique_ptr<TH3D> residual(makeInverseWeight(eps.get(), "residual_weight_SE_over_ME", epsMin, epsMax));

    TString outFile = outDir + "/" + sp.tag + ".root";
    TFile fout(outFile, "RECREATE");
    eps->Write("ccdb_object");
    eps->Write("epsilon_ME_over_SE");
    if (residual) residual->Write("residual_weight_SE_over_ME");
    if (hSE) hSE->Write("SE_TGT_source");
    if (hME) hME->Write("ME_REP_source");
    fout.Close();
    ++nBuilt;

    TString qaPng2D = outDir + "/qa/" + sp.tag + "_shapeQA_2D.png";
    TString qaPng1D = outDir + "/qa/" + sp.tag + "_shapeQA_1D.png";

    TString qaPdf2D = outDir + "/qa/" + sp.tag + "_shapeQA_2D.pdf";
    TString qaPdf1D = outDir + "/qa/" + sp.tag + "_shapeQA_1D.pdf";



    qaFile->cd();
    TDirectory* tagDir = qaFile->mkdir(sp.tag);

    makeQAPlot(qaPng2D, qaPng1D, sp.tag, hSE.get(), hME.get(), eps.get(), acc, qacsv, tagDir);
    makeQAPlot(qaPdf2D, qaPdf1D, sp.tag, hSE.get(), hME.get(), eps.get(), acc, qacsv, tagDir);

    TString ccdbPath = baseCcdbPath + "/" + sp.tag;
    if (uploadOne(outFile, ccdbPath, ccdbUrl, startValidity, endValidity, acc, doUpload)) {
      if (doUpload) ++nUploaded;
    }
  }

  qaFile->Close();
  delete qaFile;
  
  if (summary) fclose(summary);
  if (qacsv) fclose(qacsv);

  std::cout << "\nFinished.\n"
            << "  Weight files: " << outDir << "/*.root\n"
            << "  2D QA plots:  " << outDir << "/qa/*_shapeQA_2D.png\n"
            << "  1D QA plots:  " << outDir << "/qa/*_shapeQA_1D.png\n"
            << "  Summary:      " << summaryName << "\n"
            << "  QA summary:   " << qaCsvName << "\n"
            << "  Built:        " << nBuilt << " / " << nRequested << "\n";
  if (doUpload) std::cout << "  Uploaded:     " << nUploaded << " / " << nRequested << "\n";

  std::cout << "\nConfWeightPaths order:\n";
  std::cout << "\"ConfWeightPaths\": { \"values\": [\n";
  for (size_t i = 0; i < kSpecs.size(); ++i) {
    std::cout << "  \"" << baseCcdbPath << "/" << kSpecs[i].tag << "\"";
    if (i + 1 != kSpecs.size()) std::cout << ",";
    std::cout << "\n";
  }
  std::cout << "] },\n";
  std::cout << "\"useweight\": { \"value\": \"true\" },\n";
  std::cout << "\"useFixedWeight\": { \"value\": \"false\" },\n";
  std::cout << "\"cfgMixLegMode\": { \"value\": \"" << legMode << "\" }\n";
  std::cout << "\nWeight construction mode used here: "
            << (useRaw3D ? "direct normalized 3D"
                         : "normalized phi * 2D(y/eta,pT) decomposition")
            << "\n";

  if (!doUpload) {
    std::cout << "\nNOTE: doUpload=false, so this was a dry run. Pass doUpload=true to upload.\n";
  }
}
