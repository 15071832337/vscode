//*task for extracting weights//////////
//prottay/////////////////////

void Myweight()
{
  TString name_same[8];
  TString name_mix[8];

  name_same[0]= "TGT_LL_leg1";
  name_same[1]= "TGT_LAL_leg1";
  name_same[2]= "TGT_ALL_leg1";
  name_same[3]= "TGT_ALAL_leg1";

  name_same[4]= "TGT_LL_leg2";
  name_same[5]= "TGT_LAL_leg2";
  name_same[6]= "TGT_ALL_leg2";
  name_same[7]= "TGT_ALAL_leg2";

  name_mix[0]= "REP_LL_leg1";
  name_mix[1]= "REP_LAL_leg1";
  name_mix[2]= "REP_ALL_leg1";
  name_mix[3]= "REP_ALAL_leg1";

  name_mix[4]= "REP_LL_leg2";
  name_mix[5]= "REP_LAL_leg2";
  name_mix[6]= "REP_ALL_leg2";
  name_mix[7]= "REP_ALAL_leg2";

  TFile* fin = new TFile("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/AnalysisResults_unweight_setting5.root");
  TString outDir = "weight_maps";
  gSystem->mkdir(outDir.Data(),kTRUE);
  const double minContent = 50.0; // or 3, tune as you like

  for (int i = 0; i < 8; i++) {

    TH3D* hse = (TH3D*)fin->Get(Form("lambdaspincorrderived/%s", name_same[i].Data()));
    TH3D* hme = (TH3D*)fin->Get(Form("lambdaspincorrderived/%s", name_mix[i].Data()));

    //hme->Project3D("xz")->Draw();

    if (!hse || !hme) {
      std::cout << "Missing hist for index " << i << "  " << name_same[i]
                << " / " << name_mix[i] << std::endl;
      continue;
    }

    TH3D* hR = (TH3D*)hme->Clone("ccdb_object");
    hR->Reset(); // important: start from empty

    TH1D* hRatio = new TH1D("hRatio", "ME/SE per bin", 100, 0.0, 2.0);
    TH1D* hWeight = new TH1D("hWeight", "Weights per bin", 100, 0.0, 2.0);

    const int nx = hR->GetNbinsX();
    const int ny = hR->GetNbinsY();
    const int nz = hR->GetNbinsZ();

    // ---- PASS 1: build global ME/SE scale from good bins only ----
    double sumN = 0.0;
    double sumD = 0.0;

    for (int ix = 1; ix <= nx; ++ix) {
      for (int iy = 1; iy <= ny; ++iy) {
        for (int iz = 1; iz <= nz; ++iz) {

          const double n = hme->GetBinContent(ix, iy, iz);
          const double d = hse->GetBinContent(ix, iy, iz);

          if (n >= minContent) {
            const double r = n / d;
            if (std::isfinite(r) && r > 0.0) {
              hRatio->Fill(r);
              sumN += n;
              sumD += d;
            }
          }
        }
      }
    }
   

    double mean = 1.0;
    if (sumD > 0.0) {
      mean = sumN / sumD; // global ME/SE from "good" bins
    }
    std::cout << "Hist " << name_mix[i] << "  mean(ME/SE) = " << mean << std::endl;

    // ---- PASS 2: build weight map, normalised, with low-stat bins set to 1 ----
    for (int ix = 1; ix <= nx; ++ix) {
      for (int iy = 1; iy <= ny; ++iy) {
        for (int iz = 1; iz <= nz; ++iz) {

          const double n = hme->GetBinContent(ix, iy, iz);
          const double d = hse->GetBinContent(ix, iy, iz);

          double w = 1.0; // default neutral weight

          if (n >= minContent && sumD > 0.0) {
            double r = n / d;
            if (std::isfinite(r) && r > 0.0) {
              w = r / mean; // shape-only correction
              hWeight->Fill(w);
            } else {
              w = 1.0;
              hWeight->Fill(w);
            }
          }
          // optional: clamp extremes to avoid crazy weights
          // const double wmin = 0.2, wmax = 5.0;
          // if (w < wmin) w = wmin;
          // if (w > wmax) w = wmax;

          hR->SetBinContent(ix, iy, iz, w);
          hR->SetBinError(ix, iy, iz, 0.0);
        }
      }
    }

    std::cout << "  -> weight range: " << hR->GetMinimum()
              << "  to  " << hR->GetMaximum() << std::endl;

    
    TFile* fout = new TFile(Form("/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/Sourav_code/%s/pp_%s.root", outDir.Data(), name_mix[i].Data()), "RECREATE");
    fout->cd();
    hR->Write();
   
    hRatio->Write("ratio_distribution"); // optional debug
    hWeight->Write("weight_distribution"); // optional debug
    fout->Close();

    delete hRatio;
    delete hWeight;
  }

  fin->Close();
}
