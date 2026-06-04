// upload_fixweight_samevalid.C
// ROOT usage:
//   .L upload_fixweight_samevalid.C+
//   upload();                                  // uses default [1, 2524604400000]
//   upload(1, 2524604400000);                   // explicitly
//
// NOTE: The main fix vs your version is removing ".root" correctly
//       (Chop() removes only ONE character -> ".roo" bug).

#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH3.h"
#include "TSystem.h"
#include "TString.h"

#include "CCDB/CcdbApi.h"  // o2::ccdb::CcdbApi

// CCDB "upper bound" validity used by the service (≈ 2050-12-31 23:00 UTC)
static constexpr long kCCDBMaxValidityMS = 2524604400000L; // milliseconds since epoch
static constexpr long kEpochMS           = 1L;              // 1970-01-01 00:00:00

void upload(long startValidity = kEpochMS,
            long endValidity   = kCCDBMaxValidityMS)
{
  using namespace std;

  // Files you want to upload (already produced by your weighting macro)
  vector<string> filePaths = {
    "weight_maps/pp_REP_LL_leg1.root",
    "weight_maps/pp_REP_LAL_leg1.root",
    "weight_maps/pp_REP_ALL_leg1.root",
    "weight_maps/pp_REP_ALAL_leg1.root",
    "weight_maps/pp_REP_LL_leg2.root",
    "weight_maps/pp_REP_LAL_leg2.root",
    "weight_maps/pp_REP_ALL_leg2.root",
    "weight_maps/pp_REP_ALAL_leg2.root"
  };

  // Base CCDB path (no .root at the end)
  //
  //Users/y/yosu/test/cent0"
  const string baseCcdbPath = "Users/y/yosu/v6_May27_Data_fixedweight";

  o2::ccdb::CcdbApi api;
  api.init("https://alice-ccdb.cern.ch");

  cout << "Uploading " << filePaths.size()
       << " files with validity [" << startValidity << ", " << endValidity << "]\n";

  for (const auto& fp : filePaths) {

    // Key/tag = basename without ".root"
    TString tag = gSystem->BaseName(fp.c_str());

    // FIX: Chop() removes only ONE character. Remove ".root" properly.
    if (tag.EndsWith(".root")) {
      tag.Remove(tag.Length() - 5); // 5 chars in ".root"
      // alternatively: tag.ReplaceAll(".root", "");
    }

    cout << "\nOpening: " << fp << endl;

    std::unique_ptr<TFile> f(TFile::Open(fp.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      cerr << "  [ERROR] Cannot open " << fp << endl;
      continue;
    }

    // Fetch the object to upload
    auto* objIn = dynamic_cast<TH3D*>(f->Get("ccdb_object"));
    if (!objIn) {
      cerr << "  [ERROR] 'ccdb_object' (TH3D) not found in " << fp << endl;
      continue;
    }

    // Clone & detach from file
    auto* hOut = dynamic_cast<TH3D*>(objIn->Clone(Form("ccdb_object_%s", tag.Data())));
    if (!hOut) {
      cerr << "  [ERROR] Clone() failed for " << fp << endl;
      continue;
    }
    hOut->SetDirectory(nullptr);

    // Minimal metadata (feel free to add more)
    std::map<std::string, std::string> meta;
    meta["Description"] = Form("Fixed 3D weight (tag=%s)", tag.Data());
    meta["Author"]      = "Prottay Das";
    meta["partName"]    = "blob";

    // Final CCDB path is base path + '/' + tag (no .root suffix)
    const std::string ccdbPath = baseCcdbPath + "/" + tag.Data();

    cout << "  -> CCDB path:      " << ccdbPath << "\n"
         << "  -> Valid from:     " << startValidity << "\n"
         << "  -> Valid until:    " << endValidity   << "\n";

    // Upload as ROOT object with explicit validity window
    try {
      // storeAsTFileAny(object, path, metadata, start, end)
      int rc = api.storeAsTFileAny(hOut, ccdbPath, meta, startValidity, endValidity);

      if (rc == 0) {
        cout << "  [OK] Uploaded: " << tag << endl;
      } else {
        cerr << "  [ERROR] storeAsTFileAny() returned rc=" << rc
             << " for tag=" << tag << endl;
      }
    } catch (const std::exception& e) {
      cerr << "  [FATAL] Upload failed for tag=" << tag << " : " << e.what() << endl;
    }

    delete hOut;
  }

  cout << "\nDone.\n";
}
