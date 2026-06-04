// event_quality_yfp_ypfp_onefoil.C
//
// Event-level Yfp/YpFp polygon-selection bookkeeping for one run, one foil,
// and one DAT-defined delta slice.
//
// What this does:
//   - Auto-loads ROOT input file using the same candidate paths as the ytar macro.
//   - Parses DATfiles/list_of_optics_run.dat for OpticsID, foils, and delta slices.
//   - Loads the autonomous ytar multifoil cut:
//       cuts/ytar_delta_<run>_<ytarTag>_multifoil_cut.root
//   - Loads existing Yfp/YpFp hand polygon cuts using the repo naming convention:
//       cuts/YpFpYFp_<OpticsID>_-1_cut.root
//       cuts/YpFpYFp_<run>_-1_cut.root
//   - Applies the SAME base population used in ridge_width_central_guard_multifoil_cut.C:
//       H.cer.npeSum > 2.0
//       -10 < H.gtr.dp < 10
//       finite ytar, delta
//     plus finite FP/PFP and sieve variables needed for this stage.
//   - Keeps event-level records for future quality scoring / round-robin filling.
//   - Produces per-yscol Xs/Ys validation histograms so polygon numbering can be checked.
//
// Usage:
//   root -l
//   .x event_quality_yfp_ypfp_onefoil.C(1544,0,2,"multifoil_test1","fp_pfp_v1")
//
// Args:
//   nrun             run number
//   foilIndex        autonomous ytar foil index: delta_vs_ytar_cut_foil<foilIndex>
//   deltaSliceIndex  index into DAT delcut array: delcut[i] <= delta < delcut[i+1]
//   ytarTag          tag used by ridge_width_central_guard_multifoil_cut.C output
//   outTag           tag for this macro's output

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCutG.h"
#include "TKey.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

struct OpticsRunInfo {
  Int_t RunNum = 0;
  TString OpticsID = "";
  Double_t CentAngle = 0.0;
  Int_t NumFoil = 0;
  Int_t SieveFlag = 0;
  Int_t ndelcut = 0;
  Double_t ymis = 0.0;
  vector<Double_t> ztar_foil;
  vector<Double_t> delcut;
};

struct AxisRangeInfo {
  Double_t yfpMin = -40.0;
  Double_t yfpMax =  40.0;
  Double_t ypfpMin = -0.040;
  Double_t ypfpMax =  0.040;
  bool loaded = false;
};

struct RegionCut {
  Int_t yscol = -1;
  TCutG *cut = nullptr;
};

struct EventRecord {
  Long64_t entry = -1;

  Double_t sumnpe = 0.0;
  Double_t delta = 0.0;
  Double_t ytar = 0.0;

  Double_t yfp = 0.0;
  Double_t ypfp = 0.0;
  Double_t xfp = 0.0;
  Double_t xpfp = 0.0;

  Double_t xsieve = 0.0;
  Double_t ysieve = 0.0;

  Int_t pass_base = 0;
  Int_t pass_ytar_foil = 0;
  Int_t pass_delta_slice = 0;
  Int_t pass_yfp_ypfp_region = 0;

  Int_t region_index = -1;   // vector index in regions
  Int_t yscol = -1;          // physical Y-sieve row/column index from cut name

  Double_t density_score = 0.0;
  Double_t region_distance_score = 0.0;
  Double_t expected_ys_line_score = 0.0;
  Double_t nearest_ys_line_score = 0.0;
  Double_t total_quality_score = 0.0;

  Int_t nearest_ys_line = -999;
  Int_t nearest_xs_line = -999;
  Double_t dist_to_expected_ys = 999.0;
  Double_t dist_to_nearest_ys = 999.0;
  Double_t dist_to_nearest_xs = 999.0;
};

struct RegionSummary {
  Int_t region_index = -1;
  Int_t yscol = -1;
  Int_t n_events = 0;

  Double_t yfp_mean = 0.0;
  Double_t ypfp_mean = 0.0;
  Double_t xs_mean = 0.0;
  Double_t ys_mean = 0.0;

  Double_t yfp_rms = 0.0;
  Double_t ypfp_rms = 0.0;
  Double_t xs_rms = 0.0;
  Double_t ys_rms = 0.0;

  Int_t nearest_ys_line = -999;
  Int_t nearest_xs_line = -999;
  Double_t dist_to_expected_ys = 999.0;
  Double_t dist_to_nearest_ys = 999.0;
  Double_t dist_to_nearest_xs = 999.0;
};

bool HasBranch(TTree *T, const char *bname)
{
  return T && T->GetBranch(bname);
}

vector<TString> SplitCSVLine(const TString &line)
{
  vector<TString> fields;
  TObjArray *arr = line.Tokenize(",");
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    TString s = ((TObjString*)arr->At(i))->GetString();
    s = s.Strip(TString::kBoth);
    fields.push_back(s);
  }
  delete arr;
  return fields;
}

bool LoadOpticsRunInfo(Int_t nrun,
                       OpticsRunInfo &info,
                       TString opticsFile="DATfiles/list_of_optics_run.dat")
{
  ifstream file(opticsFile.Data());

  if (!file.is_open()) {
    cerr << "ERROR: could not open " << opticsFile << endl;
    return false;
  }

  TString line;
  bool found = false;

  while (line.ReadLine(file)) {
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0) continue;
    if (line.BeginsWith("RUN:")) continue;
    if (line.BeginsWith("RunNum")) continue;
    if (line.BeginsWith("ztar_foil")) continue;
    if (line.BeginsWith("delcut")) continue;
    if (!line.Contains(",")) continue;

    vector<TString> fields = SplitCSVLine(line);
    if (fields.size() < 7) continue;

    Int_t thisRun = fields[0].Atoi();
    if (thisRun != nrun) continue;

    info.RunNum = thisRun;
    info.OpticsID = fields[1];
    info.CentAngle = fields[2].Atof();
    info.NumFoil = fields[3].Atoi();
    info.SieveFlag = fields[4].Atoi();
    info.ndelcut = fields[5].Atoi();
    info.ymis = fields[6].Atof();

    TString foilLine;
    TString delLine;
    if (!foilLine.ReadLine(file)) break;
    if (!delLine.ReadLine(file)) break;

    vector<TString> foilFields = SplitCSVLine(foilLine);
    vector<TString> delFields  = SplitCSVLine(delLine);

    info.ztar_foil.clear();
    info.delcut.clear();

    for (auto &s : foilFields) {
      if (s.Length() > 0) info.ztar_foil.push_back(s.Atof());
    }

    for (auto &s : delFields) {
      if (s.Length() > 0) info.delcut.push_back(s.Atof());
    }

    found = true;
    break;
  }

  if (!found) {
    cerr << "ERROR: run " << nrun << " not found in " << opticsFile << endl;
    return false;
  }

  if ((Int_t)info.ztar_foil.size() != info.NumFoil) {
    cerr << "ERROR: parsed " << info.ztar_foil.size()
         << " foil positions but NumFoil=" << info.NumFoil << endl;
    return false;
  }

  if ((Int_t)info.delcut.size() != info.ndelcut + 1) {
    cerr << "ERROR: parsed " << info.delcut.size()
         << " delta boundaries but ndelcut+1=" << info.ndelcut+1 << endl;
    return false;
  }

  cout << "Loaded optics run info from " << opticsFile << endl;
  cout << "  RunNum   = " << info.RunNum << endl;
  cout << "  OpticsID = " << info.OpticsID << endl;
  cout << "  NumFoil  = " << info.NumFoil << endl;
  cout << "  ndelcut  = " << info.ndelcut << endl;

  for (Int_t i = 0; i < info.NumFoil; i++) {
    cout << "  foil " << i << " ztar = " << info.ztar_foil[i] << endl;
  }

  for (Int_t i = 0; i < info.ndelcut; i++) {
    cout << "  delta slice " << i << " = ["
         << info.delcut[i] << ", " << info.delcut[i+1]
         << "] center=" << 0.5*(info.delcut[i] + info.delcut[i+1]) << endl;
  }

  return true;
}

AxisRangeInfo LoadYfpYpfpAxisRange(Int_t deltaSliceIndex,
                                   TString filename="DATfiles/AxisRange_ypfp_yfp.dat")
{
  AxisRangeInfo out;

  ifstream file(filename.Data());
  if (!file.is_open()) {
    cerr << "WARNING: could not open " << filename
         << "; using default Yfp/YpFp ranges." << endl;
    return out;
  }

  TString line;
  Int_t idx = 0;
  while (line.ReadLine(file)) {
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0) continue;

    if (idx == deltaSliceIndex) {
      vector<TString> f = SplitCSVLine(line);
      if (f.size() >= 4) {
        out.yfpMin = f[0].Atof();
        out.yfpMax = f[1].Atof();
        out.ypfpMin = f[2].Atof();
        out.ypfpMax = f[3].Atof();
        out.loaded = true;
      }
      break;
    }
    idx++;
  }

  if (out.loaded) {
    cout << "Loaded Yfp/YpFp axis range for delta slice " << deltaSliceIndex << endl;
    cout << "  yfp  = [" << out.yfpMin  << ", " << out.yfpMax  << "]" << endl;
    cout << "  ypfp = [" << out.ypfpMin << ", " << out.ypfpMax << "]" << endl;
  } else {
    cerr << "WARNING: no axis range found for delta slice " << deltaSliceIndex
         << "; using defaults." << endl;
  }

  return out;
}

TString BuildInputRootPath(Int_t nrun, TString inputFileID="-1")
{
  vector<TString> candidates;

  candidates.push_back(Form("ROOTfiles/OPTICS/nps_hms_optics_%d_1_%s.root",
                            nrun, inputFileID.Data()));

  for (auto &p : candidates) {
    if (!gSystem->AccessPathName(p)) {
      cout << "Found input ROOT file: " << p << endl;
      return p;
    }
  }

  cerr << "ERROR: Could not auto-find input ROOT file. Tried:" << endl;
  for (auto &p : candidates) cerr << "  " << p << endl;

  return "";
}

TString BuildYtarCutPath(Int_t nrun, TString ytarTag)
{
  TString p = Form("cuts/ytar_delta_%d_%s_multifoil_cut.root",
                   nrun, ytarTag.Data());

  if (!gSystem->AccessPathName(p)) {
    cout << "Found ytar cut file: " << p << endl;
    return p;
  }

  cerr << "ERROR: ytar cut file not found: " << p << endl;
  return "";
}

TString BuildYfpYpfpCutPath(Int_t nrun, TString opticsID, TString handID="-1")
{
  vector<TString> candidates;

  candidates.push_back(Form("cuts/YpFpYFp_%s_%s_cut.root",
                            opticsID.Data(), handID.Data()));
  candidates.push_back(Form("cuts/YpFpYFp_%d_%s_cut.root",
                            nrun, handID.Data()));
  candidates.push_back(Form("cuts/YpFpYFp_%s_-1_cut.root",
                            opticsID.Data()));
  candidates.push_back(Form("cuts/YpFpYFp_%d_-1_cut.root",
                            nrun));

  for (auto &p : candidates) {
    if (!gSystem->AccessPathName(p)) {
      cout << "Found Yfp/YpFp cut file: " << p << endl;
      return p;
    }
  }

  cerr << "ERROR: could not find Yfp/YpFp cut file. Tried:" << endl;
  for (auto &p : candidates) cerr << "  " << p << endl;

  return "";
}

vector<RegionCut> LoadYfpYpfpRegions(TString filename,
                                     Int_t foilIndex,
                                     Int_t deltaSliceIndex)
{
  vector<RegionCut> regions;

  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cerr << "ERROR: could not open Yfp/YpFp cut file: " << filename << endl;
    return regions;
  }

  for (Int_t yscol = 0; yscol < 9; yscol++) {
    TString cname = Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",
                         yscol, foilIndex, deltaSliceIndex);

    TCutG *c = (TCutG*)f->Get(cname);
    if (!c) continue;

    TCutG *clone = (TCutG*)c->Clone(Form("%s_active", cname.Data()));
    clone->SetDirectory(nullptr);
    clone->SetLineColor(kRed);
    clone->SetLineWidth(2);

    RegionCut rc;
    rc.yscol = yscol;
    rc.cut = clone;
    regions.push_back(rc);
  }

  f->Close();

  cout << "Loaded " << regions.size()
       << " Yfp/YpFp yscol cuts for foil " << foilIndex
       << " delta slice " << deltaSliceIndex << endl;

  for (auto &r : regions) {
    cout << "  yscol " << r.yscol << " : " << r.cut->GetName() << endl;
  }

  return regions;
}

Int_t NearestLineIndex(Double_t val, const vector<Double_t> &lines, Double_t &dist)
{
  Int_t best = -1;
  dist = 999.0;

  for (size_t i = 0; i < lines.size(); i++) {
    Double_t d = std::abs(val - lines[i]);
    if (d < dist) {
      dist = d;
      best = (Int_t)i;
    }
  }

  return best;
}

Double_t Clamp01(Double_t x)
{
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

Double_t ExpScore(Double_t dist, Double_t sigma)
{
  if (sigma <= 0.0) return 0.0;
  return std::exp(-0.5 * std::pow(dist / sigma, 2));
}

void DrawSieveLines(const vector<Double_t> &xs_cent,
                    const vector<Double_t> &ys_cent,
                    Double_t ysMin,
                    Double_t ysMax,
                    Double_t xsMin,
                    Double_t xsMax,
                    Int_t highlightYs=-1)
{
  for (size_t i = 0; i < ys_cent.size(); i++) {
    TLine *l = new TLine(ys_cent[i], xsMin, ys_cent[i], xsMax);
    l->SetLineColor(kRed);
    l->SetLineWidth((Int_t)i == highlightYs ? 4 : 1);
    l->SetLineStyle((Int_t)i == highlightYs ? 1 : 2);
    l->Draw("same");

    TText *txt = new TText(ys_cent[i], xsMax * 0.92, Form("Y%d", (Int_t)i));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.025);
    txt->Draw("same");
  }

  for (size_t i = 0; i < xs_cent.size(); i++) {
    TLine *l = new TLine(ysMin, xs_cent[i], ysMax, xs_cent[i]);
    l->SetLineColor(kRed+1);
    l->SetLineWidth(1);
    l->SetLineStyle(3);
    l->Draw("same");
  }
}

Int_t FindRegionIndex(Double_t ypfp, Double_t yfp,
                      const vector<RegionCut> &regions)
{
  for (size_t i = 0; i < regions.size(); i++) {
    if (!regions[i].cut) continue;
    if (regions[i].cut->IsInside(ypfp, yfp)) return (Int_t)i;
  }
  return -1;
}

void event_quality_yfp_ypfp_onefoil(
    Int_t nrun = 1544,
    Int_t foilIndex = 0,
    Int_t deltaSliceIndex = 2,
    const char *ytarTag = "multifoil_test1",
    const char *outTag = "fp_pfp_test"
)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TString inputID = "-1";
  TString handID  = "-1";
  TString ytarID  = ytarTag;
  TString tagID   = outTag;

  OpticsRunInfo runInfo;
  if (!LoadOpticsRunInfo(nrun, runInfo)) return;

  if (foilIndex < 0 || foilIndex >= runInfo.NumFoil) {
    cerr << "ERROR: foilIndex=" << foilIndex
         << " outside [0," << runInfo.NumFoil-1 << "]" << endl;
    return;
  }

  if (deltaSliceIndex < 0 || deltaSliceIndex >= runInfo.ndelcut) {
    cerr << "ERROR: deltaSliceIndex=" << deltaSliceIndex
         << " outside [0," << runInfo.ndelcut-1 << "]" << endl;
    return;
  }

  Double_t delMin = runInfo.delcut[deltaSliceIndex];
  Double_t delMax = runInfo.delcut[deltaSliceIndex + 1];
  Double_t delCenter = 0.5 * (delMin + delMax);

  AxisRangeInfo axis = LoadYfpYpfpAxisRange(deltaSliceIndex);

  // Sieve line locations used by existing Yfp/YpFp and Xfp/XpFp cut macros.
  vector<Double_t> xs_cent;
  vector<Double_t> ys_cent;
  for (Int_t i = 0; i < 9; i++) {
    xs_cent.push_back((i - 4) * 2.54);
    ys_cent.push_back((i - 4) * 0.6 * 2.54);
  }

  TString inputroot = BuildInputRootPath(nrun, inputID);
  if (inputroot == "") return;

  TString ytarCutFile = BuildYtarCutPath(nrun, ytarID);
  if (ytarCutFile == "") return;

  TString yfpCutFile = BuildYfpYpfpCutPath(nrun, runInfo.OpticsID, handID);
  if (yfpCutFile == "") return;

  vector<RegionCut> regions = LoadYfpYpfpRegions(yfpCutFile, foilIndex, deltaSliceIndex);
  if (regions.empty()) {
    cerr << "ERROR: no Yfp/YpFp region cuts found for foil " << foilIndex
         << " delta slice " << deltaSliceIndex << endl;
    cerr << "Expected names like: hYpFpYFp_cut_yscol_<0-8>_nfoil_"
         << foilIndex << "_ndel_" << deltaSliceIndex << endl;
    return;
  }

  TString outRootName =
    Form("hist/event_quality_yfp_ypfp_run%d_foil%d_del%d_%s.root",
         nrun, foilIndex, deltaSliceIndex, tagID.Data());

  TString outPdfName =
    Form("plots/event_quality_yfp_ypfp_run%d_foil%d_del%d_%s.pdf",
         nrun, foilIndex, deltaSliceIndex, tagID.Data());

  TString outCsvName =
    Form("plots/event_quality_yfp_ypfp_run%d_foil%d_del%d_%s_region_summary.csv",
         nrun, foilIndex, deltaSliceIndex, tagID.Data());

  gSystem->mkdir("hist", kTRUE);
  gSystem->mkdir("plots", kTRUE);

  cout << "\n=== event_quality_yfp_ypfp_onefoil ===" << endl;
  cout << "Run:           " << nrun << endl;
  cout << "OpticsID:      " << runInfo.OpticsID << endl;
  cout << "Input ROOT:    " << inputroot << endl;
  cout << "Ytar cuts:     " << ytarCutFile << endl;
  cout << "Yfp/YpFp cuts: " << yfpCutFile << endl;
  cout << "Foil index:    " << foilIndex
       << " ztar=" << runInfo.ztar_foil[foilIndex] << endl;
  cout << "Delta slice:   " << deltaSliceIndex
       << " [" << delMin << ", " << delMax << "]"
       << " center=" << delCenter << endl;
  cout << "Output ROOT:   " << outRootName << endl;
  cout << "Output PDF:    " << outPdfName << endl;
  cout << "Output CSV:    " << outCsvName << endl;

  TFile *fin = TFile::Open(inputroot, "READ");
  if (!fin || fin->IsZombie()) {
    cerr << "ERROR: could not open input ROOT file: " << inputroot << endl;
    return;
  }

  TTree *T = (TTree*)fin->Get("T");
  if (!T) {
    cerr << "ERROR: could not find tree T in " << inputroot << endl;
    fin->Close();
    return;
  }

  vector<const char*> requiredBranches = {
    "H.cer.npeSum",
    "H.gtr.dp",
    "H.gtr.y",
    "H.dc.y_fp",
    "H.dc.yp_fp",
    "H.dc.x_fp",
    "H.dc.xp_fp",
    "H.extcor.xsieve",
    "H.extcor.ysieve"
  };

  bool missingBranch = false;
  for (auto b : requiredBranches) {
    if (!HasBranch(T, b)) {
      cerr << "ERROR: missing required branch: " << b << endl;
      missingBranch = true;
    }
  }

  if (missingBranch) {
    fin->Close();
    return;
  }

  Double_t sumnpe = 0.0;
  Double_t delta = 0.0;
  Double_t ytar = 0.0;
  Double_t yfp = 0.0;
  Double_t ypfp = 0.0;
  Double_t xfp = 0.0;
  Double_t xpfp = 0.0;
  Double_t xsieve = 0.0;
  Double_t ysieve = 0.0;

  T->SetBranchStatus("*", 0);
  for (auto b : requiredBranches) T->SetBranchStatus(b, 1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchAddress("H.dc.x_fp", &xfp);
  T->SetBranchAddress("H.dc.xp_fp", &xpfp);
  T->SetBranchAddress("H.extcor.xsieve", &xsieve);
  T->SetBranchAddress("H.extcor.ysieve", &ysieve);

  TFile *fytar = TFile::Open(ytarCutFile, "READ");
  if (!fytar || fytar->IsZombie()) {
    cerr << "ERROR: could not open ytar cut file: " << ytarCutFile << endl;
    fin->Close();
    return;
  }

  TString ytarCutName = Form("delta_vs_ytar_cut_foil%d", foilIndex);
  TCutG *ytarCutRaw = (TCutG*)fytar->Get(ytarCutName);
  if (!ytarCutRaw) {
    cerr << "ERROR: could not find " << ytarCutName
         << " in " << ytarCutFile << endl;
    fytar->Close();
    fin->Close();
    return;
  }

  TCutG *ytarCut = (TCutG*)ytarCutRaw->Clone("ytarCut_active");
  ytarCut->SetDirectory(nullptr);
  fytar->Close();

  const Int_t nRegions = (Int_t)regions.size();

  TH2D *hYfpYpfp_all = new TH2D(
    "hYfpYpfp_all",
    Form("Run %d foil %d delta slice %d;YpFp;Yfp",
         nrun, foilIndex, deltaSliceIndex),
    240, axis.ypfpMin, axis.ypfpMax,
    240, axis.yfpMin, axis.yfpMax
  );

  TH2D *hYfpYpfp_selected = new TH2D(
    "hYfpYpfp_selected",
    Form("Selected Yfp/YpFp Run %d foil %d delta slice %d;YpFp;Yfp",
         nrun, foilIndex, deltaSliceIndex),
    240, axis.ypfpMin, axis.ypfpMax,
    240, axis.yfpMin, axis.yfpMax
  );

  TH2D *hXsYs_selected = new TH2D(
    "hXsYs_selected",
    Form("Combined selected Xs/Ys Run %d foil %d delta slice %d;Ysieve;Xsieve",
         nrun, foilIndex, deltaSliceIndex),
    240, -8.5, 8.5,
    240, -13.0, 13.0
  );

  vector<TH2D*> hYfpYpfp_region(nRegions, nullptr);
  vector<TH2D*> hXsYs_region(nRegions, nullptr);
  vector<TH1D*> hQuality_region(nRegions, nullptr);

  for (Int_t ir = 0; ir < nRegions; ir++) {
    Int_t yscol = regions[ir].yscol;

    hYfpYpfp_region[ir] = new TH2D(
      Form("hYfpYpfp_region_yscol%d", yscol),
      Form("Yfp/YpFp selected yscol %d Run %d foil %d del %d;YpFp;Yfp",
           yscol, nrun, foilIndex, deltaSliceIndex),
      240, axis.ypfpMin, axis.ypfpMax,
      240, axis.yfpMin, axis.yfpMax
    );

    hXsYs_region[ir] = new TH2D(
      Form("hXsYs_region_yscol%d", yscol),
      Form("Xs/Ys validation yscol %d Run %d foil %d del %d;Ysieve;Xsieve",
           yscol, nrun, foilIndex, deltaSliceIndex),
      240, -8.5, 8.5,
      240, -13.0, 13.0
    );

    hQuality_region[ir] = new TH1D(
      Form("hQuality_region_yscol%d", yscol),
      Form("Quality score yscol %d;quality;events", yscol),
      100, 0.0, 1.0
    );
  }

  // First pass: select events, fill density/all histograms, store event records.
  vector<EventRecord> records;

  Long64_t nentries = T->GetEntries();
  Long64_t nPassBase = 0;
  Long64_t nPassYtar = 0;
  Long64_t nPassDelta = 0;
  Long64_t nPassYfpRegion = 0;

  for (Long64_t iev = 0; iev < nentries; iev++) {
    T->GetEntry(iev);

    bool passBase =
      (sumnpe > 2.0) &&
      (delta > -10.0 && delta < 10.0) &&
      std::isfinite(delta) &&
      std::isfinite(ytar) &&
      std::isfinite(yfp) &&
      std::isfinite(ypfp) &&
      std::isfinite(xfp) &&
      std::isfinite(xpfp) &&
      std::isfinite(xsieve) &&
      std::isfinite(ysieve);

    if (!passBase) continue;
    nPassBase++;

    bool passYtarFoil = ytarCut->IsInside(ytar, delta);
    if (!passYtarFoil) continue;
    nPassYtar++;

    bool passDeltaSlice = (delta >= delMin && delta < delMax);
    if (!passDeltaSlice) continue;
    nPassDelta++;

    hYfpYpfp_all->Fill(ypfp, yfp);

    Int_t regionIndex = FindRegionIndex(ypfp, yfp, regions);
    if (regionIndex < 0) continue;
    nPassYfpRegion++;

    EventRecord rec;
    rec.entry = iev;
    rec.sumnpe = sumnpe;
    rec.delta = delta;
    rec.ytar = ytar;
    rec.yfp = yfp;
    rec.ypfp = ypfp;
    rec.xfp = xfp;
    rec.xpfp = xpfp;
    rec.xsieve = xsieve;
    rec.ysieve = ysieve;
    rec.pass_base = 1;
    rec.pass_ytar_foil = 1;
    rec.pass_delta_slice = 1;
    rec.pass_yfp_ypfp_region = 1;
    rec.region_index = regionIndex;
    rec.yscol = regions[regionIndex].yscol;

    records.push_back(rec);

    hYfpYpfp_selected->Fill(ypfp, yfp);
    hYfpYpfp_region[regionIndex]->Fill(ypfp, yfp);
    hXsYs_selected->Fill(ysieve, xsieve);
    hXsYs_region[regionIndex]->Fill(ysieve, xsieve);
  }

  cout << "\nEvent counts:" << endl;
  cout << "  pass base:        " << nPassBase << endl;
  cout << "  pass ytar foil:   " << nPassYtar << endl;
  cout << "  pass delta slice: " << nPassDelta << endl;
  cout << "  pass Yfp region:  " << nPassYfpRegion << endl;
  cout << "  stored records:   " << records.size() << endl;

  if (records.empty()) {
    cerr << "ERROR: no selected records. Stop." << endl;
    fin->Close();
    return;
  }

  // Smoothed density map. This is intentionally isolated so it can be swapped
  // for true KDE later without changing the event bookkeeping/output structure.
  TH2D *hDensity = (TH2D*)hYfpYpfp_all->Clone("hYfpYpfp_density_smooth");
  hDensity->SetTitle(Form("Smoothed Yfp/YpFp density Run %d foil %d del %d;YpFp;Yfp",
                          nrun, foilIndex, deltaSliceIndex));
  hDensity->Smooth(1);
  hDensity->Smooth(1);
  Double_t densityMax = hDensity->GetMaximum();
  if (densityMax <= 0.0) densityMax = 1.0;

  // Region summaries.
  vector<RegionSummary> summaries(nRegions);
  for (Int_t ir = 0; ir < nRegions; ir++) {
    summaries[ir].region_index = ir;
    summaries[ir].yscol = regions[ir].yscol;
  }

  vector<Double_t> sum_yfp(nRegions,0), sum_ypfp(nRegions,0), sum_xs(nRegions,0), sum_ys(nRegions,0);
  vector<Double_t> sum2_yfp(nRegions,0), sum2_ypfp(nRegions,0), sum2_xs(nRegions,0), sum2_ys(nRegions,0);

  for (auto &rec : records) {
    Int_t ir = rec.region_index;
    summaries[ir].n_events++;
    sum_yfp[ir]  += rec.yfp;
    sum_ypfp[ir] += rec.ypfp;
    sum_xs[ir]   += rec.xsieve;
    sum_ys[ir]   += rec.ysieve;

    sum2_yfp[ir]  += rec.yfp * rec.yfp;
    sum2_ypfp[ir] += rec.ypfp * rec.ypfp;
    sum2_xs[ir]   += rec.xsieve * rec.xsieve;
    sum2_ys[ir]   += rec.ysieve * rec.ysieve;
  }

  for (Int_t ir = 0; ir < nRegions; ir++) {
    Int_t n = summaries[ir].n_events;
    if (n <= 0) continue;

    summaries[ir].yfp_mean  = sum_yfp[ir] / n;
    summaries[ir].ypfp_mean = sum_ypfp[ir] / n;
    summaries[ir].xs_mean   = sum_xs[ir] / n;
    summaries[ir].ys_mean   = sum_ys[ir] / n;

    summaries[ir].yfp_rms  = std::sqrt(std::max(0.0, sum2_yfp[ir]/n  - summaries[ir].yfp_mean*summaries[ir].yfp_mean));
    summaries[ir].ypfp_rms = std::sqrt(std::max(0.0, sum2_ypfp[ir]/n - summaries[ir].ypfp_mean*summaries[ir].ypfp_mean));
    summaries[ir].xs_rms   = std::sqrt(std::max(0.0, sum2_xs[ir]/n   - summaries[ir].xs_mean*summaries[ir].xs_mean));
    summaries[ir].ys_rms   = std::sqrt(std::max(0.0, sum2_ys[ir]/n   - summaries[ir].ys_mean*summaries[ir].ys_mean));

    summaries[ir].nearest_ys_line = NearestLineIndex(summaries[ir].ys_mean, ys_cent,
                                                     summaries[ir].dist_to_nearest_ys);
    summaries[ir].nearest_xs_line = NearestLineIndex(summaries[ir].xs_mean, xs_cent,
                                                     summaries[ir].dist_to_nearest_xs);

    if (summaries[ir].yscol >= 0 && summaries[ir].yscol < (Int_t)ys_cent.size()) {
      summaries[ir].dist_to_expected_ys = std::abs(summaries[ir].ys_mean - ys_cent[summaries[ir].yscol]);
    }
  }

  // Score each selected event.
  const Double_t sigmaExpectedYs = 0.65; // cm-ish, loose first-pass diagnostic
  const Double_t sigmaNearestYs  = 0.65;
  const Double_t sigmaNearestXs  = 1.20;

  for (auto &rec : records) {
    Int_t bx = hDensity->GetXaxis()->FindBin(rec.ypfp);
    Int_t by = hDensity->GetYaxis()->FindBin(rec.yfp);
    rec.density_score = Clamp01(hDensity->GetBinContent(bx, by) / densityMax);

    RegionSummary &rs = summaries[rec.region_index];
    Double_t zy = 0.0;
    Double_t zp = 0.0;
    if (rs.yfp_rms > 0.0)  zy = (rec.yfp  - rs.yfp_mean)  / rs.yfp_rms;
    if (rs.ypfp_rms > 0.0) zp = (rec.ypfp - rs.ypfp_mean) / rs.ypfp_rms;
    Double_t d2 = zy*zy + zp*zp;
    rec.region_distance_score = std::exp(-0.5 * d2);

    if (rec.yscol >= 0 && rec.yscol < (Int_t)ys_cent.size()) {
      rec.dist_to_expected_ys = std::abs(rec.ysieve - ys_cent[rec.yscol]);
      rec.expected_ys_line_score = ExpScore(rec.dist_to_expected_ys, sigmaExpectedYs);
    }

    rec.nearest_ys_line = NearestLineIndex(rec.ysieve, ys_cent, rec.dist_to_nearest_ys);
    rec.nearest_xs_line = NearestLineIndex(rec.xsieve, xs_cent, rec.dist_to_nearest_xs);
    rec.nearest_ys_line_score = ExpScore(rec.dist_to_nearest_ys, sigmaNearestYs);

    // First-pass score. Keep components in tree so the formula can change later.
    rec.total_quality_score =
      0.40 * rec.density_score +
      0.30 * rec.region_distance_score +
      0.30 * rec.expected_ys_line_score;

    rec.total_quality_score = Clamp01(rec.total_quality_score);

    hQuality_region[rec.region_index]->Fill(rec.total_quality_score);
  }

  // Output ROOT file and trees.
  TFile *fout = TFile::Open(outRootName, "RECREATE");
  if (!fout || fout->IsZombie()) {
    cerr << "ERROR: could not create output ROOT file: " << outRootName << endl;
    fin->Close();
    return;
  }

  TTree *Tout = new TTree("event_quality", "event-level Yfp/YpFp selection and scoring");

  Int_t out_run = nrun;
  TString out_opticsID = runInfo.OpticsID;
  Int_t out_foilIndex = foilIndex;
  Double_t out_foilZ = runInfo.ztar_foil[foilIndex];
  Int_t out_deltaSliceIndex = deltaSliceIndex;
  Double_t out_delMin = delMin;
  Double_t out_delMax = delMax;
  Double_t out_delCenter = delCenter;

  EventRecord out;

  Tout->Branch("run", &out_run, "run/I");
  Tout->Branch("foilIndex", &out_foilIndex, "foilIndex/I");
  Tout->Branch("foilZ", &out_foilZ, "foilZ/D");
  Tout->Branch("deltaSliceIndex", &out_deltaSliceIndex, "deltaSliceIndex/I");
  Tout->Branch("delMin", &out_delMin, "delMin/D");
  Tout->Branch("delMax", &out_delMax, "delMax/D");
  Tout->Branch("delCenter", &out_delCenter, "delCenter/D");

  Tout->Branch("entry", &out.entry, "entry/L");
  Tout->Branch("sumnpe", &out.sumnpe, "sumnpe/D");
  Tout->Branch("delta", &out.delta, "delta/D");
  Tout->Branch("ytar", &out.ytar, "ytar/D");
  Tout->Branch("yfp", &out.yfp, "yfp/D");
  Tout->Branch("ypfp", &out.ypfp, "ypfp/D");
  Tout->Branch("xfp", &out.xfp, "xfp/D");
  Tout->Branch("xpfp", &out.xpfp, "xpfp/D");
  Tout->Branch("xsieve", &out.xsieve, "xsieve/D");
  Tout->Branch("ysieve", &out.ysieve, "ysieve/D");

  Tout->Branch("pass_base", &out.pass_base, "pass_base/I");
  Tout->Branch("pass_ytar_foil", &out.pass_ytar_foil, "pass_ytar_foil/I");
  Tout->Branch("pass_delta_slice", &out.pass_delta_slice, "pass_delta_slice/I");
  Tout->Branch("pass_yfp_ypfp_region", &out.pass_yfp_ypfp_region, "pass_yfp_ypfp_region/I");
  Tout->Branch("region_index", &out.region_index, "region_index/I");
  Tout->Branch("yscol", &out.yscol, "yscol/I");

  Tout->Branch("density_score", &out.density_score, "density_score/D");
  Tout->Branch("region_distance_score", &out.region_distance_score, "region_distance_score/D");
  Tout->Branch("expected_ys_line_score", &out.expected_ys_line_score, "expected_ys_line_score/D");
  Tout->Branch("nearest_ys_line_score", &out.nearest_ys_line_score, "nearest_ys_line_score/D");
  Tout->Branch("total_quality_score", &out.total_quality_score, "total_quality_score/D");

  Tout->Branch("nearest_ys_line", &out.nearest_ys_line, "nearest_ys_line/I");
  Tout->Branch("nearest_xs_line", &out.nearest_xs_line, "nearest_xs_line/I");
  Tout->Branch("dist_to_expected_ys", &out.dist_to_expected_ys, "dist_to_expected_ys/D");
  Tout->Branch("dist_to_nearest_ys", &out.dist_to_nearest_ys, "dist_to_nearest_ys/D");
  Tout->Branch("dist_to_nearest_xs", &out.dist_to_nearest_xs, "dist_to_nearest_xs/D");

  for (auto &rec : records) {
    out = rec;
    Tout->Fill();
  }

  TTree *Treg = new TTree("region_summary", "Yfp/YpFp region summary and Xs/Ys validation");

  RegionSummary rsout;
  Treg->Branch("run", &out_run, "run/I");
  Treg->Branch("foilIndex", &out_foilIndex, "foilIndex/I");
  Treg->Branch("foilZ", &out_foilZ, "foilZ/D");
  Treg->Branch("deltaSliceIndex", &out_deltaSliceIndex, "deltaSliceIndex/I");
  Treg->Branch("delMin", &out_delMin, "delMin/D");
  Treg->Branch("delMax", &out_delMax, "delMax/D");
  Treg->Branch("region_index", &rsout.region_index, "region_index/I");
  Treg->Branch("yscol", &rsout.yscol, "yscol/I");
  Treg->Branch("n_events", &rsout.n_events, "n_events/I");
  Treg->Branch("yfp_mean", &rsout.yfp_mean, "yfp_mean/D");
  Treg->Branch("ypfp_mean", &rsout.ypfp_mean, "ypfp_mean/D");
  Treg->Branch("xs_mean", &rsout.xs_mean, "xs_mean/D");
  Treg->Branch("ys_mean", &rsout.ys_mean, "ys_mean/D");
  Treg->Branch("yfp_rms", &rsout.yfp_rms, "yfp_rms/D");
  Treg->Branch("ypfp_rms", &rsout.ypfp_rms, "ypfp_rms/D");
  Treg->Branch("xs_rms", &rsout.xs_rms, "xs_rms/D");
  Treg->Branch("ys_rms", &rsout.ys_rms, "ys_rms/D");
  Treg->Branch("nearest_ys_line", &rsout.nearest_ys_line, "nearest_ys_line/I");
  Treg->Branch("nearest_xs_line", &rsout.nearest_xs_line, "nearest_xs_line/I");
  Treg->Branch("dist_to_expected_ys", &rsout.dist_to_expected_ys, "dist_to_expected_ys/D");
  Treg->Branch("dist_to_nearest_ys", &rsout.dist_to_nearest_ys, "dist_to_nearest_ys/D");
  Treg->Branch("dist_to_nearest_xs", &rsout.dist_to_nearest_xs, "dist_to_nearest_xs/D");

  ofstream csv(outCsvName.Data());
  csv << "run,foilIndex,foilZ,deltaSliceIndex,delMin,delMax,region_index,yscol,n_events,"
      << "yfp_mean,ypfp_mean,xs_mean,ys_mean,yfp_rms,ypfp_rms,xs_rms,ys_rms,"
      << "nearest_ys_line,nearest_xs_line,dist_to_expected_ys,dist_to_nearest_ys,dist_to_nearest_xs\n";

  for (auto &rs : summaries) {
    rsout = rs;
    Treg->Fill();

    csv << nrun << ","
        << foilIndex << ","
        << out_foilZ << ","
        << deltaSliceIndex << ","
        << delMin << ","
        << delMax << ","
        << rs.region_index << ","
        << rs.yscol << ","
        << rs.n_events << ","
        << rs.yfp_mean << ","
        << rs.ypfp_mean << ","
        << rs.xs_mean << ","
        << rs.ys_mean << ","
        << rs.yfp_rms << ","
        << rs.ypfp_rms << ","
        << rs.xs_rms << ","
        << rs.ys_rms << ","
        << rs.nearest_ys_line << ","
        << rs.nearest_xs_line << ","
        << rs.dist_to_expected_ys << ","
        << rs.dist_to_nearest_ys << ","
        << rs.dist_to_nearest_xs << "\n";
  }
  csv.close();

  // Diagnostic PDF.
  TCanvas *c = new TCanvas("c_event_quality", "event quality", 1300, 950);
  c->Print(outPdfName + "[");

  c->Clear();
  c->Divide(2,2);

  c->cd(1);
  gPad->SetLogz(1);
  hYfpYpfp_all->Draw("colz");
  for (auto &r : regions) {
    if (r.cut) r.cut->Draw("L same");
  }
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.12, 0.93, Form("Run %d foil %d del %d: all + Yfp/YpFp cuts",
                                nrun, foilIndex, deltaSliceIndex));

  c->cd(2);
  gPad->SetLogz(1);
  hDensity->Draw("colz");
  tx.DrawLatex(0.12, 0.93, "Smoothed density map");

  c->cd(3);
  gPad->SetLogz(1);
  hYfpYpfp_selected->Draw("colz");
  for (auto &r : regions) {
    if (r.cut) r.cut->Draw("L same");
  }
  tx.DrawLatex(0.12, 0.93, "Selected Yfp/YpFp events");

  c->cd(4);
  gPad->SetLogz(1);
  hXsYs_selected->Draw("colz");
  DrawSieveLines(xs_cent, ys_cent, -8.5, 8.5, -13.0, 13.0, -1);
  tx.DrawLatex(0.12, 0.93, "Combined Xs/Ys validation");

  c->Print(outPdfName);

  for (Int_t ir = 0; ir < nRegions; ir++) {
    Int_t yscol = regions[ir].yscol;

    c->Clear();
    c->Divide(2,2);

    c->cd(1);
    gPad->SetLogz(1);
    hYfpYpfp_region[ir]->Draw("colz");
    regions[ir].cut->Draw("L same");
    tx.DrawLatex(0.12, 0.93, Form("Yfp/YpFp selected yscol %d", yscol));

    c->cd(2);
    gPad->SetLogz(1);
    hXsYs_region[ir]->Draw("colz");
    DrawSieveLines(xs_cent, ys_cent, -8.5, 8.5, -13.0, 13.0, yscol);
    tx.DrawLatex(0.12, 0.93, Form("Xs/Ys validation yscol %d", yscol));

    c->cd(3);
    hQuality_region[ir]->Draw();

    c->cd(4);
    TH1D *hFrame = new TH1D(Form("hFrame_region_%d", ir), "", 1, 0, 1);
    hFrame->SetMinimum(0);
    hFrame->SetMaximum(1);
    hFrame->Draw("axis");

    TLatex tinfo;
    tinfo.SetNDC();
    tinfo.SetTextSize(0.045);

    RegionSummary &rs = summaries[ir];
    tinfo.DrawLatex(0.10, 0.88, Form("Region index: %d", ir));
    tinfo.DrawLatex(0.10, 0.78, Form("yscol: %d", rs.yscol));
    tinfo.DrawLatex(0.10, 0.68, Form("N events: %d", rs.n_events));
    tinfo.DrawLatex(0.10, 0.58, Form("mean Ys: %.3f", rs.ys_mean));
    tinfo.DrawLatex(0.10, 0.48, Form("mean Xs: %.3f", rs.xs_mean));
    tinfo.DrawLatex(0.10, 0.38, Form("nearest Ys line: %d", rs.nearest_ys_line));
    tinfo.DrawLatex(0.10, 0.28, Form("nearest Xs line: %d", rs.nearest_xs_line));
    tinfo.DrawLatex(0.10, 0.18, Form("dist to expected Ys: %.3f", rs.dist_to_expected_ys));

    c->Print(outPdfName);
  }

  c->Print(outPdfName + "]");

  fout->cd();
  Tout->Write();
  Treg->Write();

  hYfpYpfp_all->Write();
  hYfpYpfp_selected->Write();
  hDensity->Write();
  hXsYs_selected->Write();

  for (Int_t ir = 0; ir < nRegions; ir++) {
    hYfpYpfp_region[ir]->Write();
    hXsYs_region[ir]->Write();
    hQuality_region[ir]->Write();
    if (regions[ir].cut) regions[ir].cut->Write(Form("active_yfp_ypfp_cut_yscol%d", regions[ir].yscol));
  }

  fout->Close();
  fin->Close();

  cout << "\nDone." << endl;
  cout << "Wrote ROOT: " << outRootName << endl;
  cout << "Wrote PDF:  " << outPdfName << endl;
  cout << "Wrote CSV:  " << outCsvName << endl;
}
