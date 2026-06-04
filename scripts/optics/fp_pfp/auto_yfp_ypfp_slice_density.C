// auto_yfp_ypfp_slice_density.C
//
// Slice-wise automatic YpFp-vs-YFp cut builder for HMS/NPS optics.
//
// Purpose:
//   Develop an unsupervised, density-based replacement for hand-drawn
//   YpFp/YFp sieve-row cuts, starting with one foil and the central delta
//   slice/pair defined in DATfiles/list_of_optics_run.dat.
//
// Method, per selected delta slice:
//   1. Read delta boundaries from list_of_optics_run.dat.
//   2. Apply base PID + auto ytar foil cut + selected delta slice.
//   3. Fill and smooth a TH2D in YpFp vs YFp.
//   4. Threshold the smoothed density and find 8-connected components.
//   5. Convert each accepted component into a rotated-rectangle TCutG.
//   6. Apply each candidate cut and label it by the strongest Ys guide-line
//      population in the Xs-vs-Ys diagnostic plane.
//   7. Write diagnostic plots: combined YpFp/YFp, combined Xs/Ys, and one
//      Xs/Ys plot per assigned yscol.
//
// Usage:
//   hcana -l
//   .x auto_yfp_ypfp_slice_density.C(1544,"y_density_test1","-1","multifoil_test1",0,-2)
//
// ndelMode:
//   -2 : central pair/slice nearest delta=0  [default for development]
//   -1 : all delta slices                    [later]
//  >=0 : one explicit delta-slice index
//
// Notes:
//   This macro writes a tagged development cut file by default:
//     cuts/YpFpYFp_<OpticsID>_<inputFileID>_<tag>_cut.root
//   It does not overwrite the legacy expert-cut file unless you copy/rename it.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"

using namespace std;

struct OpticsRunInfo {
  int run = -1;
  TString opticsID = "";
  double centAngleDeg = 0.0;
  int numFoil = 0;
  int sieveFlag = 1;
  int ndelcut = 0;
  double ymis = 0.0;
  vector<double> ztarFoil;
  vector<double> delcut;
};

struct FPEventY {
  double ypfp = 0.0;
  double yfp = 0.0;
  double xs = 0.0;
  double ys = 0.0;
  double delta = 0.0;
  double ytar = 0.0;
};

struct BinPoint {
  int ix = 0;
  int iy = 0;
  double x = 0.0;
  double y = 0.0;
  double w = 0.0;
};

struct ComponentY {
  int rawIndex = -1;
  int yscol = -1;
  double score = 0.0;
  double totalWeight = 0.0;
  int areaBins = 0;
  double mx = 0.0;
  double my = 0.0;
  double vx = 1.0;
  double vy = 0.0;
  double ux = 0.0;
  double uy = 1.0;
  double sLo = 0.0;
  double sHi = 0.0;
  double uLo = 0.0;
  double uHi = 0.0;
  Long64_t nSelected = 0;
  TCutG *cut = nullptr;
  TH2D *hYsXs = nullptr;
  vector<BinPoint> bins;
};

bool HasBranch(TTree *T, const char *bname) {
  return T && T->GetBranch(bname);
}

vector<TString> SplitCSV(const TString &line) {
  vector<TString> out;
  TObjArray *arr = line.Tokenize(",");
  for (int i = 0; i < arr->GetEntries(); i++) {
    TString tok = ((TObjString*)arr->At(i))->GetString();
    tok = tok.Strip(TString::kBoth);
    out.push_back(tok);
  }
  delete arr;
  return out;
}

bool LooksLikeRunLine(const TString &line) {
  vector<TString> v = SplitCSV(line);
  if (v.size() < 7) return false;
  return v[0].IsDigit();
}

bool ReadOpticsRunInfo(int nrun, OpticsRunInfo &info,
                       const char *fname="DATfiles/list_of_optics_run.dat") {
  ifstream in(fname);
  if (!in.is_open()) {
    cerr << "ERROR: cannot open " << fname << endl;
    return false;
  }

  string sline;
  while (getline(in, sline)) {
    TString line(sline.c_str());
    line = line.Strip(TString::kBoth);
    if (!LooksLikeRunLine(line)) continue;

    vector<TString> main = SplitCSV(line);
    int run = main[0].Atoi();

    string zlineStd, dlineStd;
    if (!getline(in, zlineStd)) break;
    if (!getline(in, dlineStd)) break;

    if (run != nrun) continue;

    info.run = run;
    info.opticsID = main[1];
    info.centAngleDeg = main[2].Atof();
    info.numFoil = main[3].Atoi();
    info.sieveFlag = main[4].Atoi();
    info.ndelcut = main[5].Atoi();
    info.ymis = main[6].Atof();

    vector<TString> ztok = SplitCSV(TString(zlineStd.c_str()));
    vector<TString> dtok = SplitCSV(TString(dlineStd.c_str()));

    for (auto &z : ztok) if (z.Length() > 0) info.ztarFoil.push_back(z.Atof());
    for (auto &d : dtok) if (d.Length() > 0) info.delcut.push_back(d.Atof());

    if ((int)info.ztarFoil.size() < info.numFoil) {
      cerr << "ERROR: run " << nrun << " has fewer ztar foil values than NumFoil." << endl;
      return false;
    }

    if ((int)info.delcut.size() < 2) {
      cerr << "ERROR: run " << nrun << " has fewer than two delta boundaries." << endl;
      return false;
    }

    int ndelcut_from_header = info.ndelcut;
    int ndelcut_from_boundaries = (int)info.delcut.size() - 1;

    if (ndelcut_from_header != ndelcut_from_boundaries) {
      cout << "WARNING: header ndelcut=" << ndelcut_from_header
           << " but delta boundary line has " << info.delcut.size()
           << " values, so using " << ndelcut_from_boundaries
           << " usable delta intervals." << endl;
      info.ndelcut = ndelcut_from_boundaries;
    }

    cout << "Loaded run metadata: run=" << info.run
         << " OpticsID=" << info.opticsID
         << " NumFoil=" << info.numFoil
         << " usable_ndelcut=" << info.ndelcut << endl;
    return true;
  }

  cerr << "ERROR: did not find run " << nrun << " in " << fname << endl;
  return false;
}

vector<int> SelectDeltaSlices(const OpticsRunInfo &info, int ndelMode) {
  vector<int> out;

  if (ndelMode >= 0) {
    if (ndelMode < info.ndelcut) out.push_back(ndelMode);
    return out;
  }

  if (ndelMode == -1) {
    for (int nd = 0; nd < info.ndelcut; nd++) out.push_back(nd);
    return out;
  }

  // ndelMode == -2: central slice or central pair.  If zero is exactly on a
  // boundary, use the two neighboring slices for early development.
  const double eps = 1e-9;
  for (int nd = 0; nd < info.ndelcut; nd++) {
    double lo = info.delcut[nd];
    double hi = info.delcut[nd+1];
    if (std::abs(hi) < eps) out.push_back(nd);
    if (std::abs(lo) < eps) out.push_back(nd);
  }

  if (!out.empty()) {
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
    return out;
  }

  int best = -1;
  double bestAbs = 1e99;
  for (int nd = 0; nd < info.ndelcut; nd++) {
    double lo = info.delcut[nd];
    double hi = info.delcut[nd+1];
    if (lo <= 0.0 && 0.0 < hi) {
      out.push_back(nd);
      return out;
    }
    double cen = 0.5 * (lo + hi);
    if (std::abs(cen) < bestAbs) {
      bestAbs = std::abs(cen);
      best = nd;
    }
  }
  if (best >= 0) out.push_back(best);
  return out;
}

TString BuildInputRootPath(Int_t nrun, TString opticsID, TString inputFileID) {
  vector<TString> candidates;

  candidates.push_back(Form("ROOTfiles/OPTICS/nps_hms_optics_%d_1_%s.root", nrun, inputFileID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/nps_hms_optics_%s_1_%s.root", opticsID.Data(), inputFileID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/nps_hms_optics_%d_1_-1.root", nrun));
  candidates.push_back(Form("ROOTfiles/OPTICS/nps_hms_optics_%s_1_-1.root", opticsID.Data()));

  candidates.push_back(Form("ROOTfiles/OPTICS/5_878GeV/nps_hms_optics_hadd_%d_1_%s.root", nrun, inputFileID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_%d_1_%s.root", nrun, inputFileID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/5_878GeV/nps_hms_optics_hadd_%s_1_%s.root", opticsID.Data(), inputFileID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_%s_1_%s.root", opticsID.Data(), inputFileID.Data()));

  candidates.push_back(Form("ROOTfiles/OPTICS/5_878GeV/nps_hms_optics_hadd_%d_1_-1.root", nrun));
  candidates.push_back(Form("ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_%d_1_-1.root", nrun));
  candidates.push_back(Form("ROOTfiles/OPTICS/5_878GeV/nps_hms_optics_hadd_%s_1_-1.root", opticsID.Data()));
  candidates.push_back(Form("ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_%s_1_-1.root", opticsID.Data()));

  for (auto &p : candidates) {
    if (!gSystem->AccessPathName(p)) {
      cout << "Found input ROOT file: " << p << endl;
      return p;
    }
  }

  cerr << "ERROR: Could not auto-find input ROOT file. Tried:" << endl;
  for (auto &p : candidates) cerr << "  " << p << endl;
  return candidates[0];
}

TCutG* LoadYtarCut(const OpticsRunInfo &info,
                   int foilIndex,
                   TString inputFileID,
                   TString ytarCutTag) {
  vector<TString> candidates;

  candidates.push_back(Form("cuts/ytar_delta_%d_%s_multifoil_cut.root",
                            info.run, ytarCutTag.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%d_%s_cut.root",
                            info.run, ytarCutTag.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%s_%s_cut.root",
                            info.opticsID.Data(), inputFileID.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%s_%d_cut.root",
                            info.opticsID.Data(), inputFileID.Atoi()));

  for (auto &path : candidates) {
    if (gSystem->AccessPathName(path)) continue;

    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) continue;

    TString cname = Form("delta_vs_ytar_cut_foil%d", foilIndex);
    TCutG *c = (TCutG*)f->Get(cname);
    if (c) {
      TCutG *clone = (TCutG*)c->Clone(Form("loaded_%s", cname.Data()));
      clone->SetLineColor(kRed);
      clone->SetLineWidth(3);
      cout << "Loaded ytar cut: " << path << " :: " << cname << endl;
      f->Close();
      return clone;
    }

    f->Close();
  }

  cerr << "ERROR: could not find ytar cut for foil " << foilIndex << endl;
  cerr << "Tried auto tag = " << ytarCutTag << " and legacy input ID = " << inputFileID << endl;
  return nullptr;
}

double Quantile(vector<double> v, double q) {
  if (v.empty()) return 0.0;
  sort(v.begin(), v.end());
  if (q <= 0.0) return v.front();
  if (q >= 1.0) return v.back();
  double pos = q * (v.size() - 1);
  int i = (int)floor(pos);
  double f = pos - i;
  if (i + 1 >= (int)v.size()) return v.back();
  return (1.0 - f) * v[i] + f * v[i+1];
}

void DrawYsGuideLines() {
  for (int nys = 0; nys < 9; nys++) {
    double pos = (nys - 4) * 0.6 * 2.54;
    TLine *line = new TLine(pos, -12.5, pos, 12.5);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw("same");

    TText *txt = new TText(pos, -12.0, Form("%d", nys));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.03);
    txt->SetTextAlign(22);
    txt->Draw("same");
  }
}

TCutG* MakeRotRectCut(const char *name,
                      double mx, double my,
                      double vx, double vy,
                      double ux, double uy,
                      double sLo, double sHi,
                      double uLo, double uHi) {
  TCutG *cut = new TCutG(name, 5);
  cut->SetTitle(Form("%s;YpFp;YFp", name));

  auto setp = [&](int i, double s, double u) {
    double x = mx + s * vx + u * ux;
    double y = my + s * vy + u * uy;
    cut->SetPoint(i, x, y);
  };

  setp(0, sLo, uLo);
  setp(1, sHi, uLo);
  setp(2, sHi, uHi);
  setp(3, sLo, uHi);
  setp(4, sLo, uLo);

  cut->SetLineColor(kRed);
  cut->SetLineWidth(3);
  return cut;
}

void CharacterizeComponent(ComponentY &c,
                           double sPadFrac,
                           double uPadFrac,
                           double sPadMin,
                           double uPadMin) {
  c.areaBins = (int)c.bins.size();
  c.totalWeight = 0.0;
  c.mx = 0.0;
  c.my = 0.0;

  for (auto &b : c.bins) {
    c.totalWeight += b.w;
    c.mx += b.w * b.x;
    c.my += b.w * b.y;
  }
  if (c.totalWeight <= 0.0) return;

  c.mx /= c.totalWeight;
  c.my /= c.totalWeight;

  double cxx = 0.0, cxy = 0.0, cyy = 0.0;
  for (auto &b : c.bins) {
    double dx = b.x - c.mx;
    double dy = b.y - c.my;
    cxx += b.w * dx * dx;
    cxy += b.w * dx * dy;
    cyy += b.w * dy * dy;
  }
  cxx /= c.totalWeight;
  cxy /= c.totalWeight;
  cyy /= c.totalWeight;

  double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
  c.vx = cos(theta);
  c.vy = sin(theta);
  if (c.vy < 0.0) { c.vx *= -1.0; c.vy *= -1.0; }
  c.ux = -c.vy;
  c.uy =  c.vx;

  vector<double> sVals, uVals;
  for (auto &b : c.bins) {
    double dx = b.x - c.mx;
    double dy = b.y - c.my;
    sVals.push_back(dx * c.vx + dy * c.vy);
    uVals.push_back(dx * c.ux + dy * c.uy);
  }

  double s0 = Quantile(sVals, 0.00);
  double s1 = Quantile(sVals, 1.00);
  double u0 = Quantile(uVals, 0.00);
  double u1 = Quantile(uVals, 1.00);

  double sPad = std::max(sPadMin, sPadFrac * (s1 - s0));
  double uPad = std::max(uPadMin, uPadFrac * (u1 - u0));

  c.sLo = s0 - sPad;
  c.sHi = s1 + sPad;
  c.uLo = u0 - uPad;
  c.uHi = u1 + uPad;
}

vector<ComponentY> FindConnectedComponents(TH2D *hSmooth,
                                           double relThreshold,
                                           double absThreshold,
                                           int minAreaBins,
                                           double minWeight) {
  vector<ComponentY> comps;
  if (!hSmooth) return comps;

  int nx = hSmooth->GetNbinsX();
  int ny = hSmooth->GetNbinsY();
  double hmax = hSmooth->GetMaximum();
  double threshold = std::max(absThreshold, relThreshold * hmax);

  vector<vector<int>> visited(nx + 1, vector<int>(ny + 1, 0));

  int dx[8] = {-1,0,1,-1,1,-1,0,1};
  int dy[8] = {-1,-1,-1,0,0,1,1,1};

  for (int ix = 1; ix <= nx; ix++) {
    for (int iy = 1; iy <= ny; iy++) {
      if (visited[ix][iy]) continue;
      double v = hSmooth->GetBinContent(ix, iy);
      if (v < threshold) continue;

      ComponentY comp;
      queue<pair<int,int>> q;
      q.push({ix, iy});
      visited[ix][iy] = 1;

      while (!q.empty()) {
        auto cur = q.front();
        q.pop();
        int cx = cur.first;
        int cy = cur.second;
        double cv = hSmooth->GetBinContent(cx, cy);
        if (cv < threshold) continue;

        BinPoint bp;
        bp.ix = cx;
        bp.iy = cy;
        bp.x = hSmooth->GetXaxis()->GetBinCenter(cx);
        bp.y = hSmooth->GetYaxis()->GetBinCenter(cy);
        bp.w = cv;
        comp.bins.push_back(bp);

        for (int k = 0; k < 8; k++) {
          int nx2 = cx + dx[k];
          int ny2 = cy + dy[k];
          if (nx2 < 1 || nx2 > nx || ny2 < 1 || ny2 > ny) continue;
          if (visited[nx2][ny2]) continue;
          if (hSmooth->GetBinContent(nx2, ny2) < threshold) continue;
          visited[nx2][ny2] = 1;
          q.push({nx2, ny2});
        }
      }

      CharacterizeComponent(comp, 0.35, 0.50, 0.0015, 1.0);
      if (comp.areaBins >= minAreaBins && comp.totalWeight >= minWeight) {
        comp.rawIndex = comps.size();
        comps.push_back(comp);
      }
    }
  }

  cout << "  density threshold = " << threshold
       << "  raw accepted components = " << comps.size() << endl;
  return comps;
}

void LabelComponentsByYs(vector<ComponentY> &components,
                         const vector<FPEventY> &events,
                         int nrun,
                         int foilIndex,
                         int ndel,
                         double ysWindow) {
  for (auto &c : components) {
    c.hYsXs = new TH2D(Form("hYsXs_auto_candidate_%d_nfoil_%d_ndel_%d", c.rawIndex, foilIndex, ndel),
                       Form("Run %d candidate %d: foil %d ndel %d;Ys;Xs", nrun, c.rawIndex, foilIndex, ndel),
                       100, -7.0, 7.0, 120, -12.5, 12.5);

    vector<int> score(9, 0);
    Long64_t nsel = 0;

    for (auto &e : events) {
      if (!c.cut) continue;
      if (!c.cut->IsInside(e.ypfp, e.yfp)) continue;
      c.hYsXs->Fill(e.ys, e.xs);
      nsel++;

      for (int iy = 0; iy < 9; iy++) {
        double ysNom = (iy - 4) * 0.6 * 2.54;
        if (std::abs(e.ys - ysNom) < ysWindow) score[iy]++;
      }
    }

    int best = -1;
    int bestScore = -1;
    for (int iy = 0; iy < 9; iy++) {
      if (score[iy] > bestScore) {
        bestScore = score[iy];
        best = iy;
      }
    }

    c.yscol = best;
    c.score = bestScore;
    c.nSelected = nsel;
  }
}

void auto_yfp_ypfp_slice_density(Int_t nrun=1544,
                                 const char *tag="y_density_test1",
                                 const char *inputFileID="-1",
                                 const char *ytarCutTag="multifoil_test1",
                                 Int_t foilIndex=0,
                                 Int_t ndelMode=-2)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TString outTag = tag;
  TString inputID = inputFileID;
  TString ytarTag = ytarCutTag;

  // Development parameters. These are intentionally visible and simple.
  const double npeMin = 2.0;
  const double baseDeltaMin = -10.0;
  const double baseDeltaMax =  10.0;

  const int nBinsYpFp = 180;
  const double ypfpMin = -0.035;
  const double ypfpMax =  0.035;
  const int nBinsYFp = 180;
  const double yfpMin = -35.0;
  const double yfpMax =  35.0;

  const int smoothPasses = 1;
  const double relThreshold = 0.16;
  const double absThreshold = 3.0;
  const int minAreaBins = 6;
  const double minComponentWeight = 25.0;
  const int maxComponentsToKeep = 12;
  const double ysWindow = 0.55;

  OpticsRunInfo info;
  if (!ReadOpticsRunInfo(nrun, info)) return;

  vector<int> deltaSlices = SelectDeltaSlices(info, ndelMode);
  if (deltaSlices.empty()) {
    cerr << "ERROR: no delta slices selected for ndelMode=" << ndelMode << endl;
    return;
  }

  cout << "Selected delta slices:";
  for (auto nd : deltaSlices) {
    cout << " ndel=" << nd << "[" << info.delcut[nd] << "," << info.delcut[nd+1] << "]";
  }
  cout << endl;

  TString inputroot = BuildInputRootPath(nrun, info.opticsID, inputID);
  TFile *fin = TFile::Open(inputroot, "READ");
  if (!fin || fin->IsZombie()) {
    cerr << "ERROR: could not open input ROOT file: " << inputroot << endl;
    return;
  }

  TTree *T = (TTree*)fin->Get("T");
  if (!T) {
    cerr << "ERROR: could not find tree T" << endl;
    return;
  }

  const char *required[] = {
    "H.cer.npeSum", "H.gtr.dp", "H.gtr.y",
    "H.dc.y_fp", "H.dc.yp_fp", "H.extcor.ysieve", "H.extcor.xsieve"
  };
  for (auto b : required) {
    if (!HasBranch(T, b)) {
      cerr << "ERROR: missing branch " << b << endl;
      return;
    }
  }

  TCutG *ytarCut = LoadYtarCut(info, foilIndex, inputID, ytarTag);
  if (!ytarCut) return;

  Double_t sumnpe = 0.0, delta = 0.0, ytar = 0.0;
  Double_t yfp = 0.0, ypfp = 0.0, ys = 0.0, xs = 0.0;

  cout << "Disabling unused branches." << endl;
  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("H.cer.npeSum", 1);
  T->SetBranchStatus("H.gtr.dp", 1);
  T->SetBranchStatus("H.gtr.y", 1);
  T->SetBranchStatus("H.dc.y_fp", 1);
  T->SetBranchStatus("H.dc.yp_fp", 1);
  T->SetBranchStatus("H.extcor.ysieve", 1);
  T->SetBranchStatus("H.extcor.xsieve", 1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchAddress("H.extcor.ysieve", &ys);
  T->SetBranchAddress("H.extcor.xsieve", &xs);

  TString outCutFile = Form("cuts/YpFpYFp_%s_%s_%s_cut.root", info.opticsID.Data(), inputID.Data(), outTag.Data());
  TString outPdf = Form("plots/auto_yfp_ypfp_slice_density_run%d_%s.pdf", nrun, outTag.Data());
  TString outCsv = Form("plots/auto_yfp_ypfp_slice_density_run%d_%s.csv", nrun, outTag.Data());

  gSystem->mkdir("cuts", kTRUE);
  gSystem->mkdir("plots", kTRUE);

  ofstream csv(outCsv.Data());
  csv << "run,tag,inputFileID,ytarCutTag,foilIndex,ndel,dlo,dhi,component,yscol,score,nSelected,areaBins,totalWeight,mx,my,vx,vy,sLo,sHi,uLo,uHi\n";

  TCanvas *c = new TCanvas("c_auto_yfp_density", "auto YpFp/YFp slice density", 1150, 850);
  c->Print(outPdf + "[");

  TFile *fout = new TFile(outCutFile, "RECREATE");

  Long64_t nentries = T->GetEntries();
  cout << "Tree entries = " << nentries << endl;

  for (auto ndel : deltaSlices) {
    double dLo = info.delcut[ndel];
    double dHi = info.delcut[ndel+1];
    double dC  = 0.5 * (dLo + dHi);

    cout << "\nProcessing ndel=" << ndel << " range=[" << dLo << ", " << dHi << "] center=" << dC << endl;

    vector<FPEventY> events;

    TH2D *hY = new TH2D(Form("hYpFpYFp_base_nfoil_%d_ndel_%d", foilIndex, ndel),
                        Form("Run %d foil %d ndel %d [%.1f, %.1f];YpFp;YFp", nrun, foilIndex, ndel, dLo, dHi),
                        nBinsYpFp, ypfpMin, ypfpMax,
                        nBinsYFp, yfpMin, yfpMax);

    TH2D *hYsXsBase = new TH2D(Form("hYsXs_base_nfoil_%d_ndel_%d", foilIndex, ndel),
                               Form("Run %d base events: foil %d ndel %d;Ys;Xs", nrun, foilIndex, ndel),
                               100, -7.0, 7.0, 120, -12.5, 12.5);

    for (Long64_t i = 0; i < nentries; i++) {
      T->GetEntry(i);

      if (!(sumnpe > npeMin)) continue;
      if (!(delta > baseDeltaMin && delta < baseDeltaMax)) continue;
      if (!(delta >= dLo && delta < dHi)) continue;
      if (!std::isfinite(ytar) || !std::isfinite(delta)) continue;
      if (!std::isfinite(yfp) || !std::isfinite(ypfp)) continue;
      if (!ytarCut->IsInside(ytar, delta)) continue;

      FPEventY ev;
      ev.ypfp = ypfp;
      ev.yfp = yfp;
      ev.xs = xs;
      ev.ys = ys;
      ev.delta = delta;
      ev.ytar = ytar;
      events.push_back(ev);
      hY->Fill(ypfp, yfp);
      hYsXsBase->Fill(ys, xs);
    }

    cout << "  selected events after ytar + delta = " << events.size() << endl;
    if (events.size() < 50) {
      cout << "  skipping ndel=" << ndel << " because too few events." << endl;
      continue;
    }

    TH2D *hSmooth = (TH2D*)hY->Clone(Form("hYpFpYFp_smooth_nfoil_%d_ndel_%d", foilIndex, ndel));
    for (int s = 0; s < smoothPasses; s++) hSmooth->Smooth(1);

    vector<ComponentY> comps = FindConnectedComponents(hSmooth, relThreshold, absThreshold,
                                                       minAreaBins, minComponentWeight);

    sort(comps.begin(), comps.end(), [](const ComponentY &a, const ComponentY &b) {
      return a.totalWeight > b.totalWeight;
    });
    if ((int)comps.size() > maxComponentsToKeep) comps.resize(maxComponentsToKeep);

    // Build tentative cuts.
    for (size_t ic = 0; ic < comps.size(); ic++) {
      comps[ic].rawIndex = ic;
      TString cnameTmp = Form("candidate_yfp_ypfp_component_%zu_nfoil_%d_ndel_%d", ic, foilIndex, ndel);
      comps[ic].cut = MakeRotRectCut(cnameTmp, comps[ic].mx, comps[ic].my,
                                     comps[ic].vx, comps[ic].vy,
                                     comps[ic].ux, comps[ic].uy,
                                     comps[ic].sLo, comps[ic].sHi,
                                     comps[ic].uLo, comps[ic].uHi);
    }

    LabelComponentsByYs(comps, events, nrun, foilIndex, ndel, ysWindow);

    // Resolve duplicate yscol assignments: keep the candidate with the highest score.
    map<int,int> bestCompForYs;
    for (int ic = 0; ic < (int)comps.size(); ic++) {
      int yscol = comps[ic].yscol;
      if (yscol < 0 || yscol > 8) continue;
      if (!bestCompForYs.count(yscol) || comps[ic].score > comps[bestCompForYs[yscol]].score) {
        bestCompForYs[yscol] = ic;
      }
    }

    vector<ComponentY> kept;
    for (auto &kv : bestCompForYs) kept.push_back(comps[kv.second]);
    sort(kept.begin(), kept.end(), [](const ComponentY &a, const ComponentY &b) {
      return a.yscol < b.yscol;
    });

    cout << "  kept assigned components = " << kept.size() << endl;
    for (auto &comp : kept) {
      TString cname = Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d", comp.yscol, foilIndex, ndel);
      comp.cut->SetName(cname);
      comp.cut->SetTitle(Form("%s;YpFp;YFp", cname.Data()));
      comp.cut->SetLineColor(kRed);
      comp.cut->SetLineWidth(3);
      cout << "    yscol=" << comp.yscol
           << " score=" << comp.score
           << " nSelected=" << comp.nSelected
           << " areaBins=" << comp.areaBins
           << " weight=" << comp.totalWeight << endl;

      csv << nrun << "," << outTag << "," << inputID << "," << ytarTag << ","
          << foilIndex << "," << ndel << "," << dLo << "," << dHi << ","
          << comp.rawIndex << "," << comp.yscol << "," << comp.score << ","
          << comp.nSelected << "," << comp.areaBins << "," << comp.totalWeight << ","
          << comp.mx << "," << comp.my << "," << comp.vx << "," << comp.vy << ","
          << comp.sLo << "," << comp.sHi << "," << comp.uLo << "," << comp.uHi << "\n";
    }

    // Combined selected Xs/Ys for this delta slice.
    TH2D *hYsXsCombined = new TH2D(Form("hYsXs_auto_combined_nfoil_%d_ndel_%d", foilIndex, ndel),
                                  Form("Run %d combined auto selections: foil %d ndel %d;Ys;Xs", nrun, foilIndex, ndel),
                                  100, -7.0, 7.0, 120, -12.5, 12.5);
    for (auto &e : events) {
      for (auto &comp : kept) {
        if (comp.cut && comp.cut->IsInside(e.ypfp, e.yfp)) {
          hYsXsCombined->Fill(e.ys, e.xs);
          break;
        }
      }
    }

    // PDF page 1: YpFp/YFp with auto components.
    c->Clear();
    gPad->SetLogz(1);
    hY->Draw("colz");
    for (auto &comp : kept) {
      if (comp.cut) comp.cut->Draw("L same");
      TText *txt = new TText(comp.mx, comp.my, Form("%d", comp.yscol));
      txt->SetTextColor(kRed);
      txt->SetTextSize(0.04);
      txt->SetTextAlign(22);
      txt->Draw("same");
    }
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.026);
    tx.DrawLatex(0.12, 0.94, Form("Run %d foil %d ndel %d [%.1f, %.1f]: YpFp/YFp density components", nrun, foilIndex, ndel, dLo, dHi));
    c->Print(outPdf);
    gPad->SetLogz(0);

    // PDF page 2: all ytar+delta events in Xs/Ys.
    c->Clear();
    gPad->SetLogz(1);
    hYsXsBase->Draw("colz");
    DrawYsGuideLines();
    tx.DrawLatex(0.12, 0.94, "All events after ytar + delta selection; red lines = nominal Ys guide positions");
    c->Print(outPdf);
    gPad->SetLogz(0);

    // PDF page 3: combined auto-selected events.
    c->Clear();
    gPad->SetLogz(1);
    hYsXsCombined->Draw("colz");
    DrawYsGuideLines();
    tx.DrawLatex(0.12, 0.94, "Combined auto YpFp/YFp selections for this delta slice");
    c->Print(outPdf);
    gPad->SetLogz(0);

    // PDF pages: one per assigned yscol.
    for (auto &comp : kept) {
      if (!comp.hYsXs) continue;
      c->Clear();
      gPad->SetLogz(1);
      comp.hYsXs->SetTitle(Form("Run %d foil %d ndel %d yscol %d only;Ys;Xs", nrun, foilIndex, ndel, comp.yscol));
      comp.hYsXs->Draw("colz");
      DrawYsGuideLines();
      tx.DrawLatex(0.12, 0.94, Form("Single auto YpFp/YFp selection: yscol %d, n=%lld, score=%.0f", comp.yscol, comp.nSelected, comp.score));
      c->Print(outPdf);
      gPad->SetLogz(0);
    }

    fout->cd();
    hY->Write("", TObject::kOverwrite);
    hSmooth->Write("", TObject::kOverwrite);
    hYsXsBase->Write("", TObject::kOverwrite);
    hYsXsCombined->Write("", TObject::kOverwrite);
    for (auto &comp : kept) {
      if (comp.cut) comp.cut->Write("", TObject::kOverwrite);
      if (comp.hYsXs) comp.hYsXs->Write("", TObject::kOverwrite);
    }
  }

  c->Print(outPdf + "]");
  fout->Close();
  csv.close();
  fin->Close();

  cout << "\nWrote tagged development cut file: " << outCutFile << endl;
  cout << "Wrote PDF: " << outPdf << endl;
  cout << "Wrote CSV: " << outCsv << endl;
  cout << "\nTo test with the ytar file you named, run:" << endl;
  cout << "  .x auto_yfp_ypfp_slice_density.C(" << nrun
       << ",\"" << outTag << "\",\"" << inputID << "\",\""
       << ytarTag << "\"," << foilIndex << ",-2)" << endl;
}
