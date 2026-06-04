// test_yfp_ypfp_one_island_localellipse.C
//
// Hyper-focused development macro for the Yfp/YpFp selection problem.
// Goal: in ONE delta slice, find ONE local high-density YpFp/YFp island
// and test whether it maps cleanly onto one Ys guide line in the Xs/Ys plane.
//
// This is not the production cut writer. It is an algorithm diagnostic.
// It writes a short PDF with:
//   1. all YpFp/YFp events with the best local candidate cut overlaid
//   2. all Xs/Ys events after ytar + delta selection
//   3. Xs/Ys events selected by the best YpFp/YFp candidate only
//   4. YpFp/YFp selected events only
//
// Usage:
//   hcana -l
//   .x test_yfp_ypfp_one_island_localellipse.C(1544,"oneIsland_test1","-1","multifoil_test1",0,2)
//
// Main knobs, if first try is too tight/loose:
//   minPeakFracOfMax : how high a local maximum must be relative to slice max
//   localRelLevel    : flood-fill threshold relative to that local peak
//   padScaleS/U      : geometric expansion of the resulting local rectangle

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <queue>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

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

struct LocalCandidate {
  int index = -1;
  int peakIx = -1;
  int peakIy = -1;
  double peakX = 0.0;
  double peakY = 0.0;
  double peakVal = 0.0;
  vector<BinPoint> bins;

  double weight = 0.0;
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

  // PCA/cuts are constructed in normalized coordinates so that
  // YpFp and YFp have comparable numerical scale.
  double xScale = 1.0;
  double yScale = 1.0;

  int yscol = -1;
  int bestScore = 0;
  int nSelected = 0;
  double purity = 0.0;
  double quality = 0.0;
  TCutG *cut = nullptr;
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

    if ((int)info.delcut.size() < 2) {
      cerr << "ERROR: run " << nrun << " has fewer than two delta boundaries." << endl;
      return false;
    }

    int ndel_from_boundaries = (int)info.delcut.size() - 1;
    if (info.ndelcut != ndel_from_boundaries) {
      cout << "WARNING: header ndelcut=" << info.ndelcut
           << " but boundary line gives " << ndel_from_boundaries
           << " usable intervals. Using boundaries." << endl;
      info.ndelcut = ndel_from_boundaries;
    }

    cout << "Loaded run metadata: run=" << info.run
         << " OpticsID=" << info.opticsID
         << " usable_ndelcut=" << info.ndelcut << endl;
    return true;
  }

  cerr << "ERROR: did not find run " << nrun << " in " << fname << endl;
  return false;
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
  candidates.push_back(Form("cuts/ytar_delta_%d_%s_multifoil_cut.root", info.run, ytarCutTag.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%d_%s_cut.root", info.run, ytarCutTag.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%s_%s_cut.root", info.opticsID.Data(), inputFileID.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%s_%d_cut.root", info.opticsID.Data(), inputFileID.Atoi()));

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

void CharacterizeCandidate(LocalCandidate &cand,
                           double padScaleS,
                           double padScaleU,
                           double padMinS,
                           double padMinU) {
  cand.weight = 0.0;
  cand.mx = 0.0;
  cand.my = 0.0;

  for (auto &b : cand.bins) {
    cand.weight += b.w;
    cand.mx += b.w * b.x;
    cand.my += b.w * b.y;
  }
  if (cand.weight <= 0.0) return;
  cand.mx /= cand.weight;
  cand.my /= cand.weight;

  double cxx = 0.0, cxy = 0.0, cyy = 0.0;
  for (auto &b : cand.bins) {
    double dx = b.x - cand.mx;
    double dy = b.y - cand.my;
    cxx += b.w * dx * dx;
    cxy += b.w * dx * dy;
    cyy += b.w * dy * dy;
  }
  cxx /= cand.weight;
  cxy /= cand.weight;
  cyy /= cand.weight;

  double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
  cand.vx = cos(theta);
  cand.vy = sin(theta);
  if (cand.vy < 0.0) { cand.vx *= -1.0; cand.vy *= -1.0; }
  cand.ux = -cand.vy;
  cand.uy =  cand.vx;

  vector<double> sVals, uVals;
  for (auto &b : cand.bins) {
    double dx = b.x - cand.mx;
    double dy = b.y - cand.my;
    sVals.push_back(dx * cand.vx + dy * cand.vy);
    uVals.push_back(dx * cand.ux + dy * cand.uy);
  }

  double s0 = Quantile(sVals, 0.00);
  double s1 = Quantile(sVals, 1.00);
  double u0 = Quantile(uVals, 0.00);
  double u1 = Quantile(uVals, 1.00);

  double sPad = std::max(padMinS, padScaleS * (s1 - s0));
  double uPad = std::max(padMinU, padScaleU * (u1 - u0));

  cand.sLo = s0 - sPad;
  cand.sHi = s1 + sPad;
  cand.uLo = u0 - uPad;
  cand.uHi = u1 + uPad;
}

TCutG* MakeCandidateCut(const char *name, const LocalCandidate &cand) {
  TCutG *cut = new TCutG(name, 5);
  cut->SetTitle(Form("%s;YpFp;YFp", name));

  auto setp = [&](int i, double s, double u) {
    double x = cand.mx + cand.xScale * (s * cand.vx + u * cand.ux);
    double y = cand.my + cand.yScale * (s * cand.vy + u * cand.uy);
    cut->SetPoint(i, x, y);
  };

  setp(0, cand.sLo, cand.uLo);
  setp(1, cand.sHi, cand.uLo);
  setp(2, cand.sHi, cand.uHi);
  setp(3, cand.sLo, cand.uHi);
  setp(4, cand.sLo, cand.uLo);

  cut->SetLineColor(kRed);
  cut->SetLineWidth(3);
  return cut;
}

vector<pair<int,int>> FindLocalMaxima(TH2D *h,
                                      double minPeakFracOfMax,
                                      int maxPeaks,
                                      double minSepX,
                                      double minSepY) {
  vector<pair<int,int>> peaks;
  vector<pair<double,pair<int,int>>> raw;
  if (!h) return peaks;

  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  double hmax = h->GetMaximum();
  double minVal = minPeakFracOfMax * hmax;

  for (int ix = 2; ix <= nx-1; ix++) {
    for (int iy = 2; iy <= ny-1; iy++) {
      double v = h->GetBinContent(ix, iy);
      if (v < minVal) continue;

      bool isMax = true;
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          if (dx == 0 && dy == 0) continue;
          if (h->GetBinContent(ix+dx, iy+dy) > v) isMax = false;
        }
      }
      if (!isMax) continue;
      raw.push_back({v, {ix, iy}});
    }
  }

  sort(raw.begin(), raw.end(), [](const auto &a, const auto &b) {
    return a.first > b.first;
  });

  for (auto &r : raw) {
    int ix = r.second.first;
    int iy = r.second.second;
    double x = h->GetXaxis()->GetBinCenter(ix);
    double y = h->GetYaxis()->GetBinCenter(iy);

    bool tooClose = false;
    for (auto &p : peaks) {
      double xp = h->GetXaxis()->GetBinCenter(p.first);
      double yp = h->GetYaxis()->GetBinCenter(p.second);
      double d2 = pow((x - xp)/minSepX, 2) + pow((y - yp)/minSepY, 2);
      if (d2 < 1.0) {
        tooClose = true;
        break;
      }
    }
    if (tooClose) continue;

    peaks.push_back({ix, iy});
    if ((int)peaks.size() >= maxPeaks) break;
  }

  return peaks;
}

LocalCandidate FloodFromPeak(TH2D *h,
                             int peakIx,
                             int peakIy,
                             int index,
                             double localRelLevel,
                             double absThreshold,
                             int maxBins) {
  LocalCandidate cand;
  cand.index = index;
  cand.peakIx = peakIx;
  cand.peakIy = peakIy;
  cand.peakX = h->GetXaxis()->GetBinCenter(peakIx);
  cand.peakY = h->GetYaxis()->GetBinCenter(peakIy);
  cand.peakVal = h->GetBinContent(peakIx, peakIy);

  double threshold = std::max(absThreshold, localRelLevel * cand.peakVal);
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();

  vector<vector<int>> visited(nx + 1, vector<int>(ny + 1, 0));
  queue<pair<int,int>> q;
  q.push({peakIx, peakIy});
  visited[peakIx][peakIy] = 1;

  int dx8[8] = {-1,0,1,-1,1,-1,0,1};
  int dy8[8] = {-1,-1,-1,0,0,1,1,1};

  while (!q.empty()) {
    auto cur = q.front();
    q.pop();
    int ix = cur.first;
    int iy = cur.second;
    double v = h->GetBinContent(ix, iy);
    if (v < threshold) continue;

    BinPoint bp;
    bp.ix = ix;
    bp.iy = iy;
    bp.x = h->GetXaxis()->GetBinCenter(ix);
    bp.y = h->GetYaxis()->GetBinCenter(iy);
    bp.w = v;
    cand.bins.push_back(bp);
    if ((int)cand.bins.size() >= maxBins) break;

    for (int k = 0; k < 8; k++) {
      int jx = ix + dx8[k];
      int jy = iy + dy8[k];
      if (jx < 1 || jx > nx || jy < 1 || jy > ny) continue;
      if (visited[jx][jy]) continue;
      if (h->GetBinContent(jx, jy) < threshold) continue;
      visited[jx][jy] = 1;
      q.push({jx, jy});
    }
  }

  return cand;
}



TCutG* MakeEllipseCut(const char *name, const LocalCandidate &cand, int npts=80) {
  TCutG *cut = new TCutG(name, npts + 1);
  cut->SetTitle(Form("%s;YpFp;YFp", name));

  for (int i = 0; i <= npts; i++) {
    double t = 2.0 * M_PI * ((double)i / (double)npts);
    double s = cand.sHi * cos(t);
    double u = cand.uHi * sin(t);
    double x = cand.mx + cand.xScale * (s * cand.vx + u * cand.ux);
    double y = cand.my + cand.yScale * (s * cand.vy + u * cand.uy);
    cut->SetPoint(i, x, y);
  }

  cut->SetLineColor(kRed);
  cut->SetLineWidth(3);
  return cut;
}

LocalCandidate BuildLocalEllipseCandidate(const vector<FPEventY> &events,
                                          TH2D *hSmooth,
                                          int peakIx,
                                          int peakIy,
                                          int index,
                                          double seedHalfYpFp,
                                          double seedHalfYFp,
                                          double nSigmaMajor,
                                          double nSigmaMinor,
                                          int minSeedEvents)
{
  LocalCandidate cand;
  cand.index = index;
  cand.peakIx = peakIx;
  cand.peakIy = peakIy;
  cand.peakX = hSmooth->GetXaxis()->GetBinCenter(peakIx);
  cand.peakY = hSmooth->GetYaxis()->GetBinCenter(peakIy);
  cand.peakVal = hSmooth->GetBinContent(peakIx, peakIy);
  cand.xScale = seedHalfYpFp;
  cand.yScale = seedHalfYFp;

  vector<const FPEventY*> seed;
  seed.reserve(2000);
  for (auto &e : events) {
    if (std::abs(e.ypfp - cand.peakX) > seedHalfYpFp) continue;
    if (std::abs(e.yfp  - cand.peakY) > seedHalfYFp) continue;
    seed.push_back(&e);
  }

  if ((int)seed.size() < minSeedEvents) {
    cand.quality = -1.0;
    return cand;
  }

  cand.weight = (double)seed.size();
  cand.mx = 0.0;
  cand.my = 0.0;
  for (auto *e : seed) {
    cand.mx += e->ypfp;
    cand.my += e->yfp;
  }
  cand.mx /= cand.weight;
  cand.my /= cand.weight;

  double cxx = 0.0, cxy = 0.0, cyy = 0.0;
  for (auto *e : seed) {
    double dx = (e->ypfp - cand.mx) / cand.xScale;
    double dy = (e->yfp  - cand.my) / cand.yScale;
    cxx += dx * dx;
    cxy += dx * dy;
    cyy += dy * dy;
  }
  cxx /= cand.weight;
  cxy /= cand.weight;
  cyy /= cand.weight;

  double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
  cand.vx = cos(theta);
  cand.vy = sin(theta);
  if (cand.vy < 0.0) { cand.vx *= -1.0; cand.vy *= -1.0; }
  cand.ux = -cand.vy;
  cand.uy =  cand.vx;

  vector<double> sVals, uVals;
  sVals.reserve(seed.size());
  uVals.reserve(seed.size());
  for (auto *e : seed) {
    double dx = (e->ypfp - cand.mx) / cand.xScale;
    double dy = (e->yfp  - cand.my) / cand.yScale;
    sVals.push_back(dx * cand.vx + dy * cand.vy);
    uVals.push_back(dx * cand.ux + dy * cand.uy);
  }

  // Robust scale from central quantiles, so a few tail events do not blow up the cut.
  double s16 = Quantile(sVals, 0.16);
  double s84 = Quantile(sVals, 0.84);
  double u16 = Quantile(uVals, 0.16);
  double u84 = Quantile(uVals, 0.84);
  double sigS = 0.5 * (s84 - s16);
  double sigU = 0.5 * (u84 - u16);

  // Hard guards in normalized space prevent the candidate from becoming a long band.
  double halfS = std::min(nSigmaMajor * sigS, 1.25);
  double halfU = std::min(nSigmaMinor * sigU, 0.85);

  // Fallbacks in normalized units.
  if (halfS < 0.20) halfS = 0.20;
  if (halfU < 0.12) halfU = 0.12;

  cand.sLo = -halfS;
  cand.sHi =  halfS;
  cand.uLo = -halfU;
  cand.uHi =  halfU;

  cand.cut = MakeEllipseCut(Form("local_ellipse_candidate_%d", index), cand);
  return cand;
}

void ScoreCandidate(LocalCandidate &cand,
                    const vector<FPEventY> &events,
                    double ysWindow,
                    int minSelected) {
  vector<int> score(9, 0);
  cand.nSelected = 0;

  for (auto &e : events) {
    if (!cand.cut) continue;
    if (!cand.cut->IsInside(e.ypfp, e.yfp)) continue;
    cand.nSelected++;

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

  cand.yscol = best;
  cand.bestScore = bestScore;
  cand.purity = (cand.nSelected > 0) ? ((double)cand.bestScore / (double)cand.nSelected) : 0.0;

  // Prefer a clean selection, but avoid choosing tiny statistical accidents.
  if (cand.nSelected < minSelected) cand.quality = -1.0;
  else cand.quality = cand.purity * sqrt((double)cand.bestScore);
}

void test_yfp_ypfp_one_island_localellipse(Int_t nrun=1544,
                              const char *tag="oneIsland_test1",
                              const char *inputFileID="-1",
                              const char *ytarCutTag="multifoil_test1",
                              Int_t foilIndex=0,
                              Int_t ndel=2,
                              Double_t minPeakFracOfMax=0.10,
                              Double_t localRelLevel=0.55,
                              Double_t seedHalfYpFp=0.0040,
                              Double_t seedHalfYFp=4.0,
                              Double_t nSigmaMajor=1.6,
                              Double_t nSigmaMinor=1.0)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  const double npeMin = 2.0;
  const int smoothPasses = 1;
  const int nBinsYpFp = 180;
  const double ypfpMin = -0.035;
  const double ypfpMax =  0.035;
  const int nBinsYFp = 180;
  const double yfpMin = -35.0;
  const double yfpMax =  35.0;

  const int maxPeaks = 20;
  const double minPeakSepYpFp = 0.0040;
  const double minPeakSepYFp  = 3.0;
  const double absFloodThreshold = 3.0;
  const int maxFloodBins = 500;
  const int minFloodBins = 4;
  const int minSelected = 100;
  const double ysWindow = 0.55;

  OpticsRunInfo info;
  if (!ReadOpticsRunInfo(nrun, info)) return;
  if (ndel < 0 || ndel >= info.ndelcut) {
    cerr << "ERROR: requested ndel=" << ndel << " outside valid range 0.." << info.ndelcut-1 << endl;
    return;
  }
  double dLo = info.delcut[ndel];
  double dHi = info.delcut[ndel+1];

  TString inputID = inputFileID;
  TString ytarTag = ytarCutTag;
  TString outTag = tag;

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

  TH2D *hYAll = new TH2D("hYpFpYFp_all",
                         Form("Run %d foil %d ndel %d [%.1f, %.1f];YpFp;YFp", nrun, foilIndex, ndel, dLo, dHi),
                         nBinsYpFp, ypfpMin, ypfpMax,
                         nBinsYFp, yfpMin, yfpMax);
  TH2D *hYsXsAll = new TH2D("hYsXs_all",
                            Form("Run %d all ytar+delta events;Ys;Xs", nrun),
                            100, -7.0, 7.0, 120, -12.5, 12.5);

  vector<FPEventY> events;
  Long64_t nentries = T->GetEntries();
  cout << "Tree entries = " << nentries << endl;
  cout << "Processing one delta slice: ndel=" << ndel << " [" << dLo << ", " << dHi << "]" << endl;

  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (!(sumnpe > npeMin)) continue;
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
    hYAll->Fill(ev.ypfp, ev.yfp);
    hYsXsAll->Fill(ev.ys, ev.xs);
  }

  cout << "Events after PID + ytar + delta = " << events.size() << endl;
  if (events.size() < 50) {
    cerr << "ERROR: too few events in this slice." << endl;
    return;
  }

  TH2D *hSmooth = (TH2D*)hYAll->Clone("hYpFpYFp_smooth");
  for (int i = 0; i < smoothPasses; i++) hSmooth->Smooth(1);

  vector<pair<int,int>> peaks = FindLocalMaxima(hSmooth, minPeakFracOfMax, maxPeaks,
                                                minPeakSepYpFp, minPeakSepYFp);
  cout << "Local maxima found = " << peaks.size() << endl;

  vector<LocalCandidate> candidates;
  for (size_t ip = 0; ip < peaks.size(); ip++) {
    LocalCandidate cand = BuildLocalEllipseCandidate(events, hSmooth,
                                                     peaks[ip].first, peaks[ip].second,
                                                     (int)ip,
                                                     seedHalfYpFp, seedHalfYFp,
                                                     nSigmaMajor, nSigmaMinor,
                                                     100);
    if (!cand.cut) continue;
    ScoreCandidate(cand, events, ysWindow, minSelected);
    candidates.push_back(cand);
  }

  sort(candidates.begin(), candidates.end(), [](const LocalCandidate &a, const LocalCandidate &b) {
    return a.quality > b.quality;
  });

  if (candidates.empty() || candidates[0].quality < 0.0) {
    cerr << "ERROR: no usable local candidate found. Try lowering minPeakFracOfMax or localRelLevel." << endl;
    return;
  }

  LocalCandidate best = candidates[0];
  best.cut->SetName(Form("best_candidate_yscol_%d_nfoil_%d_ndel_%d", best.yscol, foilIndex, ndel));

  cout << "Best candidate:" << endl;
  cout << "  peak(YpFp,YFp) = (" << best.peakX << ", " << best.peakY << ")" << endl;
  cout << "  assigned yscol = " << best.yscol << endl;
  cout << "  nSelected      = " << best.nSelected << endl;
  cout << "  bestScore      = " << best.bestScore << endl;
  cout << "  purity         = " << best.purity << endl;
  cout << "  quality        = " << best.quality << endl;
  cout << "  knobs          = minPeakFracOfMax=" << minPeakFracOfMax
       << " localRelLevel=" << localRelLevel
       << " seedHalfYpFp=" << seedHalfYpFp
       << " seedHalfYFp=" << seedHalfYFp
       << " nSigmaMajor=" << nSigmaMajor
       << " nSigmaMinor=" << nSigmaMinor << endl;

  TH2D *hYBest = new TH2D("hYpFpYFp_best_selected",
                          Form("Run %d best selected events;YpFp;YFp", nrun),
                          nBinsYpFp, ypfpMin, ypfpMax,
                          nBinsYFp, yfpMin, yfpMax);
  TH2D *hYsXsBest = new TH2D("hYsXs_best_selected",
                             Form("Run %d selected by best YpFp/YFp candidate;Ys;Xs", nrun),
                             100, -7.0, 7.0, 120, -12.5, 12.5);

  for (auto &e : events) {
    if (!best.cut->IsInside(e.ypfp, e.yfp)) continue;
    hYBest->Fill(e.ypfp, e.yfp);
    hYsXsBest->Fill(e.ys, e.xs);
  }

  gSystem->mkdir("plots", kTRUE);
  TString outPdf = Form("plots/test_yfp_ypfp_one_island_localellipse_run%d_%s.pdf", nrun, outTag.Data());

  TCanvas *c = new TCanvas("c_test_one_island", "test one normalized local YpFp/YFp island", 1100, 850);
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.026);

  c->Print(outPdf + "[");

  c->Clear();
  gPad->SetLogz(1);
  hYAll->Draw("colz");
  best.cut->SetLineColor(kRed);
  best.cut->SetLineWidth(3);
  best.cut->Draw("L same");
  TText *lab = new TText(best.mx, best.my, Form("cand -> Ys %d", best.yscol));
  lab->SetTextColor(kRed);
  lab->SetTextSize(0.035);
  lab->SetTextAlign(22);
  lab->Draw("same");
  tx.DrawLatex(0.12, 0.94, Form("All YpFp/YFp events; best local candidate overlaid | purity=%.3f n=%d", best.purity, best.nSelected));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  gPad->SetLogz(1);
  hYsXsAll->Draw("colz");
  DrawYsGuideLines();
  tx.DrawLatex(0.12, 0.94, "All events after PID + ytar cut + selected delta slice");
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  gPad->SetLogz(1);
  hYsXsBest->Draw("colz");
  DrawYsGuideLines();
  tx.DrawLatex(0.12, 0.94, Form("Selected events only: assigned Ys column %d | purity=%.3f n=%d score=%d", best.yscol, best.purity, best.nSelected, best.bestScore));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  gPad->SetLogz(1);
  hYBest->Draw("colz");
  best.cut->Draw("L same");
  tx.DrawLatex(0.12, 0.94, "YpFp/YFp events passing the best candidate cut only");
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Print(outPdf + "]");

  cout << "Wrote PDF: " << outPdf << endl;
  cout << "Suggested next tests:" << endl;
  cout << "  smaller seed: .x test_yfp_ypfp_one_island_localellipse.C(" << nrun << ",\"" << outTag << "_smallseed\",\"" << inputID << "\",\"" << ytarTag << "\"," << foilIndex << "," << ndel << "," << minPeakFracOfMax << "," << localRelLevel << ",0.0030,3.0," << nSigmaMajor << "," << nSigmaMinor << ")" << endl;
  cout << "  larger ellipse: .x test_yfp_ypfp_one_island_localellipse.C(" << nrun << ",\"" << outTag << "_larger\",\"" << inputID << "\",\"" << ytarTag << "\"," << foilIndex << "," << ndel << "," << minPeakFracOfMax << "," << localRelLevel << "," << seedHalfYpFp << "," << seedHalfYFp << ",2.0,1.7)" << endl;

  fin->Close();
}
