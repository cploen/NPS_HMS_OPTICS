// test_yfp_ypfp_axisvalley_coordinate_debug.C
//
// One-slice Yfp/YpFp island diagnostic using principal-axis coordinates
// and valley-bounded selection along the dense band.
//
// Purpose:
//   Work on the selection algorithm in the simplest possible loop:
//     one run, one foil, one delta slice, one selected Y-projection island.
//
// Method:
//   1. Apply PID + auto ytar foil cut + one delta slice from list_of_optics_run.dat.
//   2. Compute a normalized principal-axis basis in the YpFp/YFp plane.
//      This is only a coordinate system, not a global cut model.
//   3. Project events onto s = along-band coordinate and u = transverse coordinate.
//   4. Find peaks and valleys in the 1D s-density.
//   5. For each peak, define s boundaries by neighboring valleys.
//   6. Define a transverse u guard from robust quantiles inside that s interval.
//   7. Score each candidate using concentration near one nominal Ys guide line.
//   8. Choose the best candidate and make diagnostic plots.
//
// Usage:
//   hcana -l
//   .x test_yfp_ypfp_axisvalley_coordinate_debug.C(1544,"axiscoord_debug1","-1","multifoil_test1",0,2)
//
// Main knobs:
//   minPeakFracOfMax : threshold for s-density peak finding
//   minPeakSepS      : minimum separation between s peaks in normalized units
//   uSigmaScale      : transverse width multiplier from robust u sigma
//   minPurity        : minimum accepted Ys-line purity
//   maxNeighborFrac  : maximum accepted neighboring-Ys contamination

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

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
  double s = 0.0;
  double u = 0.0;
};

struct AxisBasis {
  double mx = 0.0;
  double my = 0.0;
  double xScale = 1.0;
  double yScale = 1.0;
  double vx = 1.0;
  double vy = 0.0;
  double ux = 0.0;
  double uy = 1.0;
};

struct SPeak {
  int index = -1;
  int bin = -1;
  double s = 0.0;
  double height = 0.0;
  double sLo = 0.0;
  double sHi = 0.0;
};

struct Candidate {
  int peakIndex = -1;
  double sPeak = 0.0;
  double sLo = 0.0;
  double sHi = 0.0;
  double uCenter = 0.0;
  double uHalf = 0.0;
  int yscol = -1;
  int bestScore = 0;
  int neighborScore = 0;
  int nSelected = 0;
  double purity = 0.0;
  double neighborFrac = 0.0;
  double quality = -1.0;
  TCutG *cut = nullptr;
};

bool HasBranch(TTree *T, const char *bname) { return T && T->GetBranch(bname); }

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

AxisBasis ComputeAxisBasis(const vector<FPEventY> &events) {
  AxisBasis b;
  vector<double> xs, ys;
  xs.reserve(events.size());
  ys.reserve(events.size());
  for (auto &e : events) {
    xs.push_back(e.ypfp);
    ys.push_back(e.yfp);
  }

  b.mx = Quantile(xs, 0.50);
  b.my = Quantile(ys, 0.50);
  b.xScale = 0.5 * (Quantile(xs, 0.84) - Quantile(xs, 0.16));
  b.yScale = 0.5 * (Quantile(ys, 0.84) - Quantile(ys, 0.16));
  if (b.xScale <= 1e-6) b.xScale = 0.006;
  if (b.yScale <= 1e-6) b.yScale = 6.0;

  double cxx = 0.0, cxy = 0.0, cyy = 0.0;
  int n = 0;
  for (auto &e : events) {
    double xn = (e.ypfp - b.mx) / b.xScale;
    double yn = (e.yfp  - b.my) / b.yScale;
    if (std::abs(xn) > 5.0 || std::abs(yn) > 5.0) continue;
    cxx += xn * xn;
    cxy += xn * yn;
    cyy += yn * yn;
    n++;
  }
  if (n > 0) {
    cxx /= n;
    cxy /= n;
    cyy /= n;
  }

  double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
  b.vx = cos(theta);
  b.vy = sin(theta);
  // Orient the s-axis to increase with YpFp for reproducible peak ordering.
  if (b.vx < 0.0) { b.vx *= -1.0; b.vy *= -1.0; }
  b.ux = -b.vy;
  b.uy =  b.vx;

  cout << "Axis basis:" << endl;
  cout << "  center(YpFp,YFp) = (" << b.mx << ", " << b.my << ")" << endl;
  cout << "  scales(YpFp,YFp) = (" << b.xScale << ", " << b.yScale << ")" << endl;
  cout << "  v_along = (" << b.vx << ", " << b.vy << ") in normalized coordinates" << endl;
  cout << "  u_trans = (" << b.ux << ", " << b.uy << ") in normalized coordinates" << endl;

  return b;
}

void FillSU(vector<FPEventY> &events, const AxisBasis &b) {
  for (auto &e : events) {
    double xn = (e.ypfp - b.mx) / b.xScale;
    double yn = (e.yfp  - b.my) / b.yScale;
    e.s = xn * b.vx + yn * b.vy;
    e.u = xn * b.ux + yn * b.uy;
  }
}

void PhysFromSU(const AxisBasis &b, double s, double u, double &x, double &y) {
  double xn = s * b.vx + u * b.ux;
  double yn = s * b.vy + u * b.uy;
  x = b.mx + b.xScale * xn;
  y = b.my + b.yScale * yn;
}

TCutG* MakeSUCut(const char *name, const AxisBasis &b,
                 double sLo, double sHi, double uLo, double uHi,
                 int color=kRed, int width=3) {
  TCutG *cut = new TCutG(name, 5);
  cut->SetTitle(Form("%s;YpFp;YFp", name));

  double x, y;
  PhysFromSU(b, sLo, uLo, x, y); cut->SetPoint(0, x, y);
  PhysFromSU(b, sHi, uLo, x, y); cut->SetPoint(1, x, y);
  PhysFromSU(b, sHi, uHi, x, y); cut->SetPoint(2, x, y);
  PhysFromSU(b, sLo, uHi, x, y); cut->SetPoint(3, x, y);
  PhysFromSU(b, sLo, uLo, x, y); cut->SetPoint(4, x, y);

  cut->SetLineColor(color);
  cut->SetLineWidth(width);
  return cut;
}


void DrawConstantSLine(const AxisBasis &b, double s, double uMin, double uMax,
                       int color, int style, int width) {
  double x1, y1, x2, y2;
  PhysFromSU(b, s, uMin, x1, y1);
  PhysFromSU(b, s, uMax, x2, y2);
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw("same");
}

void DrawConstantULine(const AxisBasis &b, double u, double sMin, double sMax,
                       int color, int style, int width) {
  double x1, y1, x2, y2;
  PhysFromSU(b, sMin, u, x1, y1);
  PhysFromSU(b, sMax, u, x2, y2);
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw("same");
}

void DrawAxisCoordinateGrid(const AxisBasis &basis,
                            const vector<SPeak> &peaks,
                            const vector<Candidate> &candidates,
                            double sMinDraw,
                            double sMaxDraw,
                            double uMinDraw,
                            double uMaxDraw) {
  // Red solid lines: s positions of density peaks.
  // Gray dashed lines: s interval boundaries used by candidate cuts.
  // Blue dotted lines: transverse u coordinates, to reveal whether the basis is sensible.

  for (auto &cnd : candidates) {
    DrawConstantSLine(basis, cnd.sLo, uMinDraw, uMaxDraw, kGray+2, 3, 1);
    DrawConstantSLine(basis, cnd.sHi, uMinDraw, uMaxDraw, kGray+2, 3, 1);
  }

  for (auto &p : peaks) {
    DrawConstantSLine(basis, p.s, uMinDraw, uMaxDraw, kRed, 1, 2);
    double xLab, yLab;
    PhysFromSU(basis, p.s, uMaxDraw*0.92, xLab, yLab);
    TText *txt = new TText(xLab, yLab, Form("p%d", p.index));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.026);
    txt->SetTextAlign(22);
    txt->Draw("same");
  }

  for (double u = -1.0; u <= 1.0001; u += 0.5) {
    DrawConstantULine(basis, u, sMinDraw, sMaxDraw, kBlue+1, 2, (std::abs(u) < 1e-9 ? 2 : 1));
    double xLab, yLab;
    PhysFromSU(basis, sMinDraw, u, xLab, yLab);
    TText *txt = new TText(xLab, yLab, Form("u=%.1f", u));
    txt->SetTextColor(kBlue+1);
    txt->SetTextSize(0.020);
    txt->SetTextAlign(12);
    txt->Draw("same");
  }
}

vector<SPeak> FindSPeaks(TH1D *hS, double minPeakFracOfMax, double minPeakSepS, int maxPeaks) {
  vector<SPeak> peaks;
  vector<SPeak> raw;
  if (!hS) return peaks;

  int nb = hS->GetNbinsX();
  double hmax = hS->GetMaximum();
  double minVal = minPeakFracOfMax * hmax;

  for (int b = 2; b <= nb - 1; b++) {
    double v = hS->GetBinContent(b);
    if (v < minVal) continue;
    if (v >= hS->GetBinContent(b-1) && v >= hS->GetBinContent(b+1)) {
      SPeak p;
      p.bin = b;
      p.s = hS->GetBinCenter(b);
      p.height = v;
      raw.push_back(p);
    }
  }

  sort(raw.begin(), raw.end(), [](const SPeak &a, const SPeak &b) {
    return a.height > b.height;
  });

  for (auto &p : raw) {
    bool tooClose = false;
    for (auto &q : peaks) {
      if (std::abs(p.s - q.s) < minPeakSepS) {
        tooClose = true;
        break;
      }
    }
    if (tooClose) continue;
    peaks.push_back(p);
    if ((int)peaks.size() >= maxPeaks) break;
  }

  sort(peaks.begin(), peaks.end(), [](const SPeak &a, const SPeak &b) {
    return a.s < b.s;
  });

  for (size_t i = 0; i < peaks.size(); i++) peaks[i].index = (int)i;
  return peaks;
}

int FindValleyBin(TH1D *hS, int b1, int b2) {
  if (!hS) return -1;
  if (b1 > b2) std::swap(b1, b2);
  int bestBin = b1;
  double bestVal = hS->GetBinContent(b1);
  for (int b = b1; b <= b2; b++) {
    double v = hS->GetBinContent(b);
    if (v < bestVal) {
      bestVal = v;
      bestBin = b;
    }
  }
  return bestBin;
}

vector<Candidate> BuildCandidatesFromSPeaks(const vector<FPEventY> &events,
                                            const vector<SPeak> &peaks,
                                            TH1D *hS,
                                            const AxisBasis &basis,
                                            double maxHalfWidthS,
                                            double uSigmaScale,
                                            double minUHalf,
                                            double maxUHalf,
                                            double ysWindow,
                                            int minSelected,
                                            double minPurity,
                                            double maxNeighborFrac) {
  vector<Candidate> out;
  if (peaks.empty()) return out;

  vector<double> boundaries;
  boundaries.resize(peaks.size() + 1);
  double sAllMin = hS->GetXaxis()->GetXmin();
  double sAllMax = hS->GetXaxis()->GetXmax();

  for (size_t i = 0; i + 1 < peaks.size(); i++) {
    int vb = FindValleyBin(hS, peaks[i].bin, peaks[i+1].bin);
    boundaries[i+1] = hS->GetBinCenter(vb);
  }

  // Edge intervals: use half the neighboring peak spacing, but do not exceed histogram range.
  if (peaks.size() == 1) {
    boundaries[0] = sAllMin;
    boundaries[1] = sAllMax;
  } else {
    double leftWidth = boundaries[1] - peaks[0].s;
    double rightWidth = peaks.back().s - boundaries[peaks.size()-1];
    boundaries[0] = std::max(sAllMin, peaks[0].s - std::abs(leftWidth));
    boundaries[peaks.size()] = std::min(sAllMax, peaks.back().s + std::abs(rightWidth));
  }

  for (size_t ip = 0; ip < peaks.size(); ip++) {
    Candidate cand;
    cand.peakIndex = peaks[ip].index;
    cand.sPeak = peaks[ip].s;
    cand.sLo = boundaries[ip];
    cand.sHi = boundaries[ip+1];

    // Guard against a bad valley producing one over-wide merged interval.
    // This caps the valley-to-valley region around the peak in s-space.
    if (maxHalfWidthS > 0.0) {
      cand.sLo = std::max(cand.sLo, cand.sPeak - maxHalfWidthS);
      cand.sHi = std::min(cand.sHi, cand.sPeak + maxHalfWidthS);
    }

    vector<double> uVals;
    uVals.reserve(5000);
    for (auto &e : events) {
      if (e.s < cand.sLo || e.s > cand.sHi) continue;
      uVals.push_back(e.u);
    }
    if (uVals.size() < 20) continue;

    cand.uCenter = Quantile(uVals, 0.50);
    double u16 = Quantile(uVals, 0.16);
    double u84 = Quantile(uVals, 0.84);
    double sigU = 0.5 * (u84 - u16);
    cand.uHalf = uSigmaScale * sigU;
    if (cand.uHalf < minUHalf) cand.uHalf = minUHalf;
    if (cand.uHalf > maxUHalf) cand.uHalf = maxUHalf;

    cand.cut = MakeSUCut(Form("axisvalley_candidate_peak%d", cand.peakIndex), basis,
                         cand.sLo, cand.sHi,
                         cand.uCenter - cand.uHalf,
                         cand.uCenter + cand.uHalf,
                         kRed, 3);

    vector<int> score(9, 0);
    cand.nSelected = 0;
    for (auto &e : events) {
      if (e.s < cand.sLo || e.s > cand.sHi) continue;
      if (std::abs(e.u - cand.uCenter) > cand.uHalf) continue;
      cand.nSelected++;
      for (int iy = 0; iy < 9; iy++) {
        double ysNom = (iy - 4) * 0.6 * 2.54;
        if (std::abs(e.ys - ysNom) < ysWindow) score[iy]++;
      }
    }

    int best = -1;
    int bestScore = -1;
    for (int iy = 0; iy < 9; iy++) {
      if (score[iy] > bestScore) { bestScore = score[iy]; best = iy; }
    }
    int neighborScore = 0;
    if (best > 0) neighborScore += score[best-1];
    if (best < 8) neighborScore += score[best+1];

    cand.yscol = best;
    cand.bestScore = bestScore;
    cand.neighborScore = neighborScore;
    cand.purity = cand.nSelected > 0 ? (double)cand.bestScore / (double)cand.nSelected : 0.0;
    cand.neighborFrac = cand.nSelected > 0 ? (double)cand.neighborScore / (double)cand.nSelected : 0.0;
    if (cand.nSelected < minSelected) cand.quality = -1.0;
    else cand.quality = cand.purity * sqrt((double)cand.bestScore) * (1.0 - cand.neighborFrac);

    out.push_back(cand);
  }

  return out;
}

bool CandidateIsAccepted(const Candidate &c, int minSelected,
                         double minPurity, double maxNeighborFrac) {
  if (c.nSelected < minSelected) return false;
  if (c.purity < minPurity) return false;
  if (c.neighborFrac > maxNeighborFrac) return false;
  return true;
}

void test_yfp_ypfp_axisvalley_coordinate_debug(Int_t nrun=1544,
                              const char *tag="axiscoord_debug1",
                              const char *inputFileID="-1",
                              const char *ytarCutTag="multifoil_test1",
                              Int_t foilIndex=0,
                              Int_t ndel=2,
                              Double_t minPeakFracOfMax=0.055,
                              Double_t minPeakSepS=0.45,
                              Int_t smoothPassesS=1,
                              Double_t maxHalfWidthS=0.55,
                              Double_t uSigmaScale=1.55,
                              Double_t minUHalf=0.20,
                              Double_t maxUHalf=0.85,
                              Double_t minPurity=0.70,
                              Double_t maxNeighborFrac=0.18)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  const double npeMin = 2.0;
  const int nBinsYpFp = 180;
  const double ypfpMin = -0.035;
  const double ypfpMax =  0.035;
  const int nBinsYFp = 180;
  const double yfpMin = -35.0;
  const double yfpMax =  35.0;
  const int nBinsS = 220;
  const int maxPeaks = 12;
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
  if (!fin || fin->IsZombie()) { cerr << "ERROR: could not open input ROOT file: " << inputroot << endl; return; }
  TTree *T = (TTree*)fin->Get("T");
  if (!T) { cerr << "ERROR: could not find tree T" << endl; return; }

  const char *required[] = {
    "H.cer.npeSum", "H.gtr.dp", "H.gtr.y",
    "H.dc.y_fp", "H.dc.yp_fp", "H.extcor.ysieve", "H.extcor.xsieve"
  };
  for (auto b : required) {
    if (!HasBranch(T, b)) { cerr << "ERROR: missing branch " << b << endl; return; }
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
                         nBinsYpFp, ypfpMin, ypfpMax, nBinsYFp, yfpMin, yfpMax);
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
    ev.ypfp = ypfp; ev.yfp = yfp; ev.xs = xs; ev.ys = ys; ev.delta = delta; ev.ytar = ytar;
    events.push_back(ev);
    hYAll->Fill(ev.ypfp, ev.yfp);
    hYsXsAll->Fill(ev.ys, ev.xs);
  }

  cout << "Events after PID + ytar + delta = " << events.size() << endl;
  if (events.size() < 50) { cerr << "ERROR: too few events in this slice." << endl; return; }

  AxisBasis basis = ComputeAxisBasis(events);
  FillSU(events, basis);

  vector<double> sVals;
  sVals.reserve(events.size());
  for (auto &e : events) sVals.push_back(e.s);
  double sMin = Quantile(sVals, 0.005);
  double sMax = Quantile(sVals, 0.995);
  double margin = 0.05 * (sMax - sMin);
  sMin -= margin;
  sMax += margin;

  TH1D *hS = new TH1D("hS_density", Form("Run %d ndel %d: density along principal axis;s coordinate;counts", nrun, ndel), nBinsS, sMin, sMax);
  for (auto &e : events) hS->Fill(e.s);
  for (int is = 0; is < smoothPassesS; is++) hS->Smooth(1);

  vector<SPeak> peaks = FindSPeaks(hS, minPeakFracOfMax, minPeakSepS, maxPeaks);
  vector<Candidate> candidates = BuildCandidatesFromSPeaks(events, peaks, hS, basis,
                                                           maxHalfWidthS,
                                                           uSigmaScale, minUHalf, maxUHalf,
                                                           ysWindow, minSelected,
                                                           minPurity, maxNeighborFrac);

  vector<double> uValsForDraw;
  uValsForDraw.reserve(events.size());
  for (auto &e : events) uValsForDraw.push_back(e.u);
  double uMinDraw = Quantile(uValsForDraw, 0.005);
  double uMaxDraw = Quantile(uValsForDraw, 0.995);
  double uMargin = 0.08 * (uMaxDraw - uMinDraw);
  uMinDraw -= uMargin;
  uMaxDraw += uMargin;

  cout << "S-axis peaks found = " << peaks.size() << endl;
  cout << "Axis-valley candidate summary:" << endl;
  cout << "  peak sPeak sLo sHi nSelected yscol purity neighborFrac quality accepted" << endl;
  for (auto &cnd : candidates) {
    bool ok = CandidateIsAccepted(cnd, minSelected, minPurity, maxNeighborFrac);
    cout << "  " << cnd.peakIndex
         << " " << cnd.sPeak
         << " " << cnd.sLo
         << " " << cnd.sHi
         << " " << cnd.nSelected
         << " " << cnd.yscol
         << " " << cnd.purity
         << " " << cnd.neighborFrac
         << " " << cnd.quality
         << " " << (ok ? "YES" : "no")
         << endl;
  }

  if (candidates.empty()) {
    cerr << "ERROR: no candidates were built. Try lowering minPeakFracOfMax or minPeakSepS." << endl;
    return;
  }

  Candidate best = candidates[0];
  for (auto &cnd : candidates) {
    if (!CandidateIsAccepted(cnd, minSelected, minPurity, maxNeighborFrac)) continue;
    if (cnd.quality > best.quality) best = cnd;
  }
  if (!CandidateIsAccepted(best, minSelected, minPurity, maxNeighborFrac)) {
    sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) { return a.quality > b.quality; });
    best = candidates[0];
    cout << "WARNING: no candidate passed acceptance criteria; using highest-quality fallback." << endl;
  }

  cout << "Best candidate:" << endl;
  cout << "  peak index     = " << best.peakIndex << endl;
  cout << "  s interval     = [" << best.sLo << ", " << best.sHi << "]" << endl;
  cout << "  u center/half  = " << best.uCenter << " +/- " << best.uHalf << endl;
  cout << "  assigned yscol = " << best.yscol << endl;
  cout << "  nSelected      = " << best.nSelected << endl;
  cout << "  purity         = " << best.purity << endl;
  cout << "  neighborFrac   = " << best.neighborFrac << endl;
  cout << "  quality        = " << best.quality << endl;
  cout << "  knobs          = minPeakFracOfMax=" << minPeakFracOfMax
       << " minPeakSepS=" << minPeakSepS
       << " smoothPassesS=" << smoothPassesS
       << " maxHalfWidthS=" << maxHalfWidthS
       << " uSigmaScale=" << uSigmaScale
       << " minUHalf=" << minUHalf
       << " maxUHalf=" << maxUHalf
       << " minPurity=" << minPurity
       << " maxNeighborFrac=" << maxNeighborFrac << endl;

  TH2D *hYBest = new TH2D("hYpFpYFp_best_selected",
                          Form("Run %d best selected events;YpFp;YFp", nrun),
                          nBinsYpFp, ypfpMin, ypfpMax, nBinsYFp, yfpMin, yfpMax);
  TH2D *hYsXsBest = new TH2D("hYsXs_best_selected",
                             Form("Run %d selected by best axis-valley candidate;Ys;Xs", nrun),
                             100, -7.0, 7.0, 120, -12.5, 12.5);

  for (auto &e : events) {
    if (e.s < best.sLo || e.s > best.sHi) continue;
    if (std::abs(e.u - best.uCenter) > best.uHalf) continue;
    hYBest->Fill(e.ypfp, e.yfp);
    hYsXsBest->Fill(e.ys, e.xs);
  }

  gSystem->mkdir("plots", kTRUE);
  TString outPdf = Form("plots/test_yfp_ypfp_axisvalley_coordinate_debug_run%d_%s.pdf", nrun, outTag.Data());

  TCanvas *c = new TCanvas("c_axisvalley_debug", "axis-valley coordinate diagnostic", 1100, 850);
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.026);

  c->Print(outPdf + "[");

  c->Clear();
  gPad->SetLogz(1);
  hYAll->Draw("colz");
  DrawAxisCoordinateGrid(basis, peaks, candidates, sMin, sMax, uMinDraw, uMaxDraw);
  tx.DrawLatex(0.10, 0.94, "Principal-axis coordinate diagnostic in YpFp/YFp: red=s peaks, gray=s boundaries, blue=u lines");
  tx.DrawLatex(0.10, 0.90, Form("basis center=(%.5f, %.2f), scales=(%.5f, %.2f), v=(%.3f, %.3f), u=(%.3f, %.3f)", basis.mx, basis.my, basis.xScale, basis.yScale, basis.vx, basis.vy, basis.ux, basis.uy));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  gPad->SetLogz(1);
  hYAll->Draw("colz");
  for (auto &cnd : candidates) {
    if (!cnd.cut) continue;
    cnd.cut->SetLineColor(kGray+1);
    cnd.cut->SetLineWidth(1);
    cnd.cut->Draw("L same");
    double xLab, yLab;
    PhysFromSU(basis, 0.5*(cnd.sLo+cnd.sHi), cnd.uCenter, xLab, yLab);
    TText *lab = new TText(xLab, yLab, Form("%d", cnd.peakIndex));
    lab->SetTextColor(kBlack);
    lab->SetTextSize(0.030);
    lab->SetTextAlign(22);
    lab->Draw("same");
  }
  if (best.cut) {
    best.cut->SetLineColor(kRed);
    best.cut->SetLineWidth(3);
    best.cut->Draw("L same");
  }
  tx.DrawLatex(0.12, 0.94, Form("All YpFp/YFp events; s-valley candidates overlaid | peaks=%zu best peak=%d Ys=%d purity=%.3f n=%d", peaks.size(), best.peakIndex, best.yscol, best.purity, best.nSelected));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  hS->SetLineColor(kBlue+1);
  hS->SetLineWidth(2);
  hS->Draw("hist");
  double ymax = hS->GetMaximum();
  for (auto &p : peaks) {
    TLine *lp = new TLine(p.s, 0.0, p.s, ymax);
    lp->SetLineColor(kRed);
    lp->SetLineStyle(2);
    lp->SetLineWidth(2);
    lp->Draw("same");
    TText *txt = new TText(p.s, 0.88*ymax, Form("%d", p.index));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.030);
    txt->SetTextAlign(22);
    txt->Draw("same");
  }
  for (auto &cnd : candidates) {
    TLine *ll = new TLine(cnd.sLo, 0.0, cnd.sLo, 0.70*ymax);
    TLine *rr = new TLine(cnd.sHi, 0.0, cnd.sHi, 0.70*ymax);
    ll->SetLineColor(kGray+2); rr->SetLineColor(kGray+2);
    ll->SetLineStyle(3); rr->SetLineStyle(3);
    ll->Draw("same"); rr->Draw("same");
  }
  tx.DrawLatex(0.12, 0.94, "1D density along principal-axis coordinate s; red=peaks, gray=valley boundaries");
  c->Print(outPdf);

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
  tx.DrawLatex(0.12, 0.94, Form("Selected events only: best peak %d assigned Ys %d | purity=%.3f neighbor=%.3f n=%d score=%d", best.peakIndex, best.yscol, best.purity, best.neighborFrac, best.nSelected, best.bestScore));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  gPad->SetLogz(1);
  hYBest->Draw("colz");
  if (best.cut) best.cut->Draw("L same");
  tx.DrawLatex(0.12, 0.94, "YpFp/YFp events passing best s-valley candidate only");
  c->Print(outPdf);
  gPad->SetLogz(0);



  // Candidate-by-candidate diagnostic pages.
  // Each candidate gets a YpFp/YFp view and an Xs/Ys validation view.
  for (auto &cnd : candidates) {
    bool ok = CandidateIsAccepted(cnd, minSelected, minPurity, maxNeighborFrac);

    TH2D *hYCand = new TH2D(Form("hYpFpYFp_candidate_peak%d", cnd.peakIndex),
                            Form("Run %d peak %d candidate;YpFp;YFp", nrun, cnd.peakIndex),
                            nBinsYpFp, ypfpMin, ypfpMax, nBinsYFp, yfpMin, yfpMax);
    TH2D *hYsXsCand = new TH2D(Form("hYsXs_candidate_peak%d", cnd.peakIndex),
                               Form("Run %d peak %d candidate;Ys;Xs", nrun, cnd.peakIndex),
                               100, -7.0, 7.0, 120, -12.5, 12.5);

    for (auto &e : events) {
      if (e.s < cnd.sLo || e.s > cnd.sHi) continue;
      if (std::abs(e.u - cnd.uCenter) > cnd.uHalf) continue;
      hYCand->Fill(e.ypfp, e.yfp);
      hYsXsCand->Fill(e.ys, e.xs);
    }

    c->Clear();
    gPad->SetLogz(1);
    hYAll->Draw("colz");
    if (cnd.cut) {
      cnd.cut->SetLineColor(ok ? kGreen+2 : kOrange+7);
      cnd.cut->SetLineWidth(3);
      cnd.cut->Draw("L same");
    }
    // Mark the density peak location in physical coordinates.
    double xPeak, yPeak;
    PhysFromSU(basis, cnd.sPeak, cnd.uCenter, xPeak, yPeak);
    TText *peakLabel = new TText(xPeak, yPeak, Form("peak %d", cnd.peakIndex));
    peakLabel->SetTextColor(kBlack);
    peakLabel->SetTextSize(0.030);
    peakLabel->SetTextAlign(22);
    peakLabel->Draw("same");
    tx.DrawLatex(0.11, 0.94, Form("Candidate peak %d in YpFp/YFp | s=[%.3f, %.3f] u=%.3f #pm %.3f | Ys=%d purity=%.3f neighbor=%.3f n=%d %s",
                                   cnd.peakIndex, cnd.sLo, cnd.sHi, cnd.uCenter, cnd.uHalf,
                                   cnd.yscol, cnd.purity, cnd.neighborFrac, cnd.nSelected,
                                   ok ? "ACCEPT" : "REJECT"));
    c->Print(outPdf);
    gPad->SetLogz(0);

    c->Clear();
    gPad->SetLogz(1);
    hYsXsCand->Draw("colz");
    DrawYsGuideLines();
    tx.DrawLatex(0.11, 0.94, Form("Candidate peak %d projected to Xs/Ys | assigned Ys=%d score=%d n=%d purity=%.3f neighbor=%.3f %s",
                                   cnd.peakIndex, cnd.yscol, cnd.bestScore, cnd.nSelected,
                                   cnd.purity, cnd.neighborFrac,
                                   ok ? "ACCEPT" : "REJECT"));
    c->Print(outPdf);
    gPad->SetLogz(0);
  }

  c->Print(outPdf + "]");

  cout << "Wrote PDF: " << outPdf << endl;
  cout << "Suggested diagnostic variations:" << endl;
  cout << "  find more peaks: .x test_yfp_ypfp_axisvalley_coordinate_debug.C(" << nrun << ",\"" << outTag << "_morePeaks\",\"" << inputID << "\",\"" << ytarTag << "\"," << foilIndex << "," << ndel << ",0.035,0.30,0," << uSigmaScale << "," << minUHalf << "," << maxUHalf << "," << minPurity << "," << maxNeighborFrac << ")" << endl;
  cout << "  narrower transverse guard: .x test_yfp_ypfp_axisvalley_coordinate_debug.C(" << nrun << ",\"" << outTag << "_narrowU\",\"" << inputID << "\",\"" << ytarTag << "\"," << foilIndex << "," << ndel << "," << minPeakFracOfMax << "," << minPeakSepS << "," << smoothPassesS << ",1.25," << minUHalf << "," << maxUHalf << "," << minPurity << "," << maxNeighborFrac << ")" << endl;

  fin->Close();
}
