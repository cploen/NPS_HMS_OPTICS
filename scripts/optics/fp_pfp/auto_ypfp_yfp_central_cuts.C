// auto_ypfp_yfp_central_cuts.C
//
// First-pass automatic YpFp/YFp cut builder for HMS/NPS optics.
//
// Purpose:
//   Build the YpFp-vs-YFp graphical cuts that replace the hand-drawn
//   set_ypfp_yfp_cuts.C workflow, starting with one foil and the central
//   delta slice from DATfiles/list_of_optics_run.dat.
//
// Method:
//   1. Read run metadata and delta slices from list_of_optics_run.dat.
//   2. Select the central delta slice: the slice containing delta=0,
//      otherwise the slice with center closest to zero.
//   3. Load the automatic ytar-vs-delta foil cuts.
//   4. Apply base production cuts + selected foil ytar cut + central delta slice.
//   5. Treat YpFp/YFp as an ordered ridge problem:
//        - compute a PCA ridge axis in (YpFp, YFp),
//        - project events onto ridge coordinate s,
//        - find density peaks along s,
//        - assign each peak to the nearest forward-matrix predicted yscol.
//   6. Build simple rotated-rectangle TCutG objects around each accepted peak.
//   7. Write ROOT cuts using legacy-compatible names expected by make_hist_hms_optics.C:
//        hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d
//   8. Write diagnostic PDF/CSV, including one Xs/Ys selected-event page per yscol.
//
// Usage from repo top level:
//   root -l -b -q 'auto_ypfp_yfp_central_cuts.C(1544,"auto_yfp_central","-1","auto_ycut",0)'
//
// Then remake histograms using the produced YpFp/YFp cuts, e.g.:
//   root -l -b -q 'make_hist_hms_optics.C(1544,kTRUE,kTRUE,kFALSE,-1)'
//   root -l -b -q 'plot_yfp_cuts.C(1544,-1)'
//
// Notes:
//   - This is intentionally simple. It is a bridge macro, not the final quality scorer.
//   - It only builds cuts for one foil and the central delta slice by default.
//   - It uses unsupervised/statistical ridge-density discovery, then physics-guided indexing.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "TText.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"

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

struct FPEvent {
  double ypfp = 0.0;   // x-axis in old histograms
  double yfp = 0.0;    // y-axis in old histograms
  double xs = 0.0;
  double ys = 0.0;
  double delta = 0.0;
  double ytar = 0.0;
  double s = 0.0;      // along-ridge coordinate
  double u = 0.0;      // perpendicular coordinate
};

struct Peak1D {
  int bin = -1;
  double s = 0.0;
  double height = 0.0;
  double frac = 0.0;
  int yscol = -1;
  double sPred = 0.0;
  double sLo = 0.0;
  double sHi = 0.0;
  double uLo = 0.0;
  double uHi = 0.0;
  long nEvents = 0;
  TCutG *cut = nullptr;
  TH2D *hYsXs = nullptr;
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

    for (auto &z : ztok) {
      if (z.Length() > 0) info.ztarFoil.push_back(z.Atof());
    }
    for (auto &d : dtok) {
      if (d.Length() > 0) info.delcut.push_back(d.Atof());
    }

    if ((int)info.ztarFoil.size() < info.numFoil) {
      cerr << "ERROR: run " << nrun << " has fewer ztar foil values than NumFoil." << endl;
      return false;
    }
    // In the NPS_HMS_OPTICS DAT file, the header value called ndelcut is
    // not always the number of usable delta intervals. The actual usable
    // intervals are determined by the number of numeric boundaries on the
    // following line: N boundaries -> N-1 intervals.
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

int FindCentralDeltaIndex(const OpticsRunInfo &info) {
  int best = -1;
  double bestAbs = 1e99;

  for (int nd = 0; nd < info.ndelcut; nd++) {
    double lo = info.delcut[nd];
    double hi = info.delcut[nd+1];
    double cen = 0.5 * (lo + hi);

    if (lo <= 0.0 && 0.0 < hi) return nd;

    if (std::abs(cen) < bestAbs) {
      bestAbs = std::abs(cen);
      best = nd;
    }
  }
  return best;
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

  // New auto-ridge macro output style.
  candidates.push_back(Form("cuts/ytar_delta_%d_%s_multifoil_cut.root",
                            info.run, ytarCutTag.Data()));
  candidates.push_back(Form("cuts/ytar_delta_%d_%s_cut.root",
                            info.run, ytarCutTag.Data()));

  // Legacy make_hist style.
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

void ComputePCAAxes(const vector<FPEvent> &events,
                    double &mx, double &my,
                    double &vx, double &vy,
                    double &ux, double &uy) {
  mx = 0.0; my = 0.0;
  for (const auto &e : events) { mx += e.ypfp; my += e.yfp; }
  mx /= events.size();
  my /= events.size();

  double cxx = 0.0, cxy = 0.0, cyy = 0.0;
  for (const auto &e : events) {
    double dx = e.ypfp - mx;
    double dy = e.yfp - my;
    cxx += dx * dx;
    cxy += dx * dy;
    cyy += dy * dy;
  }
  cxx /= events.size();
  cxy /= events.size();
  cyy /= events.size();

  double theta = 0.5 * atan2(2.0 * cxy, cxx - cyy);
  vx = cos(theta);
  vy = sin(theta);

  // Choose sign so s increases roughly with Yfp. This keeps ordering less surprising.
  if (vy < 0.0) { vx *= -1.0; vy *= -1.0; }

  ux = -vy;
  uy =  vx;

  cout << "PCA ridge axis: mean=(" << mx << ", " << my << ")"
       << " v=(" << vx << ", " << vy << ")"
       << " u=(" << ux << ", " << uy << ")" << endl;
}

void AssignSU(vector<FPEvent> &events,
              double mx, double my,
              double vx, double vy,
              double ux, double uy) {
  for (auto &e : events) {
    double dx = e.ypfp - mx;
    double dy = e.yfp - my;
    e.s = dx * vx + dy * vy;
    e.u = dx * ux + dy * uy;
  }
}

vector<Peak1D> FindPeaks1D(TH1D *h,
                           double minPeakFrac,
                           double minPeakSeparationS) {
  vector<Peak1D> raw, out;
  if (!h) return out;

  double hmax = h->GetMaximum();
  if (hmax <= 0.0) return out;

  for (int b = 2; b <= h->GetNbinsX() - 1; b++) {
    double v  = h->GetBinContent(b);
    double vl = h->GetBinContent(b - 1);
    double vr = h->GetBinContent(b + 1);
    if (v < minPeakFrac * hmax) continue;
    if (v >= vl && v >= vr) {
      Peak1D p;
      p.bin = b;
      p.s = h->GetBinCenter(b);
      p.height = v;
      p.frac = v / hmax;
      raw.push_back(p);
    }
  }

  sort(raw.begin(), raw.end(), [](const Peak1D &a, const Peak1D &b) {
    return a.height > b.height;
  });

  for (auto &p : raw) {
    bool tooClose = false;
    for (auto &q : out) {
      if (std::abs(p.s - q.s) < minPeakSeparationS) {
        tooClose = true;
        break;
      }
    }
    if (!tooClose) out.push_back(p);
  }

  sort(out.begin(), out.end(), [](const Peak1D &a, const Peak1D &b) {
    return a.s < b.s;
  });

  return out;
}

bool LoadForwardMatrix(const char *fname,
                       vector<double> &xfpcoeffs,
                       vector<double> &xpfpcoeffs,
                       vector<double> &yfpcoeffs,
                       vector<double> &ypfpcoeffs,
                       vector<double> &lencoeffs,
                       vector<int> &xtarexpon,
                       vector<int> &xptarexpon,
                       vector<int> &ytarexpon,
                       vector<int> &yptarexpon,
                       vector<int> &deltaexpon) {
  ifstream oldcoeffsfile(fname);
  if (!oldcoeffsfile.is_open()) {
    cerr << "WARNING: cannot open forward matrix file: " << fname << endl;
    return false;
  }

  TString currentline;
  while (currentline.ReadLine(oldcoeffsfile, kFALSE) && !currentline.BeginsWith(" ----")) {
    if (currentline.Length() < 78) continue;

    TString sc1(currentline(1,14));
    TString sc2(currentline(15,14));
    TString sc3(currentline(29,14));
    TString sc4(currentline(43,14));
    TString sc5(currentline(57,14));

    xfpcoeffs.push_back(sc1.Atof());
    xpfpcoeffs.push_back(sc2.Atof());
    yfpcoeffs.push_back(sc3.Atof());
    ypfpcoeffs.push_back(sc4.Atof());
    lencoeffs.push_back(sc5.Atof());

    int expontemp[6];
    for (int expon = 0; expon < 6; expon++) {
      TString stemp(currentline(72 + expon, 1));
      expontemp[expon] = stemp.Atoi();
    }

    xtarexpon.push_back(expontemp[0]);
    xptarexpon.push_back(expontemp[1]);
    ytarexpon.push_back(expontemp[2]);
    yptarexpon.push_back(expontemp[3]);
    deltaexpon.push_back(expontemp[5]);
  }

  cout << "Loaded forward matrix terms: " << yfpcoeffs.size() << endl;
  return !yfpcoeffs.empty();
}

void PredictYFpYpFpGuides(const OpticsRunInfo &info,
                          int foilIndex,
                          int ndCentral,
                          double mx, double my,
                          double vx, double vy,
                          vector<double> &predYpFp,
                          vector<double> &predYFp,
                          vector<double> &predS) {
  predYpFp.assign(9, 0.0);
  predYFp.assign(9, 0.0);
  predS.assign(9, 0.0);

  vector<double> xfpcoeffs, xpfpcoeffs, yfpcoeffs, ypfpcoeffs, lencoeffs;
  vector<int> xtarexpon, xptarexpon, ytarexpon, yptarexpon, deltaexpon;

  bool haveMatrix = LoadForwardMatrix("hms_forward.dat",
                                      xfpcoeffs, xpfpcoeffs,
                                      yfpcoeffs, ypfpcoeffs, lencoeffs,
                                      xtarexpon, xptarexpon, ytarexpon,
                                      yptarexpon, deltaexpon);

  double zdis_sieve = 168.0;
  double theta = info.centAngleDeg * TMath::Pi() / 180.0;
  double DelCent = 0.5 * (info.delcut[ndCentral] + info.delcut[ndCentral+1]);
  double zfoil = info.ztarFoil[foilIndex];

  for (int nys = 0; nys < 9; nys++) {
    double ys_cent = (nys - 4) * 0.6 * 2.54;

    if (!haveMatrix) {
      // Fallback indexing only: monotonic placeholder in s.
      predYpFp[nys] = 0.0;
      predYFp[nys] = ys_cent;
      predS[nys] = (predYpFp[nys] - mx) * vx + (predYFp[nys] - my) * vy;
      continue;
    }

    double xtar = 0.0;
    double xptar = 0.0;
    double yptar = (ys_cent - zfoil * TMath::Sin(theta)) /
                   (zdis_sieve - zfoil * TMath::Cos(theta));
    double ytar = zfoil * (TMath::Sin(theta) - yptar * TMath::Cos(theta));

    double yfp = 0.0;
    double ypfp = 0.0;

    for (size_t i = 0; i < yfpcoeffs.size(); i++) {
      double etemp = pow(xtar, xtarexpon[i]) *
                     pow(xptar, xptarexpon[i]) *
                     pow(ytar, ytarexpon[i]) *
                     pow(yptar * 1000.0, yptarexpon[i]) *
                     pow(DelCent, deltaexpon[i]);
      yfp  += yfpcoeffs[i]  * etemp;
      ypfp += ypfpcoeffs[i] * etemp;
    }

    predYpFp[nys] = ypfp / 1000.0;  // match old macro guide-line units
    predYFp[nys]  = yfp;
    predS[nys]    = (predYpFp[nys] - mx) * vx + (predYFp[nys] - my) * vy;
  }
}

void AssignPeaksToYcols(vector<Peak1D> &peaks,
                        const vector<double> &predS,
                        double maxPredDistanceS) {
  set<int> usedPeak;
  set<int> usedY;

  struct Match { double d; int ip; int iy; };
  vector<Match> matches;

  for (int ip = 0; ip < (int)peaks.size(); ip++) {
    for (int iy = 0; iy < 9; iy++) {
      double d = std::abs(peaks[ip].s - predS[iy]);
      if (d <= maxPredDistanceS) matches.push_back({d, ip, iy});
    }
  }

  sort(matches.begin(), matches.end(), [](const Match &a, const Match &b) {
    return a.d < b.d;
  });

  for (auto &m : matches) {
    if (usedPeak.count(m.ip) || usedY.count(m.iy)) continue;
    peaks[m.ip].yscol = m.iy;
    peaks[m.ip].sPred = predS[m.iy];
    usedPeak.insert(m.ip);
    usedY.insert(m.iy);
  }

  sort(peaks.begin(), peaks.end(), [](const Peak1D &a, const Peak1D &b) {
    return a.yscol < b.yscol;
  });
}

TCutG* MakeRotatedRectangleCut(const char *name,
                               double mx, double my,
                               double vx, double vy,
                               double ux, double uy,
                               double sLo, double sHi,
                               double uLo, double uHi) {
  TCutG *cut = new TCutG(name, 5);
  cut->SetTitle(Form("%s;YpFp;YFp", name));

  auto setPoint = [&](int i, double s, double u) {
    double x = mx + s * vx + u * ux;
    double y = my + s * vy + u * uy;
    cut->SetPoint(i, x, y);
  };

  setPoint(0, sLo, uLo);
  setPoint(1, sHi, uLo);
  setPoint(2, sHi, uHi);
  setPoint(3, sLo, uHi);
  setPoint(4, sLo, uLo);

  cut->SetLineColor(kBlack);
  cut->SetLineWidth(3);
  return cut;
}

void DrawYsGuideLines() {
  for (int nys = 0; nys < 9; nys++) {
    double pos = (nys - 4) * 0.6 * 2.54;
    TLine *line = new TLine(pos, -12.5, pos, 12.5);
    line->SetLineColor(kRed);
    line->SetLineWidth(1);
    line->Draw("same");

    TText *txt = new TText(pos, -14.2, Form("%d", nys));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.03);
    txt->Draw("same");
  }
}

void auto_ypfp_yfp_central_cuts(Int_t nrun = 1544,
                                const char *tag = "auto_yfp_central",
                                const char *inputFileID = "-1",
                                const char *ytarCutTag = "auto_ycut",
                                Int_t foilIndex = 0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TString outTag = tag;
  TString inputID = inputFileID;
  TString ytarTag = ytarCutTag;

  const double baseDeltaMin = -10.0;
  const double baseDeltaMax =  10.0;
  const double npeMin = 2.0;

  const int    sBins = 240;
  const int    smoothPasses = 3;
  const double minPeakFrac = 0.08;
  const double minPeakSeparationS = 1.0;
  const double maxPredDistanceS = 4.5;
  const double sHalfWidthFallback = 0.70;
  const double sGuardScale = 1.10;
  const double uQuantileLow = 0.02;
  const double uQuantileHigh = 0.98;
  const double uPadMin = 0.0008;
  const double uPadFrac = 0.20;
  const int minEventsPerCut = 25;

  OpticsRunInfo info;
  if (!ReadOpticsRunInfo(nrun, info)) return;

  if (foilIndex < 0 || foilIndex >= info.numFoil) {
    cerr << "ERROR: foilIndex " << foilIndex << " outside NumFoil=" << info.numFoil << endl;
    return;
  }

  int ndCentral = FindCentralDeltaIndex(info);
  if (ndCentral < 0) {
    cerr << "ERROR: could not determine central delta slice." << endl;
    return;
  }

  double dLo = info.delcut[ndCentral];
  double dHi = info.delcut[ndCentral+1];
  double dC  = 0.5 * (dLo + dHi);

  cout << "Using central delta slice nd=" << ndCentral
       << " range=[" << dLo << ", " << dHi << "] center=" << dC << endl;

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

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchAddress("H.extcor.ysieve", &ys);
  T->SetBranchAddress("H.extcor.xsieve", &xs);

  vector<FPEvent> events;

  TH2D *hYpFpYFp = new TH2D("hYpFpYFp_selected_base",
                            Form("Run %d foil %d central #delta;YpFp;YFp", nrun, foilIndex),
                            120, -0.035, 0.035, 180, -40.0, 40.0);
  TH2D *hYsXsBase = new TH2D("hYsXs_base",
                             Form("Run %d foil %d central #delta;Ys;Xs", nrun, foilIndex),
                             100, -7.0, 7.0, 120, -12.5, 12.5);

  Long64_t nentries = T->GetEntries();

  cout << "Tree entries = " << nentries << endl;
  cout << "Starting central-delta + ytar-selected event loop..." << endl;

  Long64_t nPrintEvery = 100000;

  for (Long64_t i = 0; i < nentries; i++) {
    if (i > 0 && i % nPrintEvery == 0) {
      cout << "  processed event " << i << " / " << nentries << endl;
    }
  cout << "Finished event loop." << endl;
  cout << "Events after central-delta + ytar cut = " << events.size() << endl;
    T->GetEntry(i);

    if (!(sumnpe > npeMin)) continue;
    if (!(delta > baseDeltaMin && delta < baseDeltaMax)) continue;
    if (!(delta >= dLo && delta < dHi)) continue;
    if (!std::isfinite(ytar) || !std::isfinite(delta)) continue;
    if (!std::isfinite(yfp) || !std::isfinite(ypfp)) continue;
    if (!ytarCut->IsInside(ytar, delta)) continue;

    FPEvent ev;
    ev.ypfp = ypfp;
    ev.yfp = yfp;
    ev.xs = xs;
    ev.ys = ys;
    ev.delta = delta;
    ev.ytar = ytar;
    events.push_back(ev);
    hYpFpYFp->Fill(ypfp, yfp);
    hYsXsBase->Fill(ys, xs);
    nBase++;
  }

  cout << "Selected central FP/PFP events: " << nBase << endl;
  if (events.size() < 100) {
    cerr << "ERROR: too few events for automatic FP/PFP selection." << endl;
    return;
  }

  double mx, my, vx, vy, ux, uy;
  ComputePCAAxes(events, mx, my, vx, vy, ux, uy);
  AssignSU(events, mx, my, vx, vy, ux, uy);

  vector<double> sVals;
  sVals.reserve(events.size());
  for (const auto &e : events) sVals.push_back(e.s);
  double sMin = Quantile(sVals, 0.005);
  double sMax = Quantile(sVals, 0.995);
  double sPad = 0.05 * (sMax - sMin);
  sMin -= sPad;
  sMax += sPad;

  TH1D *hS = new TH1D("hS_ridge_coordinate",
                      Form("Run %d foil %d central #delta;ridge coordinate s;counts", nrun, foilIndex),
                      sBins, sMin, sMax);
  for (const auto &e : events) hS->Fill(e.s);

  TH1D *hSSmooth = (TH1D*)hS->Clone("hS_ridge_coordinate_smooth");
  for (int i = 0; i < smoothPasses; i++) hSSmooth->Smooth(1);

  vector<Peak1D> peaks = FindPeaks1D(hSSmooth, minPeakFrac, minPeakSeparationS);
  cout << "Detected 1D ridge peaks: " << peaks.size() << endl;

  vector<double> predYpFp, predYFp, predS;
  PredictYFpYpFpGuides(info, foilIndex, ndCentral, mx, my, vx, vy,
                       predYpFp, predYFp, predS);

  AssignPeaksToYcols(peaks, predS, maxPredDistanceS);

  // Remove unassigned peaks.
  vector<Peak1D> assigned;
  for (auto &p : peaks) {
    if (p.yscol >= 0) assigned.push_back(p);
  }
  peaks = assigned;

  cout << "Assigned peaks:" << endl;
  for (auto &p : peaks) {
    cout << "  yscol=" << p.yscol
         << " peak_s=" << p.s
         << " pred_s=" << p.sPred
         << " height=" << p.height
         << " frac=" << p.frac << endl;
  }

  if (peaks.empty()) {
    cerr << "ERROR: no peaks could be assigned to yscol indices." << endl;
    return;
  }

  // Build s boundaries from neighboring assigned peaks.
  sort(peaks.begin(), peaks.end(), [](const Peak1D &a, const Peak1D &b) { return a.s < b.s; });
  for (int i = 0; i < (int)peaks.size(); i++) {
    double lo = peaks[i].s - sHalfWidthFallback;
    double hi = peaks[i].s + sHalfWidthFallback;

    if (i > 0) lo = 0.5 * (peaks[i-1].s + peaks[i].s);
    if (i + 1 < (int)peaks.size()) hi = 0.5 * (peaks[i].s + peaks[i+1].s);

    double width = hi - lo;
    double center = peaks[i].s;
    double half = 0.5 * width * sGuardScale;
    peaks[i].sLo = center - half;
    peaks[i].sHi = center + half;

    vector<double> localU;
    for (const auto &e : events) {
      if (e.s >= peaks[i].sLo && e.s < peaks[i].sHi) localU.push_back(e.u);
    }

    if ((int)localU.size() < minEventsPerCut) {
      cout << "Skipping yscol=" << peaks[i].yscol
           << " because only " << localU.size() << " events in local s window." << endl;
      peaks[i].cut = nullptr;
      continue;
    }

    double qlo = Quantile(localU, uQuantileLow);
    double qhi = Quantile(localU, uQuantileHigh);
    double upad = std::max(uPadMin, uPadFrac * (qhi - qlo));
    peaks[i].uLo = qlo - upad;
    peaks[i].uHi = qhi + upad;
    peaks[i].nEvents = localU.size();

    TString cname = Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",
                         peaks[i].yscol, foilIndex, ndCentral);
    peaks[i].cut = MakeRotatedRectangleCut(cname, mx, my, vx, vy, ux, uy,
                                           peaks[i].sLo, peaks[i].sHi,
                                           peaks[i].uLo, peaks[i].uHi);
    peaks[i].hYsXs = new TH2D(Form("hYsXs_auto_yscol_%d_nfoil_%d_ndel_%d",
                                   peaks[i].yscol, foilIndex, ndCentral),
                              Form("Run %d auto YpFp/YFp cut: foil %d nd %d yscol %d;Ys;Xs",
                                   nrun, foilIndex, ndCentral, peaks[i].yscol),
                              100, -7.0, 7.0, 120, -12.5, 12.5);
  }

  // Fill selected Xs/Ys histograms using cut polygons.
  for (const auto &e : events) {
    for (auto &p : peaks) {
      if (!p.cut || !p.hYsXs) continue;
      if (p.cut->IsInside(e.ypfp, e.yfp)) p.hYsXs->Fill(e.ys, e.xs);
    }
  }

  TString outCutFile = Form("cuts/YpFpYFp_%s_%s_cut.root", info.opticsID.Data(), inputID.Data());
  TString outPdf = Form("plots/auto_ypfp_yfp_central_run%d_%s.pdf", nrun, outTag.Data());
  TString outCsv = Form("plots/auto_ypfp_yfp_central_run%d_%s.csv", nrun, outTag.Data());

  gSystem->mkdir("cuts", kTRUE);
  gSystem->mkdir("plots", kTRUE);

  ofstream csv(outCsv.Data());
  csv << "run,tag,inputFileID,ytarCutTag,foilIndex,ndel,dlo,dhi,yscol,peak_s,pred_s,peak_height,peak_frac,sLo,sHi,uLo,uHi,nEvents\n";
  for (auto &p : peaks) {
    if (!p.cut) continue;
    csv << nrun << "," << outTag << "," << inputID << "," << ytarTag << ","
        << foilIndex << "," << ndCentral << "," << dLo << "," << dHi << ","
        << p.yscol << "," << p.s << "," << p.sPred << ","
        << p.height << "," << p.frac << ","
        << p.sLo << "," << p.sHi << "," << p.uLo << "," << p.uHi << ","
        << p.nEvents << "\n";
  }
  csv.close();

  TCanvas *c = new TCanvas("c", "auto ypfp yfp central cuts", 1100, 850);
  c->Print(outPdf + "[");

  c->Clear();
  gPad->SetLogz(1);
  hYpFpYFp->Draw("colz");
  for (int iy = 0; iy < 9; iy++) {
    TLine *ln = new TLine(predYpFp[iy] - 0.002, predYFp[iy],
                          predYpFp[iy] + 0.002, predYFp[iy]);
    ln->SetLineColor(kRed);
    ln->SetLineWidth(2);
    ln->Draw("same");
    TText *txt = new TText(predYpFp[iy], predYFp[iy], Form("%d", iy));
    txt->SetTextColor(kRed);
    txt->SetTextSize(0.03);
    txt->Draw("same");
  }
  for (auto &p : peaks) if (p.cut) p.cut->Draw("L same");
  TLatex tx;
  tx.SetNDC(); tx.SetTextSize(0.026);
  tx.DrawLatex(0.12, 0.94, Form("run %d foil %d central delta nd=%d [%.1f, %.1f]", nrun, foilIndex, ndCentral, dLo, dHi));
  c->Print(outPdf);
  gPad->SetLogz(0);

  c->Clear();
  hS->SetLineColor(kGray+2);
  hS->Draw("hist");
  hSSmooth->SetLineColor(kBlue);
  hSSmooth->SetLineWidth(3);
  hSSmooth->Draw("hist same");
  for (auto &p : peaks) {
    TLine *ln = new TLine(p.s, 0.0, p.s, hSSmooth->GetMaximum());
    ln->SetLineColor(kRed);
    ln->SetLineWidth(2);
    ln->Draw("same");
    tx.DrawLatexNDC(0.13, 0.88 - 0.035 * (&p - &peaks[0]),
                    Form("yscol %d: s=%.3f pred=%.3f", p.yscol, p.s, p.sPred));
  }
  c->Print(outPdf);

  c->Clear();
  gPad->SetLogz(1);
  hYsXsBase->Draw("colz");
  DrawYsGuideLines();
  tx.DrawLatexNDC(0.13, 0.94, "Base selected events after ytar + central delta; red = nominal Ys guide lines");
  c->Print(outPdf);
  gPad->SetLogz(0);

  for (auto &p : peaks) {
    if (!p.cut || !p.hYsXs) continue;
    c->Clear();
    gPad->SetLogz(1);
    p.hYsXs->Draw("colz");
    DrawYsGuideLines();
    tx.DrawLatexNDC(0.13, 0.94,
                    Form("yscol %d selected by auto YpFp/YFp cut | n=%ld", p.yscol, (long)p.hYsXs->GetEntries()));
    c->Print(outPdf);
    gPad->SetLogz(0);
  }

  c->Print(outPdf + "]");

  TFile *fout = new TFile(outCutFile, "UPDATE");
  hYpFpYFp->Write("hYpFpYFp_auto_base", TObject::kOverwrite);
  hYsXsBase->Write("hYsXs_auto_base", TObject::kOverwrite);
  hS->Write("hS_auto", TObject::kOverwrite);
  hSSmooth->Write("hS_auto_smooth", TObject::kOverwrite);

  for (auto &p : peaks) {
    if (p.cut) p.cut->Write("", TObject::kOverwrite);
    if (p.hYsXs) p.hYsXs->Write("", TObject::kOverwrite);
  }
  fout->Close();

  cout << "Wrote cut ROOT file: " << outCutFile << endl;
  cout << "Wrote PDF: " << outPdf << endl;
  cout << "Wrote CSV: " << outCsv << endl;

  fin->Close();
}
