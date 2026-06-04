// assign_yfp_ypfp_angleScanBands_qpStripe_with_pages.C
//
// Diagnostic: choose the projection angle that maximizes resolved 1D peak structure.
// This directly tests: "rotate until the yfp/ypfp islands separate best."
//
// Coordinates:
//   x = ypfp
//   y = yfp
//   xz = (x - mean_x)/sigma_x
//   yz = (y - mean_y)/sigma_y
//
// For each theta in [0,180):
//   q = xz*cos(theta) + yz*sin(theta)        // candidate separation coordinate
//   p = -xz*sin(theta) + yz*cos(theta)       // perpendicular coordinate
//   histogram q, smooth, find peaks, score
//
// The best theta is the one with the largest number of accepted q peaks.
// Tie-breakers prefer stronger and better-separated peaks.
//
// Run from repo top directory, e.g.
//   hcana -l -q 'assign_yfp_ypfp_angleScanBands_strongWeak_with_pages.C(1544,-10,-8,"auto_ycut")'
//
// Useful test:
//   hcana -l -q 'assign_yfp_ypfp_angleScanBands_strongWeak_with_pages.C(1544,-8,-5,"auto_ycut",9,1.0,0.18,0.05,0.15,2,0.005,0.04,0.45)'

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using std::cout;
using std::endl;

struct OpticsRunInfo {
  int run = -1;
  TString opticsID = "";
  double centAngle = 0.0;
  int numFoil = 0;
  int sieveFlag = 0;
  int ndelcut = 0;
  std::vector<double> zfoil;
  std::vector<double> delcut;
};

struct QPeak {
  double q = 0.0;
  double height = 0.0;
  double prominence = 0.0;
  int bin = 0;
};

struct AngleResult {
  double thetaDeg = 0.0;
  double thetaRad = 0.0;
  int nPeaks = 0;
  double score = 0.0;
  double totalPeakHeight = 0.0;
  double meanValleyDrop = 0.0;
  std::vector<QPeak> peaks;
  std::vector<double> bounds;
};

static std::vector<TString> SplitCSV_angleScan(const TString& line) {
  std::vector<TString> out;
  TString s(line);
  TObjArray* arr = s.Tokenize(",");
  for (int i = 0; i < arr->GetEntries(); ++i) {
    TString tok = ((TObjString*)arr->At(i))->GetString();
    tok = tok.Strip(TString::kBoth);
    out.push_back(tok);
  }
  delete arr;
  return out;
}

static bool ReadOpticsRunInfo_angleScan(int nrun, OpticsRunInfo& info,
                                         const char* metaFile="DATfiles/list_of_optics_run.dat") {
  std::ifstream fin(metaFile);
  if (!fin.is_open()) {
    cout << "ERROR: cannot open metadata file: " << metaFile << endl;
    return false;
  }

  std::string raw;
  while (std::getline(fin, raw)) {
    TString line(raw.c_str());
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0 || line.BeginsWith("#")) continue;

    auto tok = SplitCSV_angleScan(line);
    if (tok.size() < 6) continue;
    if (tok[0].Atoi() != nrun) continue;

    info.run       = tok[0].Atoi();
    info.opticsID  = tok[1];
    info.centAngle = tok[2].Atof();
    info.numFoil   = tok[3].Atoi();
    info.sieveFlag = tok[4].Atoi();
    info.ndelcut   = tok[5].Atoi();

    int idx = 6;
    for (int i = 0; i < info.numFoil && idx < (int)tok.size(); ++i, ++idx)
      info.zfoil.push_back(tok[idx].Atof());
    for (int i = 0; i < info.ndelcut + 1 && idx < (int)tok.size(); ++i, ++idx)
      info.delcut.push_back(tok[idx].Atof());

    return true;
  }

  cout << "ERROR: run " << nrun << " not found in " << metaFile << endl;
  return false;
}

static TCutG* GetFirstCutG_angleScan(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj && obj->InheritsFrom(TCutG::Class())) return (TCutG*)obj;
  }
  return nullptr;
}

static TCutG* LoadYtarCut_angleScan(int nrun, const TString& ytarTag, int foilIndex=0, int numFoil=1) {
  // Do not change the yfp/ypfp angle-scan logic here.
  // This helper only chooses the input ytar TCutG.

  std::vector<TString> fnames;

  // Christine's current ytar auto-cut files are often written with a blank tag:
  //   cuts/ytar_ridge_cut__run1540.root
  // Keep that as the first choice so the existing working input is preserved.
  fnames.push_back(Form("cuts/ytar_ridge_cut__run%d.root", nrun));

  // Also support the older tagged naming convention.
  if (ytarTag.Length() > 0) {
    TString tagged = Form("cuts/ytar_ridge_cut_%s_run%d.root", ytarTag.Data(), nrun);
    bool already=false;
    for (const auto& f0 : fnames) { if (f0 == tagged) already=true; }
    if (!already) fnames.push_back(tagged);
  }

  std::vector<TString> names;
  names.push_back(Form("delta_vs_ytar_cut_foil%d", foilIndex));
  names.push_back(Form("ytar_delta_cut_foil%d", foilIndex));
  names.push_back(Form("foil%d", foilIndex));
  names.push_back(Form("cut_foil%d", foilIndex));

  // Only allow generic/fallback cut names for single-foil files.
  // For multi-foil runs, silently grabbing a generic/first TCutG is dangerous.
  if (numFoil <= 1) {
    names.push_back("delta_vs_ytar_cut");
    names.push_back("ytar_delta_cut");
  }

  for (const auto& fname : fnames) {
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
      cout << "WARNING: could not open ytar cut file: " << fname << endl;
      continue;
    }

    for (auto& name : names) {
      TCutG* c = (TCutG*)f->Get(name);
      if (c) {
        cout << "Loaded ytar cut: " << fname << " :: " << name
             << "  (foilIndex=" << foilIndex << ")" << endl;
        return c;
      }
    }

    if (numFoil <= 1) {
      TCutG* first = GetFirstCutG_angleScan(f);
      if (first) {
        cout << "Loaded first available TCutG from " << fname << ": " << first->GetName() << endl;
        return first;
      }
    }

    cout << "WARNING: no matching TCutG found in " << fname
         << " for foilIndex=" << foilIndex << endl;
  }

  cout << "ERROR: no usable ytar TCutG found for run " << nrun
       << ", foilIndex=" << foilIndex << endl;
  return nullptr;
}

static double LocalValleyMin_angleScan(TH1D* h, int b1, int b2) {
  if (b2 < b1) std::swap(b1,b2);
  double valley = std::numeric_limits<double>::max();
  for (int b=b1; b<=b2; ++b)
    valley = std::min(valley, h->GetBinContent(b));
  return valley;
}

static std::vector<QPeak> FindPeaks1D_angleScan(TH1D* hSmooth,
                                                 int maxPeaks,
                                                 double minPeakSep,
                                                 double minPeakFraction,
                                                 double minProminenceFraction) {
  struct CandPeak {
    double q = 0.0;
    double height = 0.0;
    double prominence = 0.0;
    int bin = 0;
  };

  std::vector<CandPeak> localMaxima;
  const double maxContent = hSmooth->GetMaximum();
  const double minHeight = minPeakFraction * maxContent;

  // First collect all local maxima above height threshold.
  for (int ib=2; ib<hSmooth->GetNbinsX(); ++ib) {
    const double ym = hSmooth->GetBinContent(ib-1);
    const double y0 = hSmooth->GetBinContent(ib);
    const double yp = hSmooth->GetBinContent(ib+1);
    if (y0 > ym && y0 >= yp && y0 > minHeight) {
      localMaxima.push_back({hSmooth->GetBinCenter(ib), y0, 0.0, ib});
    }
  }

  // Estimate each maximum's valley prominence relative to neighboring maxima.
  // This rejects shoulder/tail ripples that are not separated by a real dip.
  std::sort(localMaxima.begin(), localMaxima.end(),
            [](const CandPeak& a, const CandPeak& b){ return a.q < b.q; });

  std::vector<QPeak> cands;
  for (size_t i=0; i<localMaxima.size(); ++i) {
    const double h0 = localMaxima[i].height;

    double leftProm = 1.0;
    double rightProm = 1.0;

    if (i > 0) {
      const double valley = LocalValleyMin_angleScan(hSmooth, localMaxima[i-1].bin, localMaxima[i].bin);
      leftProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    if (i+1 < localMaxima.size()) {
      const double valley = LocalValleyMin_angleScan(hSmooth, localMaxima[i].bin, localMaxima[i+1].bin);
      rightProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    // For interior peaks, require separation on both sides. For edge peaks, one side is enough.
    double prom = 0.0;
    if (i == 0 && localMaxima.size() > 1) prom = rightProm;
    else if (i+1 == localMaxima.size() && localMaxima.size() > 1) prom = leftProm;
    else if (localMaxima.size() == 1) prom = 1.0;
    else prom = std::min(leftProm, rightProm);

    if (prom >= minProminenceFraction) {
      cands.push_back({localMaxima[i].q, h0, prom, localMaxima[i].bin});
    }
  }

  // Keep strongest/prominent peaks first, then suppress duplicates too close in q.
  std::sort(cands.begin(), cands.end(),
            [](const QPeak& a, const QPeak& b){
              if (a.prominence != b.prominence) return a.prominence > b.prominence;
              return a.height > b.height;
            });

  std::vector<QPeak> peaks;
  for (const auto& cand : cands) {
    bool tooClose=false;
    for (const auto& pk : peaks) {
      if (std::fabs(cand.q - pk.q) < minPeakSep) { tooClose=true; break; }
    }
    if (tooClose) continue;
    peaks.push_back(cand);
    if ((int)peaks.size() >= maxPeaks) break;
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const QPeak& a, const QPeak& b){ return a.q < b.q; });

  return peaks;
}

static std::vector<double> BoundariesFromPeaks_angleScan(const std::vector<QPeak>& peaks) {
  std::vector<double> bounds;
  for (size_t i=0; i+1<peaks.size(); ++i)
    bounds.push_back(0.5*(peaks[i].q + peaks[i+1].q));
  return bounds;
}

static double MeanValleyDrop_angleScan(TH1D* hSmooth, const std::vector<QPeak>& peaks) {
  if (peaks.size() < 2) return 0.0;

  double sumDrop = 0.0;
  int n = 0;

  for (size_t i=0; i+1<peaks.size(); ++i) {
    int b1 = peaks[i].bin;
    int b2 = peaks[i+1].bin;
    if (b2 < b1) std::swap(b1,b2);

    double valley = std::numeric_limits<double>::max();
    for (int b=b1; b<=b2; ++b)
      valley = std::min(valley, hSmooth->GetBinContent(b));

    const double lowPeak = std::min(peaks[i].height, peaks[i+1].height);
    if (lowPeak > 0 && std::isfinite(valley)) {
      sumDrop += (lowPeak - valley) / lowPeak; // 0=no valley, 1=deep valley
      n++;
    }
  }

  return (n > 0) ? sumDrop / (double)n : 0.0;
}

static AngleResult ScoreAngle_angleScan(const std::vector<double>& xzvals,
                                         const std::vector<double>& yzvals,
                                         double thetaDeg,
                                         int maxBands,
                                         double minPeakSep,
                                         double minPeakFraction,
                                         double minProminenceFraction,
                                         int smoothPasses,
                                         int nBinsQ = 240) {
  const Long64_t N = (Long64_t)xzvals.size();
  const double th = thetaDeg * TMath::Pi()/180.0;
  const double c = std::cos(th);
  const double s = std::sin(th);

  std::vector<double> qvals;
  qvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) qvals.push_back(xzvals[i]*c + yzvals[i]*s);

  double qmin=*std::min_element(qvals.begin(), qvals.end());
  double qmax=*std::max_element(qvals.begin(), qvals.end());
  double qpad=0.05*(qmax-qmin);
  qmin-=qpad; qmax+=qpad;

  TH1D* hQ = new TH1D(Form("hQ_scan_tmp_%g", thetaDeg), "temporary q scan", nBinsQ, qmin, qmax);
  hQ->SetDirectory(nullptr);
  for (double q : qvals) hQ->Fill(q);

  TH1D* hS = (TH1D*)hQ->Clone(Form("hQ_scan_tmp_smooth_%g", thetaDeg));
  hS->SetDirectory(nullptr);
  for (int i=0; i<smoothPasses; ++i) hS->Smooth(1);

  AngleResult res;
  res.thetaDeg = thetaDeg;
  res.thetaRad = th;
  res.peaks = FindPeaks1D_angleScan(hS, maxBands, minPeakSep, minPeakFraction, minProminenceFraction);
  res.bounds = BoundariesFromPeaks_angleScan(res.peaks);
  res.nPeaks = (int)res.peaks.size();
  for (const auto& pk : res.peaks) res.totalPeakHeight += pk.height;
  res.meanValleyDrop = MeanValleyDrop_angleScan(hS, res.peaks);

  // Strongly prefer more peaks. Tie-break by valley depth and total height.
  // This is intentionally simple and transparent.
  double meanProm = 0.0;
  for (const auto& pk : res.peaks) meanProm += pk.prominence;
  if (!res.peaks.empty()) meanProm /= (double)res.peaks.size();

  res.score = 1000000.0 * res.nPeaks
            + 100000.0 * meanProm
            + 10000.0 * res.meanValleyDrop
            + 0.1 * res.totalPeakHeight;

  delete hQ;
  delete hS;
  return res;
}

static AngleResult ScoreStripePhi_angleScan(const std::vector<double>& qvals,
                                           const std::vector<double>& pvals,
                                           double phiDeg,
                                           int maxBands,
                                           double minPeakSep,
                                           double minPeakFraction,
                                           double minProminenceFraction,
                                           int smoothPasses,
                                           int nBinsR = 240) {
  // Second-stage scan in the already-computed q,p plane.
  // We look for the coordinate r that is perpendicular to the visible stripe family:
  //   u = q*cos(phi) - p*sin(phi)       // along-stripe coordinate
  //   r = q*sin(phi) + p*cos(phi)       // stripe-separation coordinate
  // Bands are then assigned in r, not in q.
  const Long64_t N = (Long64_t)qvals.size();
  const double ph = phiDeg * TMath::Pi()/180.0;
  const double c = std::cos(ph);
  const double s = std::sin(ph);

  std::vector<double> rvals;
  rvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) rvals.push_back(qvals[i]*s + pvals[i]*c);

  double rmin=*std::min_element(rvals.begin(), rvals.end());
  double rmax=*std::max_element(rvals.begin(), rvals.end());
  double rpad=0.05*(rmax-rmin);
  rmin-=rpad; rmax+=rpad;

  TH1D* hR = new TH1D(Form("hR_stripe_scan_tmp_%g", phiDeg), "temporary r stripe scan", nBinsR, rmin, rmax);
  hR->SetDirectory(nullptr);
  for (double r : rvals) hR->Fill(r);

  TH1D* hS = (TH1D*)hR->Clone(Form("hR_stripe_scan_tmp_smooth_%g", phiDeg));
  hS->SetDirectory(nullptr);
  for (int i=0; i<smoothPasses; ++i) hS->Smooth(1);

  AngleResult res;
  res.thetaDeg = phiDeg;
  res.thetaRad = ph;
  res.peaks = FindPeaks1D_angleScan(hS, maxBands, minPeakSep, minPeakFraction, minProminenceFraction);
  res.bounds = BoundariesFromPeaks_angleScan(res.peaks);
  res.nPeaks = (int)res.peaks.size();
  for (const auto& pk : res.peaks) res.totalPeakHeight += pk.height;
  res.meanValleyDrop = MeanValleyDrop_angleScan(hS, res.peaks);

  double meanProm = 0.0;
  for (const auto& pk : res.peaks) meanProm += pk.prominence;
  if (!res.peaks.empty()) meanProm /= (double)res.peaks.size();

  res.score = 1000000.0 * res.nPeaks
            + 100000.0 * meanProm
            + 10000.0 * res.meanValleyDrop
            + 0.1 * res.totalPeakHeight;

  delete hR;
  delete hS;
  return res;
}

void assign_yfp_ypfp_angleScanBands_qpStripe_with_pages(Int_t nrun=1544,
                                                Double_t deltaMin=-10.0,
                                                Double_t deltaMax=-8.0,
                                                TString ytarTag="auto_ycut",
                                                Int_t maxBands=9,
                                                Double_t thetaStepDeg=1.0,
                                                Double_t minPeakSep=0.18,
                                                Double_t minPeakFraction=0.05,
                                                Double_t minProminenceFraction=0.15,
                                                Int_t smoothPasses=2,
                                                Double_t minWeakPeakFraction=0.005,
                                                Double_t minWeakProminenceFraction=0.04,
                                                Double_t weakMinSepFactor=0.45,
                                                Int_t foilIndex=0,
                                                Bool_t useYtarCut=true,
                                                Long64_t maxEvents=-1) {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  OpticsRunInfo info;
  if (!ReadOpticsRunInfo_angleScan(nrun, info)) return;

  if (info.numFoil <= 0) {
    cout << "ERROR: metadata says NumFoil <= 0 for run " << nrun << endl;
    return;
  }
  if (foilIndex < 0 || foilIndex >= info.numFoil) {
    cout << "ERROR: requested foilIndex=" << foilIndex
         << " but metadata NumFoil=" << info.numFoil
         << " for run " << nrun << endl;
    return;
  }

  TString inroot = Form("ROOTfiles/OPTICS/nps_hms_optics_%d_1_-1.root", nrun);
  TFile* fin = TFile::Open(inroot, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file: " << inroot << endl;
    return;
  }

  TTree* T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: tree T not found in " << inroot << endl;
    return;
  }

  TCutG* ytarCut = nullptr;
  if (useYtarCut) {
    ytarCut = LoadYtarCut_angleScan(nrun, ytarTag, foilIndex, info.numFoil);
    if (!ytarCut) {
      cout << "ERROR: useYtarCut=true but no ytar cut was loaded. Stopping." << endl;
      return;
    }
  }

  Double_t sumnpe=0, etracknorm=0;
  Double_t ytar=0, delta=0, yfp=0, ypfp=0;

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("H.cer.npeSum",1);
  T->SetBranchStatus("H.cal.etottracknorm",1);
  T->SetBranchStatus("H.gtr.y",1);
  T->SetBranchStatus("H.gtr.dp",1);
  T->SetBranchStatus("H.dc.y_fp",1);
  T->SetBranchStatus("H.dc.yp_fp",1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.cal.etottracknorm", &etracknorm);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);

  // Plot convention: x = ypfp, y = yfp.
  std::vector<double> xvals; // ypfp
  std::vector<double> yvals; // yfp

  Long64_t nentries = T->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

  Long64_t nPassBasic=0, nPassYtar=0, nPassDelta=0;
  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    if (!(sumnpe > 6.0 && etracknorm > 0.65)) continue;
    nPassBasic++;

    if (useYtarCut && ytarCut) {
      if (!ytarCut->IsInside(ytar, delta)) continue; // x=ytar, y=delta
    }
    nPassYtar++;

    if (!(delta >= deltaMin && delta < deltaMax)) continue;
    nPassDelta++;

    xvals.push_back(ypfp);
    yvals.push_back(yfp);
  }

  const Long64_t N = xvals.size();

  cout << "\n=== Angle-scan yfp/ypfp band diagnostic ===" << endl;
  cout << "Run: " << nrun << endl;
  cout << "Input: " << inroot << endl;
  cout << "Metadata opticsID: " << info.opticsID << ", NumFoil: " << info.numFoil << endl;
  cout << "Requested foilIndex: " << foilIndex;
  if (foilIndex >= 0 && foilIndex < (int)info.zfoil.size()) cout << "  zfoil=" << info.zfoil[foilIndex] << " cm";
  cout << endl;
  cout << "Delta slice: [" << deltaMin << ", " << deltaMax << ") %" << endl;
  cout << "Max allowed bands: " << maxBands << endl;
  cout << "Theta step: " << thetaStepDeg << " deg" << endl;
  cout << "Min projected peak separation: " << minPeakSep << endl;
  cout << "Min projected peak height fraction: " << minPeakFraction << endl;
  cout << "Min projected peak prominence fraction: " << minProminenceFraction << endl;
  cout << "Smoothing passes: " << smoothPasses << endl;
  cout << "Weak peak fraction after angle choice: " << minWeakPeakFraction << endl;
  cout << "Weak peak prominence after angle choice: " << minWeakProminenceFraction << endl;
  cout << "Weak peak min separation factor after angle choice: " << weakMinSepFactor << endl;
  cout << "Events after PID/basic: " << nPassBasic << endl;
  cout << "Events after ytar cut:  " << nPassYtar << endl;
  cout << "Events in delta slice: " << nPassDelta << endl;
  cout << "Events used:           " << N << endl;

  if (N < 100) {
    cout << "ERROR: too few events." << endl;
    return;
  }

  double mx=0, my=0;
  for (Long64_t i=0; i<N; ++i) { mx += xvals[i]; my += yvals[i]; }
  mx /= (double)N;
  my /= (double)N;

  double raw_sxx=0, raw_syy=0, raw_sxy=0;
  for (Long64_t i=0; i<N; ++i) {
    const double dx = xvals[i] - mx;
    const double dy = yvals[i] - my;
    raw_sxx += dx*dx;
    raw_syy += dy*dy;
    raw_sxy += dx*dy;
  }
  raw_sxx /= (double)(N-1);
  raw_syy /= (double)(N-1);
  raw_sxy /= (double)(N-1);

  const double sx = std::sqrt(raw_sxx);
  const double sy = std::sqrt(raw_syy);
  if (sx <= 0 || sy <= 0) {
    cout << "ERROR: zero RMS in one coordinate." << endl;
    return;
  }

  std::vector<double> xzvals, yzvals;
  xzvals.reserve(N); yzvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    xzvals.push_back((xvals[i] - mx)/sx);
    yzvals.push_back((yvals[i] - my)/sy);
  }

  // Scan all projection angles.
  std::vector<AngleResult> scan;
  for (double th=0.0; th<180.0; th += thetaStepDeg) {
    scan.push_back(ScoreAngle_angleScan(xzvals, yzvals, th,
                                         maxBands, minPeakSep,
                                         minPeakFraction, minProminenceFraction,
                                         smoothPasses));
  }

  if (scan.empty()) {
    cout << "ERROR: no scan results." << endl;
    return;
  }

  auto bestIt = std::max_element(scan.begin(), scan.end(),
    [](const AngleResult& a, const AngleResult& b){
      if (a.score != b.score) return a.score < b.score;
      return a.nPeaks < b.nPeaks;
    });
  AngleResult best = *bestIt;

  cout << "\n=== Best angle-scan result ===" << endl;
  cout << "Best theta: " << best.thetaDeg << " deg in standardized xz-yz space" << endl;
  cout << "Accepted peaks: " << best.nPeaks << endl;
  cout << "Score: " << best.score << endl;
  cout << "Mean valley drop: " << best.meanValleyDrop << endl;
  cout << "Total peak height: " << best.totalPeakHeight << endl;
  for (size_t i=0; i<best.peaks.size(); ++i) {
    cout << "  band " << i << ": r peak=" << best.peaks[i].q
         << ", height=" << best.peaks[i].height
         << ", prominence=" << best.peaks[i].prominence << endl;
  }
  cout << "Boundaries:";
  for (double b : best.bounds) cout << " " << b;
  cout << endl;

  // Compute final q,p coordinates using best theta.
  const double c = std::cos(best.thetaRad);
  const double s = std::sin(best.thetaRad);
  std::vector<double> qvals, pvals;
  qvals.reserve(N); pvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    const double q =  xzvals[i]*c + yzvals[i]*s;
    const double p = -xzvals[i]*s + yzvals[i]*c;
    qvals.push_back(q);
    pvals.push_back(p);
  }

  // -------------------------------------------------------------------------
  // Second-stage stripe scan.
  // The first scan finds a useful q,p frame, but the visible islands can still
  // appear as tilted parallel stripes in that frame.  This scan rotates only
  // within q,p and assigns bands using r, the coordinate perpendicular to the
  // stripe direction.
  // -------------------------------------------------------------------------
  std::vector<AngleResult> stripeScan;
  for (double ph=-45.0; ph<=45.0 + 1e-9; ph += thetaStepDeg) {
    stripeScan.push_back(ScoreStripePhi_angleScan(qvals, pvals, ph,
                                                  maxBands, minPeakSep,
                                                  minPeakFraction, minProminenceFraction,
                                                  smoothPasses));
  }

  auto stripeBestIt = std::max_element(stripeScan.begin(), stripeScan.end(),
    [](const AngleResult& a, const AngleResult& b){
      if (a.score != b.score) return a.score < b.score;
      return a.nPeaks < b.nPeaks;
    });
  AngleResult stripeBest = *stripeBestIt;

  cout << "\n=== Best q-p stripe result ===" << endl;
  cout << "Best phi: " << stripeBest.thetaDeg << " deg in q-p space" << endl;
  cout << "Accepted r peaks: " << stripeBest.nPeaks << endl;
  cout << "Score: " << stripeBest.score << endl;
  cout << "Mean valley drop: " << stripeBest.meanValleyDrop << endl;
  cout << "Total peak height: " << stripeBest.totalPeakHeight << endl;
  for (size_t i=0; i<stripeBest.peaks.size(); ++i) {
    cout << "  stripe band " << i << ": r peak=" << stripeBest.peaks[i].q
         << ", height=" << stripeBest.peaks[i].height
         << ", prominence=" << stripeBest.peaks[i].prominence << endl;
  }

  const double cp = std::cos(stripeBest.thetaRad);
  const double sp = std::sin(stripeBest.thetaRad);
  std::vector<double> uvals, rvals;
  uvals.reserve(N); rvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    const double u = qvals[i]*cp - pvals[i]*sp;
    const double r = qvals[i]*sp + pvals[i]*cp;
    uvals.push_back(u);
    rvals.push_back(r);
  }

  double xmin=*std::min_element(xvals.begin(), xvals.end());
  double xmax=*std::max_element(xvals.begin(), xvals.end());
  double ymin=*std::min_element(yvals.begin(), yvals.end());
  double ymax=*std::max_element(yvals.begin(), yvals.end());
  double xpad=0.10*(xmax-xmin), ypad=0.10*(ymax-ymin);
  xmin-=xpad; xmax+=xpad; ymin-=ypad; ymax+=ypad;

  double qmin=*std::min_element(qvals.begin(), qvals.end());
  double qmax=*std::max_element(qvals.begin(), qvals.end());
  double pmin=*std::min_element(pvals.begin(), pvals.end());
  double pmax=*std::max_element(pvals.begin(), pvals.end());
  double qpad=0.10*(qmax-qmin), ppad=0.10*(pmax-pmin);
  qmin-=qpad; qmax+=qpad; pmin-=ppad; pmax+=ppad;

  double umin=*std::min_element(uvals.begin(), uvals.end());
  double umax=*std::max_element(uvals.begin(), uvals.end());
  double rmin=*std::min_element(rvals.begin(), rvals.end());
  double rmax=*std::max_element(rvals.begin(), rvals.end());
  double upad=0.10*(umax-umin), rpad=0.10*(rmax-rmin);
  umin-=upad; umax+=upad; rmin-=rpad; rmax+=rpad;

  TH2D* hOrig = new TH2D("hOrig", Form("Run %d, %.1f < #delta < %.1f;ypfp;yfp", nrun, deltaMin, deltaMax),
                         220, xmin, xmax, 220, ymin, ymax);
  TH2D* hQP = new TH2D("hQP", Form("First-stage angle-scan coords, Run %d;q;p", nrun),
                       220, qmin, qmax, 220, pmin, pmax);
  TH2D* hUR = new TH2D("hUR", Form("Second-stage stripe coords, Run %d;u = along-stripe coordinate;r = stripe-separation coordinate", nrun),
                       220, umin, umax, 220, rmin, rmax);
  TH1D* hQ = new TH1D("hR", Form("Run %d;r = stripe-separation coordinate;counts", nrun), 240, rmin, rmax);
  TH1D* hP = new TH1D("hU", Form("Run %d;u = along-stripe coordinate;counts", nrun), 220, umin, umax);

  for (Long64_t i=0; i<N; ++i) {
    hOrig->Fill(xvals[i], yvals[i]);
    hQP->Fill(qvals[i], pvals[i]);
    hUR->Fill(uvals[i], rvals[i]);
    hQ->Fill(rvals[i]);
    hP->Fill(uvals[i]);
  }

  TH1D* hQSmooth = (TH1D*)hQ->Clone("hQSmooth");
  for (int i=0; i<smoothPasses; ++i) hQSmooth->Smooth(1);

  // -------------------------------------------------------------------------
  // Final band-center selection:
  //   * The angle scan above uses ONLY strong/prominent peaks. This protects
  //     the rotation from tiny edge ripples.
  //   * After theta is fixed, we allow weak-but-separated peaks to become
  //     their own diagnostic bands. These are tagged as weak so a later stage
  //     can keep, reject, or down-weight them.
  // -------------------------------------------------------------------------
  std::vector<QPeak> bandPeaks = stripeBest.peaks;
  std::vector<int> bandIsWeak(bandPeaks.size(), 0);

  struct LocalQMax {
    double q = 0.0;
    double height = 0.0;
    double prominence = 0.0;
    int bin = 0;
  };

  std::vector<LocalQMax> localMaximaFinal;
  const double finalMaxContent = hQSmooth->GetMaximum();
  const double weakMinHeight = minWeakPeakFraction * finalMaxContent;

  for (int ib=2; ib<hQSmooth->GetNbinsX(); ++ib) {
    const double ym = hQSmooth->GetBinContent(ib-1);
    const double y0 = hQSmooth->GetBinContent(ib);
    const double yp = hQSmooth->GetBinContent(ib+1);
    if (y0 > ym && y0 >= yp && y0 > weakMinHeight) {
      localMaximaFinal.push_back({hQSmooth->GetBinCenter(ib), y0, 0.0, ib});
    }
  }

  std::sort(localMaximaFinal.begin(), localMaximaFinal.end(),
            [](const LocalQMax& a, const LocalQMax& b){ return a.q < b.q; });

  std::vector<QPeak> weakCandidates;
  for (size_t i=0; i<localMaximaFinal.size(); ++i) {
    const double h0 = localMaximaFinal[i].height;
    double leftProm = 1.0;
    double rightProm = 1.0;

    if (i > 0) {
      const double valley = LocalValleyMin_angleScan(hQSmooth, localMaximaFinal[i-1].bin, localMaximaFinal[i].bin);
      leftProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    if (i+1 < localMaximaFinal.size()) {
      const double valley = LocalValleyMin_angleScan(hQSmooth, localMaximaFinal[i].bin, localMaximaFinal[i+1].bin);
      rightProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    double prom = 0.0;
    if (i == 0 && localMaximaFinal.size() > 1) prom = rightProm;
    else if (i+1 == localMaximaFinal.size() && localMaximaFinal.size() > 1) prom = leftProm;
    else if (localMaximaFinal.size() == 1) prom = 1.0;
    else prom = std::min(leftProm, rightProm);

    if (prom < minWeakProminenceFraction) continue;

    QPeak pk{localMaximaFinal[i].q, h0, prom, localMaximaFinal[i].bin};

    // Skip if this weak candidate is already represented by a strong peak.
    bool alreadyStrong = false;
    for (const auto& strong : stripeBest.peaks) {
      if (std::fabs(pk.q - strong.q) < weakMinSepFactor * minPeakSep) {
        alreadyStrong = true;
        break;
      }
    }
    if (!alreadyStrong) weakCandidates.push_back(pk);
  }

  // Add weak candidates by prominence/height, without allowing tiny duplicate ripples.
  std::sort(weakCandidates.begin(), weakCandidates.end(),
            [](const QPeak& a, const QPeak& b){
              if (a.prominence != b.prominence) return a.prominence > b.prominence;
              return a.height > b.height;
            });

  int nWeakAdded = 0;
  for (const auto& cand : weakCandidates) {
    if ((int)bandPeaks.size() >= maxBands) break;

    bool tooClose = false;
    for (const auto& pk : bandPeaks) {
      if (std::fabs(cand.q - pk.q) < weakMinSepFactor * minPeakSep) {
        tooClose = true;
        break;
      }
    }
    if (tooClose) continue;

    bandPeaks.push_back(cand);
    bandIsWeak.push_back(1);
    nWeakAdded++;
  }

  // Sort final peak list and carry weak/strong labels along with it.
  std::vector<size_t> order(bandPeaks.size());
  for (size_t i=0; i<order.size(); ++i) order[i]=i;
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b){ return bandPeaks[a].q < bandPeaks[b].q; });

  std::vector<QPeak> sortedBandPeaks;
  std::vector<int> sortedBandIsWeak;
  for (size_t idx : order) {
    sortedBandPeaks.push_back(bandPeaks[idx]);
    sortedBandIsWeak.push_back(bandIsWeak[idx]);
  }
  bandPeaks.swap(sortedBandPeaks);
  bandIsWeak.swap(sortedBandIsWeak);

  std::vector<double> bandBounds = BoundariesFromPeaks_angleScan(bandPeaks);

  cout << "\n=== Final strong+weak q-band centers ===" << endl;
  cout << "Strong peaks from q-p stripe scan: " << stripeBest.peaks.size() << endl;
  cout << "Weak candidates considered after fixed angle: " << weakCandidates.size() << endl;
  cout << "Weak peaks added after fixed angle: " << nWeakAdded << endl;
  cout << "Final diagnostic bands: " << bandPeaks.size() << endl;
  for (size_t i=0; i<bandPeaks.size(); ++i) {
    cout << "  rband " << i
         << " [" << (bandIsWeak[i] ? "weak" : "strong") << "]"
         << ": peak r=" << bandPeaks[i].q
         << ", height=" << bandPeaks[i].height
         << ", prominence=" << bandPeaks[i].prominence << endl;
  }

  std::vector<int> bandIndex(N, -1);
  std::vector<int> bandCounts(bandPeaks.size(), 0);
  for (Long64_t i=0; i<N; ++i) {
    int b=0;
    while (b < (int)bandBounds.size() && rvals[i] > bandBounds[b]) b++;
    if (b >= 0 && b < (int)bandPeaks.size()) {
      bandIndex[i]=b;
      bandCounts[b]++;
    }
  }

  cout << "\n=== Final q-band assignment ===" << endl;
  for (size_t i=0; i<bandPeaks.size(); ++i) {
    cout << "  rband " << i << " [" << (bandIsWeak[i] ? "weak" : "strong") << "]: peak r=" << bandPeaks[i].q
         << ", assigned events=" << bandCounts[i] << endl;
  }

  TString outbase = Form("plots/yfp_ypfp_qpStripeStrongWeakBand_pages_run%d_foil%d_delta_%g_to_%g", nrun, foilIndex, deltaMin, deltaMax);
  outbase.ReplaceAll("-", "m");
  outbase.ReplaceAll(".", "p");
  TString outpdf = outbase + ".pdf";
  TString outroot = outbase + ".root";
  gSystem->mkdir("plots", kTRUE);

  std::vector<int> colors = {kRed+1, kOrange+7, kSpring+5, kGreen+2, kAzure+7, kBlue+1, kMagenta+1, kViolet+7, kCyan+2, kGray+2};
  std::vector<TGraph*> graphs(bandPeaks.size(), nullptr);
  for (size_t b=0; b<bandPeaks.size(); ++b) graphs[b] = new TGraph();

  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)graphs.size()) continue;
    int pnt = graphs[b]->GetN();
    graphs[b]->SetPoint(pnt, xvals[i], yvals[i]);
  }

  std::vector<double> bandLow(bandPeaks.size(), qmin), bandHigh(bandPeaks.size(), qmax);
  for (size_t b=0; b<bandPeaks.size(); ++b) {
    bandLow[b]  = (b==0) ? qmin : bandBounds[b-1];
    bandHigh[b] = (b+1==bandPeaks.size()) ? qmax : bandBounds[b];
  }

  auto qpFromUR = [&](double u0, double r0, double& q0, double& p0) {
    q0 = u0*cp + r0*sp;
    p0 = -u0*sp + r0*cp;
  };

  auto rawFromQP = [&](double q0, double p0, double& xraw, double& yraw) {
    const double xz0 = q0*c - p0*s;
    const double yz0 = q0*s + p0*c;
    xraw = mx + sx*xz0;
    yraw = my + sy*yz0;
  };

  auto drawConstRLine = [&](double r0, int color, int style, int width) {
    // Constant-r line: this is the boundary/center line for stripe assignment.
    const double uA = umin;
    const double uB = umax;
    double qA, pA, qB, pB;
    qpFromUR(uA, r0, qA, pA);
    qpFromUR(uB, r0, qB, pB);
    double xA, yA, xB, yB;
    rawFromQP(qA, pA, xA, yA);
    rawFromQP(qB, pB, xB, yB);
    TLine* line = new TLine(xA, yA, xB, yB);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("same");
  };

  auto drawConstRLineQP = [&](double r0, int color, int style, int width) {
    const double uA = umin;
    const double uB = umax;
    double qA, pA, qB, pB;
    qpFromUR(uA, r0, qA, pA);
    qpFromUR(uB, r0, qB, pB);
    TLine* line = new TLine(qA, pA, qB, pB);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("same");
  };

  auto drawConstRLineUR = [&](double r0, int color, int style, int width) {
    TLine* line = new TLine(umin, r0, umax, r0);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("same");
  };

  TCanvas* c1 = new TCanvas("c_angleScanBand_assign", "Angle-scan q-band assignment", 1200, 950);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetLogz();
  hOrig->Draw("COLZ");
  for (const auto& pk : bandPeaks) drawConstRLine(pk.q, kRed+1, 1, 2);
  for (double b : bandBounds) drawConstRLine(b, kRed+1, 2, 2);
  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.030);
  lat.DrawLatex(0.12,0.92,Form("Best #theta = %.2f deg z-space; bands=%zu, weak=%d", best.thetaDeg, bandPeaks.size(), nWeakAdded));

  c1->cd(2);
  gPad->SetLogz();
  hQP->Draw("COLZ");
  for (const auto& pk : bandPeaks) drawConstRLineQP(pk.q, kRed+1, 1, 2);
  for (double b : bandBounds) drawConstRLineQP(b, kRed+1, 2, 2);
  lat.DrawLatex(0.12,0.92,Form("q-p stripe correction #phi = %.2f deg; r-bands=%zu", stripeBest.thetaDeg, bandPeaks.size()));

  c1->cd(3);
  hQ->Draw("HIST");
  hQSmooth->SetLineColor(kRed+1);
  hQSmooth->SetLineWidth(2);
  hQSmooth->Draw("HIST SAME");
  for (const auto& pk : bandPeaks) {
    TLine* lpk = new TLine(pk.q, 0, pk.q, hQ->GetMaximum()*1.05);
    lpk->SetLineColor(kRed+1); lpk->SetLineWidth(2); lpk->Draw("same");
  }
  for (double b : bandBounds) {
    TLine* lb = new TLine(b, 0, b, hQ->GetMaximum()*1.05);
    lb->SetLineColor(kRed+1); lb->SetLineStyle(2); lb->Draw("same");
  }

  c1->cd(4);
  TH2D* hFrame = new TH2D("hFrameAngleBands", Form("angle-scan q-band assignments, Run %d;ypfp;yfp", nrun), 10, xmin, xmax, 10, ymin, ymax);
  hFrame->SetMinimum(0); hFrame->SetMaximum(1);
  hFrame->Draw("AXIS");
  for (size_t b=0; b<graphs.size(); ++b) {
    graphs[b]->SetMarkerStyle(20);
    graphs[b]->SetMarkerSize(0.18);
    graphs[b]->SetMarkerColor(colors[b % colors.size()]);
    graphs[b]->Draw("P SAME");
  }

  // Score-vs-angle page
  TGraph* gScore = new TGraph();
  TGraph* gNPeaks = new TGraph();
  TGraph* gValley = new TGraph();
  for (size_t i=0; i<scan.size(); ++i) {
    gScore->SetPoint(i, scan[i].thetaDeg, scan[i].score);
    gNPeaks->SetPoint(i, scan[i].thetaDeg, scan[i].nPeaks);
    gValley->SetPoint(i, scan[i].thetaDeg, scan[i].meanValleyDrop);
  }

  c1->SaveAs(outpdf + "[");
  c1->SaveAs(outpdf);

  TCanvas* cScore = new TCanvas("c_angle_scan_score", "Angle scan score", 1200, 900);
  cScore->Divide(1,3);

  cScore->cd(1);
  gScore->SetTitle(Form("Angle scan score, Run %d, %.1f<#delta<%.1f;theta (deg);score", nrun, deltaMin, deltaMax));
  gScore->SetLineWidth(2);
  gScore->Draw("AL");
  TLine* bestScoreLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestScoreLine->SetLineColor(kRed+1); bestScoreLine->SetLineWidth(2); bestScoreLine->Draw("same");

  cScore->cd(2);
  gNPeaks->SetTitle("Accepted q peaks vs angle;theta (deg);accepted peaks");
  gNPeaks->SetLineWidth(2);
  gNPeaks->Draw("AL");
  TLine* bestPeakLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestPeakLine->SetLineColor(kRed+1); bestPeakLine->SetLineWidth(2); bestPeakLine->Draw("same");

  cScore->cd(3);
  gValley->SetTitle("Mean valley drop vs angle;theta (deg);mean valley drop");
  gValley->SetLineWidth(2);
  gValley->Draw("AL");
  TLine* bestValleyLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestValleyLine->SetLineColor(kRed+1); bestValleyLine->SetLineWidth(2); bestValleyLine->Draw("same");

  cScore->SaveAs(outpdf);

  std::vector<TH2D*> hBandOrig(bandPeaks.size(), nullptr);
  std::vector<TH2D*> hBandQP(bandPeaks.size(), nullptr);
  std::vector<TH1D*> hBandQ(bandPeaks.size(), nullptr);

  for (size_t b=0; b<bandPeaks.size(); ++b) {
    hBandOrig[b] = new TH2D(Form("hBandOrig_%zu", b), Form("Band %zu only in original coords, Run %d;ypfp;yfp", b, nrun), 220, xmin, xmax, 220, ymin, ymax);
    hBandQP[b] = new TH2D(Form("hBandQP_%zu", b), Form("Band %zu only in q-p coords, Run %d;q;p", b, nrun), 220, qmin, qmax, 220, pmin, pmax);
    hBandQ[b] = new TH1D(Form("hBandR_%zu", b), Form("Band %zu r distribution, Run %d;r;counts", b, nrun), 240, rmin, rmax);
  }

  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)bandPeaks.size()) continue;
    hBandOrig[b]->Fill(xvals[i], yvals[i]);
    hBandQP[b]->Fill(qvals[i], pvals[i]);
    hBandQ[b]->Fill(rvals[i]);
  }

  for (size_t b=0; b<bandPeaks.size(); ++b) {
    TCanvas* cb = new TCanvas(Form("c_angle_band_%zu", b), Form("Angle-scan band %zu diagnostics", b), 1200, 950);
    cb->Divide(2,2);

    cb->cd(1);
    gPad->SetLogz();
    hOrig->Draw("COLZ");
    TGraph* gall = new TGraph();
    TGraph* gsel = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int pnt = gall->GetN();
      gall->SetPoint(pnt, xvals[i], yvals[i]);
      if (bandIndex[i] == (int)b) {
        int qn = gsel->GetN();
        gsel->SetPoint(qn, xvals[i], yvals[i]);
      }
    }
    gall->SetMarkerStyle(20); gall->SetMarkerSize(0.10); gall->SetMarkerColor(kGray+1); gall->Draw("P SAME");
    gsel->SetMarkerStyle(20); gsel->SetMarkerSize(0.22); gsel->SetMarkerColor(colors[b % colors.size()]); gsel->Draw("P SAME");
    for (const auto& pk : bandPeaks) drawConstRLine(pk.q, kRed+1, 2, 1);
    drawConstRLine(bandPeaks[b].q, colors[b % colors.size()], 1, 3);
    drawConstRLine(bandLow[b], colors[b % colors.size()], 2, 3);
    drawConstRLine(bandHigh[b], colors[b % colors.size()], 2, 3);
    lat.DrawLatex(0.12,0.92,Form("Angle band %zu [%s] highlighted in original coords", b, bandIsWeak[b] ? "weak" : "strong"));
    lat.DrawLatex(0.12,0.88,Form("peak r=%.3f, range [%.3f, %.3f], N=%d", bandPeaks[b].q, bandLow[b], bandHigh[b], bandCounts[b]));

    cb->cd(2);
    gPad->SetLogz();
    hQP->Draw("COLZ");
    TGraph* gallqp = new TGraph();
    TGraph* gselqp = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int pnt = gallqp->GetN();
      gallqp->SetPoint(pnt, qvals[i], pvals[i]);
      if (bandIndex[i] == (int)b) {
        int qn = gselqp->GetN();
        gselqp->SetPoint(qn, qvals[i], pvals[i]);
      }
    }
    gallqp->SetMarkerStyle(20); gallqp->SetMarkerSize(0.10); gallqp->SetMarkerColor(kGray+1); gallqp->Draw("P SAME");
    gselqp->SetMarkerStyle(20); gselqp->SetMarkerSize(0.22); gselqp->SetMarkerColor(colors[b % colors.size()]); gselqp->Draw("P SAME");
    drawConstRLineQP(bandLow[b], colors[b % colors.size()], 2, 3);
    drawConstRLineQP(bandHigh[b], colors[b % colors.size()], 2, 3);
    drawConstRLineQP(bandPeaks[b].q, colors[b % colors.size()], 1, 3);

    cb->cd(3);
    hQ->Draw("HIST");
    hQSmooth->SetLineColor(kRed+1); hQSmooth->SetLineWidth(2); hQSmooth->Draw("HIST SAME");
    double ymaxQ = hQ->GetMaximum()*1.08;
    for (const auto& pk : bandPeaks) {
      TLine* lpk = new TLine(pk.q, 0, pk.q, ymaxQ);
      lpk->SetLineColor(kRed+1); lpk->SetLineWidth(pk.q == bandPeaks[b].q ? 3 : 1); lpk->SetLineStyle(pk.q == bandPeaks[b].q ? 1 : 2); lpk->Draw("same");
    }
    TLine* bq1 = new TLine(bandLow[b], 0, bandLow[b], ymaxQ);
    TLine* bq2 = new TLine(bandHigh[b], 0, bandHigh[b], ymaxQ);
    bq1->SetLineColor(colors[b % colors.size()]); bq2->SetLineColor(colors[b % colors.size()]);
    bq1->SetLineWidth(3); bq2->SetLineWidth(3); bq1->SetLineStyle(2); bq2->SetLineStyle(2);
    bq1->Draw("same"); bq2->Draw("same");
    hBandQ[b]->SetLineColor(colors[b % colors.size()]); hBandQ[b]->SetLineWidth(2); hBandQ[b]->Draw("HIST SAME");

    cb->cd(4);
    gPad->SetLogz();
    hBandOrig[b]->Draw("COLZ");
    for (const auto& pk : bandPeaks) drawConstRLine(pk.q, kRed+1, 2, 1);
    drawConstRLine(bandPeaks[b].q, colors[b % colors.size()], 1, 3);
    TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.035);
    lat2.DrawLatex(0.12,0.92,Form("Angle band %zu [%s] only", b, bandIsWeak[b] ? "weak" : "strong"));
    lat2.DrawLatex(0.12,0.87,Form("N=%d events", bandCounts[b]));
    lat2.DrawLatex(0.12,0.82,Form("r peak=%.3f", bandPeaks[b].q));
    lat2.DrawLatex(0.12,0.77,Form("r range [%.3f, %.3f]", bandLow[b], bandHigh[b]));

    cb->SaveAs(outpdf);
  }

  c1->SaveAs(outpdf + "]");

  TFile fout(outroot, "RECREATE");
  hOrig->Write();
  hQP->Write();
  hUR->Write();
  hQ->Write();
  hP->Write();
  hQSmooth->Write();
  gScore->Write("g_angle_score");
  gNPeaks->Write("g_angle_npeaks");
  gValley->Write("g_angle_valley");
  for (size_t b=0; b<graphs.size(); ++b) graphs[b]->Write(Form("g_rband_%zu", b));
  for (size_t b=0; b<bandPeaks.size(); ++b) {
    hBandOrig[b]->Write();
    hBandQP[b]->Write();
    hBandQ[b]->Write();
  }
  fout.WriteObject(c1, "c_angleScanBand_assign_overview");
  fout.WriteObject(cScore, "c_angle_scan_score");
  fout.Close();

  cout << "Wrote: " << outpdf << endl;
  cout << "Wrote: " << outroot << endl;
}
