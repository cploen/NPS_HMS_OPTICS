// assign_yfp_ypfp_vbands_autoPeakSep_with_pages.C
// Diagnostic: use pooled-local PCA to estimate shared island slant, then separate islands using 1D peaks in v = across-PC1 coordinate.
//
// This is intentionally different from global PCA:
//   1) select events for one run and delta slice
//   2) standardize x=ypfp and y=yfp
//   3) find the N strongest density peaks in standardized 2D space
//   4) take only core events around each peak
//   5) compute one covariance matrix per core
//   6) average those covariance matrices equally
//   7) PC1 of the pooled covariance = shared local ellipse direction
//
// Run from repo top directory, e.g.
//   hcana -l -q 'assign_yfp_ypfp_vbands_autoPeakSep_with_pages.C(1544,-5,0,"auto_ycut")'
//
// Useful tuning, if needed:
//   hcana -l -q 'assign_yfp_ypfp_vbands_autoPeakSep_with_pages.C(1544,-5,0,"auto_ycut",7,0.55,0.45)' 
// where the last numbers are: nLocalPcaPeaks, maxVBands, peakSepZ, coreRadiusZ, minVPeakSep, minPeakFraction, autoPeakSepFrac, autoPeakSepMin, autoPeakSepMax, autoPeakSepRoughN.
// If peakSepZ < 0, the macro auto-estimates it from rough 2D density peak nearest-neighbor spacing.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TKey.h>
#include <TClass.h>
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
#include <TROOT.h>

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

struct Peak2D {
  double xz = 0.0;
  double yz = 0.0;
  double count = 0.0;
};

static std::vector<TString> SplitCSV_vbands(const TString& line) {
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

static bool ReadOpticsRunInfo_vbands(int nrun, OpticsRunInfo& info,
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

    auto tok = SplitCSV_vbands(line);
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

static TCutG* GetFirstCutG_vbands(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj && obj->InheritsFrom(TCutG::Class())) return (TCutG*)obj;
  }
  return nullptr;
}

static TCutG* LoadYtarCut_vbands(int nrun, const TString& ytarTag, int foilIndex=0) {
  TString fname = Form("cuts/ytar_delta_%d_%s_multifoil_cut.root", nrun, ytarTag.Data());
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open ytar cut file: " << fname << endl;
    return nullptr;
  }

  std::vector<TString> names = {
    Form("delta_vs_ytar_cut_foil%d", foilIndex),
    Form("ytar_delta_cut_foil%d", foilIndex),
    Form("foil%d", foilIndex),
    Form("cut_foil%d", foilIndex),
    "delta_vs_ytar_cut",
    "ytar_delta_cut"
  };

  for (auto& name : names) {
    TCutG* c = (TCutG*)f->Get(name);
    if (c) {
      cout << "Loaded ytar cut: " << fname << " :: " << name << endl;
      return c;
    }
  }

  TCutG* first = GetFirstCutG_vbands(f);
  if (first) cout << "Loaded first available TCutG from " << fname << ": " << first->GetName() << endl;
  else       cout << "WARNING: no TCutG found in " << fname << endl;
  return first;
}

static void Eigen2x2Symmetric_vbands(double sxx, double syy, double sxy,
                                       double& lam1, double& lam2,
                                       double& c, double& s,
                                       double& theta) {
  const double tr = sxx + syy;
  const double detTerm = std::sqrt((sxx - syy)*(sxx - syy) + 4.0*sxy*sxy);
  lam1 = 0.5*(tr + detTerm);
  lam2 = 0.5*(tr - detTerm);
  theta = 0.5 * std::atan2(2.0*sxy, sxx - syy);
  c = std::cos(theta);
  s = std::sin(theta);
}


static double Median_vbands(std::vector<double> vals) {
  if (vals.empty()) return 0.0;
  std::sort(vals.begin(), vals.end());
  const size_t n = vals.size();
  if (n % 2 == 1) return vals[n/2];
  return 0.5 * (vals[n/2 - 1] + vals[n/2]);
}

static double Clamp_vbands(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

// Auto-estimate a 2D non-maximum-suppression radius in standardized (xz,yz) space.
// Procedure:
//   1) take the strongest occupied histogram bins as rough candidate density centers,
//   2) for each rough candidate, compute nearest-neighbor distance,
//   3) use a fraction of the median nearest-neighbor spacing,
//   4) clamp to sane bounds so sparse/noisy slices do not go feral.
static double EstimatePeakSepZ_vbands(const std::vector<Peak2D>& sortedCandidates,
                                      int maxRoughPeaks = 25,
                                      double fracOfMedianNN = 0.75,
                                      double minSep = 0.60,
                                      double maxSep = 1.50) {
  std::vector<Peak2D> rough;
  rough.reserve(maxRoughPeaks);

  for (const auto& cand : sortedCandidates) {
    rough.push_back(cand);
    if ((int)rough.size() >= maxRoughPeaks) break;
  }

  if (rough.size() < 3) {
    cout << "WARNING: too few rough 2D candidates for auto peakSepZ; using fallback 1.15" << endl;
    return 1.15;
  }

  std::vector<double> nn;
  nn.reserve(rough.size());

  for (size_t i=0; i<rough.size(); ++i) {
    double best = std::numeric_limits<double>::max();
    for (size_t j=0; j<rough.size(); ++j) {
      if (i == j) continue;
      const double dx = rough[i].xz - rough[j].xz;
      const double dy = rough[i].yz - rough[j].yz;
      const double d = std::sqrt(dx*dx + dy*dy);
      if (d < best) best = d;
    }
    if (std::isfinite(best)) nn.push_back(best);
  }

  const double medNN = Median_vbands(nn);
  double autoSep = fracOfMedianNN * medNN;
  autoSep = Clamp_vbands(autoSep, minSep, maxSep);

  cout << "\n=== Auto peakSepZ estimate ===" << endl;
  cout << "Rough 2D candidate peaks used: " << rough.size() << endl;
  cout << "Median rough nearest-neighbor distance in z-space: " << medNN << endl;
  cout << "Fraction of median NN used: " << fracOfMedianNN << endl;
  cout << "Clamp range: [" << minSep << ", " << maxSep << "]" << endl;
  cout << "Auto-selected peakSepZ: " << autoSep << endl;

  return autoSep;
}

void assign_yfp_ypfp_vbands_autoPeakSep_with_pages(Int_t nrun=1544,
                                               Double_t deltaMin=-5.0,
                                               Double_t deltaMax=0.0,
                                               TString ytarTag="auto_ycut",
                                               Int_t nLocalPcaPeaks=7,
                                               Int_t maxVBands=9,
                                               Double_t peakSepZ=1.15,
                                               Double_t coreRadiusZ=0.30,
                                               Double_t minVPeakSep=0.18,
                                               Double_t minPeakFraction=0.05,
                                               Double_t autoPeakSepFrac=0.75,
                                               Double_t autoPeakSepMin=0.60,
                                               Double_t autoPeakSepMax=1.50,
                                               Int_t autoPeakSepRoughN=25,
                                               Bool_t useYtarCut=true,
                                               Int_t foilIndex=0,
                                               Long64_t maxEvents=-1) {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  OpticsRunInfo info;
  ReadOpticsRunInfo_vbands(nrun, info);

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
  if (useYtarCut) ytarCut = LoadYtarCut_vbands(nrun, ytarTag, foilIndex);

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
      // ytar-delta cuts are x=ytar, y=delta.
      if (!ytarCut->IsInside(ytar, delta)) continue;
    }
    nPassYtar++;

    if (!(delta >= deltaMin && delta < deltaMax)) continue;
    nPassDelta++;

    xvals.push_back(ypfp);
    yvals.push_back(yfp);
  }

  const Long64_t N = xvals.size();
  cout << "\n=== Pooled local PC1 + v-band diagnostic ===" << endl;
  cout << "Run: " << nrun << endl;
  cout << "Input: " << inroot << endl;
  cout << "Metadata opticsID: " << info.opticsID << ", NumFoil: " << info.numFoil << endl;
  cout << "Delta slice: [" << deltaMin << ", " << deltaMax << ") %" << endl;
  cout << "Requested local PCA density peaks: " << nLocalPcaPeaks << endl;
  cout << "Peak suppression radius in z-space request: " << peakSepZ << endl;
  cout << "Core radius in z-space: " << coreRadiusZ << endl;
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

  double xzmin=*std::min_element(xzvals.begin(), xzvals.end());
  double xzmax=*std::max_element(xzvals.begin(), xzvals.end());
  double yzmin=*std::min_element(yzvals.begin(), yzvals.end());
  double yzmax=*std::max_element(yzvals.begin(), yzvals.end());
  double xzpad=0.05*(xzmax-xzmin), yzpad=0.05*(yzmax-yzmin);
  xzmin-=xzpad; xzmax+=xzpad; yzmin-=yzpad; yzmax+=yzpad;

  const int nBinsPeak = 260;
  TH2D* hZ = new TH2D("hZ", Form("Standardized density, Run %d;x_{z}=(ypfp-mean)/#sigma; y_{z}=(yfp-mean)/#sigma", nrun),
                      nBinsPeak, xzmin, xzmax, nBinsPeak, yzmin, yzmax);
  for (Long64_t i=0; i<N; ++i) hZ->Fill(xzvals[i], yzvals[i]);

  // Density peak detection: sort bins by occupancy, then non-maximum suppress nearby bins.
  std::vector<Peak2D> candidates;
  candidates.reserve(nBinsPeak*nBinsPeak);
  for (int ix=1; ix<=hZ->GetNbinsX(); ++ix) {
    for (int iy=1; iy<=hZ->GetNbinsY(); ++iy) {
      const double cnt = hZ->GetBinContent(ix,iy);
      if (cnt <= 0) continue;
      candidates.push_back({hZ->GetXaxis()->GetBinCenter(ix), hZ->GetYaxis()->GetBinCenter(iy), cnt});
    }
  }
  std::sort(candidates.begin(), candidates.end(), [](const Peak2D& a, const Peak2D& b){ return a.count > b.count; });

  if (peakSepZ < 0.0) {
    peakSepZ = EstimatePeakSepZ_vbands(candidates,
                                       autoPeakSepRoughN,
                                       autoPeakSepFrac,
                                       autoPeakSepMin,
                                       autoPeakSepMax);
  } else {
    cout << "\n=== Fixed peakSepZ ===" << endl;
    cout << "Using user-supplied peakSepZ: " << peakSepZ << endl;
  }

  std::vector<Peak2D> peaks;
  for (const auto& cand : candidates) {
    bool tooClose = false;
    for (const auto& pk : peaks) {
      const double dx = cand.xz - pk.xz;
      const double dy = cand.yz - pk.yz;
      if (std::sqrt(dx*dx + dy*dy) < peakSepZ) { tooClose = true; break; }
    }
    if (tooClose) continue;
    peaks.push_back(cand);
    if ((int)peaks.size() >= nLocalPcaPeaks) break;
  }

  if ((int)peaks.size() < nLocalPcaPeaks) {
    cout << "WARNING: found only " << peaks.size() << " peaks." << endl;
  }

  cout << "\nSelected density peaks in z-space using peakSepZ=" << peakSepZ << ":" << endl;
  for (size_t k=0; k<peaks.size(); ++k) {
    const double xraw = mx + sx*peaks[k].xz;
    const double yraw = my + sy*peaks[k].yz;
    cout << "  peak " << k << ": z=(" << peaks[k].xz << ", " << peaks[k].yz
         << "), raw=(ypfp " << xraw << ", yfp " << yraw << "), bin count=" << peaks[k].count << endl;
  }

  // Assign core events to nearest peak, but only inside coreRadiusZ.
  std::vector< std::vector<int> > coreIdx(peaks.size());
  for (Long64_t i=0; i<N; ++i) {
    int best = -1;
    double bestR2 = std::numeric_limits<double>::max();
    for (size_t k=0; k<peaks.size(); ++k) {
      const double dx = xzvals[i] - peaks[k].xz;
      const double dy = yzvals[i] - peaks[k].yz;
      const double r2 = dx*dx + dy*dy;
      if (r2 < bestR2) { bestR2 = r2; best = (int)k; }
    }
    if (best >= 0 && bestR2 <= coreRadiusZ*coreRadiusZ) coreIdx[best].push_back((int)i);
  }

  // Pooled local covariance in standardized coordinates, equal weight per valid island core.
  double pool_xx=0, pool_yy=0, pool_xy=0;
  int nValidCores = 0;
  const int minCoreEvents = 50;
  cout << "\nCore covariances in z-space:" << endl;
  for (size_t k=0; k<coreIdx.size(); ++k) {
    const int nk = (int)coreIdx[k].size();
    if (nk < minCoreEvents) {
      cout << "  core " << k << ": N=" << nk << " SKIPPED, too few events" << endl;
      continue;
    }
    double cx=0, cy=0;
    for (int idx : coreIdx[k]) { cx += xzvals[idx]; cy += yzvals[idx]; }
    cx /= (double)nk;
    cy /= (double)nk;

    double cxx=0, cyy=0, cxy=0;
    for (int idx : coreIdx[k]) {
      const double dx = xzvals[idx] - cx;
      const double dy = yzvals[idx] - cy;
      cxx += dx*dx;
      cyy += dy*dy;
      cxy += dx*dy;
    }
    cxx /= (double)(nk-1);
    cyy /= (double)(nk-1);
    cxy /= (double)(nk-1);

    double lamA, lamB, cc, ss, th;
    Eigen2x2Symmetric_vbands(cxx, cyy, cxy, lamA, lamB, cc, ss, th);
    double thDeg = th*180.0/TMath::Pi();
    cout << "  core " << k << ": N=" << nk
         << " cov=[" << cxx << ", " << cxy << "; " << cxy << ", " << cyy << "]"
         << " local angle=" << thDeg << " deg" << endl;

    pool_xx += cxx;
    pool_yy += cyy;
    pool_xy += cxy;
    nValidCores++;
  }

  if (nValidCores < 2) {
    cout << "ERROR: too few valid cores for pooled local PCA. Try increasing coreRadiusZ." << endl;
    return;
  }

  pool_xx /= (double)nValidCores;
  pool_yy /= (double)nValidCores;
  pool_xy /= (double)nValidCores;

  double lam1, lam2, c, s, theta;
  Eigen2x2Symmetric_vbands(pool_xx, pool_yy, pool_xy, lam1, lam2, c, s, theta);

  // Choose sign for prettier convention: positive slope in z-space if possible.
  if (c*s < 0) { c = -c; s = -s; }

  const double thetaDegZ = std::atan2(s,c) * 180.0 / TMath::Pi();

  double dxRaw = sx * c;
  double dyRaw = sy * s;
  const double normRaw = std::sqrt(dxRaw*dxRaw + dyRaw*dyRaw);
  dxRaw /= normRaw;
  dyRaw /= normRaw;
  const double thetaDegRaw = std::atan2(dyRaw, dxRaw) * 180.0 / TMath::Pi();

  // PC2 raw-space direction corresponding to z-space (-s,c).
  double dxRaw2 = -sx * s;
  double dyRaw2 =  sy * c;
  const double normRaw2 = std::sqrt(dxRaw2*dxRaw2 + dyRaw2*dyRaw2);
  dxRaw2 /= normRaw2;
  dyRaw2 /= normRaw2;

  cout << "\n=== Pooled local PC1 result ===" << endl;
  cout << "Mean x=<ypfp>: " << mx << endl;
  cout << "Mean y=<yfp>:  " << my << endl;
  cout << "Raw RMS: sigma_ypfp=" << sx << ", sigma_yfp=" << sy << endl;
  cout << "Valid cores used: " << nValidCores << endl;
  cout << "Pooled covariance in z-space:" << endl;
  cout << "  [ " << pool_xx << "   " << pool_xy << " ]" << endl;
  cout << "  [ " << pool_xy << "   " << pool_yy << " ]" << endl;
  cout << "lambda1/lambda2 pooled: " << lam1 << " / " << lam2 << endl;
  cout << "Pooled local PC1 angle in z-space: " << thetaDegZ << " deg" << endl;
  cout << "Pooled local PC1 vector in z-space: (" << c << ", " << s << ")" << endl;
  cout << "Mapped raw-space line direction: (" << dxRaw << ", " << dyRaw << ") in (ypfp,yfp)" << endl;
  cout << "Mapped raw-space line angle from +ypfp axis: " << thetaDegRaw << " deg" << endl;

  double xmin=*std::min_element(xvals.begin(), xvals.end());
  double xmax=*std::max_element(xvals.begin(), xvals.end());
  double ymin=*std::min_element(yvals.begin(), yvals.end());
  double ymax=*std::max_element(yvals.begin(), yvals.end());
  double xpad=0.10*(xmax-xmin), ypad=0.10*(ymax-ymin);
  xmin-=xpad; xmax+=xpad; ymin-=ypad; ymax+=ypad;

  TH2D* hOrig = new TH2D("hOrig", Form("Run %d, %.1f < #delta < %.1f;ypfp;yfp", nrun, deltaMin, deltaMax),
                         220, xmin, xmax, 220, ymin, ymax);

  std::vector<double> uvals, vvals;
  uvals.reserve(N); vvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    const double u =  xzvals[i]*c + yzvals[i]*s;
    const double v = -xzvals[i]*s + yzvals[i]*c;
    uvals.push_back(u);
    vvals.push_back(v);
    hOrig->Fill(xvals[i], yvals[i]);
  }

  double umin=*std::min_element(uvals.begin(), uvals.end());
  double umax=*std::max_element(uvals.begin(), uvals.end());
  double vmin=*std::min_element(vvals.begin(), vvals.end());
  double vmax=*std::max_element(vvals.begin(), vvals.end());
  double upad=0.10*(umax-umin), vpad=0.10*(vmax-vmin);
  umin-=upad; umax+=upad; vmin-=vpad; vmax+=vpad;

  TH2D* hRot = new TH2D("hRot", Form("Pooled-local-PCA rotated coords, Run %d;u = pooled local PC1 coordinate;v = pooled local PC2 coordinate", nrun),
                        220, umin, umax, 220, vmin, vmax);
  TH1D* hU = new TH1D("hU", Form("Run %d;u = along pooled local PC1;counts", nrun), 220, umin, umax);
  TH1D* hV = new TH1D("hV", Form("Run %d;v = across pooled local PC1;counts", nrun), 220, vmin, vmax);

  for (Long64_t i=0; i<N; ++i) {
    hRot->Fill(uvals[i], vvals[i]);
    hU->Fill(uvals[i]);
    hV->Fill(vvals[i]);
  }

  // --- Find 1D peaks in v, then assign v-bands by midpoint boundaries. ---
  TH1D* hVSmooth = (TH1D*)hV->Clone("hVSmooth");
  hVSmooth->Smooth(2);

  struct VPeak { double v=0.0; double height=0.0; int bin=0; };
  std::vector<VPeak> vcands;
  const double vmaxContent = hVSmooth->GetMaximum();
  const double minHeight = minPeakFraction * vmaxContent;
  for (int ib=2; ib<hVSmooth->GetNbinsX(); ++ib) {
    const double ym = hVSmooth->GetBinContent(ib-1);
    const double y0 = hVSmooth->GetBinContent(ib);
    const double yp = hVSmooth->GetBinContent(ib+1);
    if (y0 > ym && y0 >= yp && y0 > minHeight) {
      vcands.push_back({hVSmooth->GetBinCenter(ib), y0, ib});
    }
  }
  std::sort(vcands.begin(), vcands.end(), [](const VPeak& a, const VPeak& b){ return a.height > b.height; });

  std::vector<VPeak> vpeaks;
  int nRejectedTooClose = 0;
  for (const auto& cand : vcands) {
    bool tooClose=false;
    for (const auto& pk : vpeaks) {
      if (std::fabs(cand.v - pk.v) < minVPeakSep) { tooClose=true; break; }
    }
    if (tooClose) { nRejectedTooClose++; continue; }
    vpeaks.push_back(cand);
    if ((int)vpeaks.size() >= maxVBands) break;
  }
  std::sort(vpeaks.begin(), vpeaks.end(), [](const VPeak& a, const VPeak& b){ return a.v < b.v; });

  std::vector<double> vbounds;
  for (size_t i=0; i+1<vpeaks.size(); ++i) vbounds.push_back(0.5*(vpeaks[i].v + vpeaks[i+1].v));

  std::vector<int> bandIndex(N, -1);
  std::vector<int> bandCounts(vpeaks.size(), 0);
  for (Long64_t i=0; i<N; ++i) {
    int b=0;
    while (b < (int)vbounds.size() && vvals[i] > vbounds[b]) b++;
    if (b >= 0 && b < (int)vpeaks.size()) { bandIndex[i]=b; bandCounts[b]++; }
  }

  cout << "\n=== v-band peak assignment ===" << endl;
  cout << "Max allowed v peaks: " << maxVBands << endl;
  cout << "Minimum v-peak separation: " << minVPeakSep << endl;
  cout << "Minimum v-peak height fraction: " << minPeakFraction << endl;
  cout << "Candidate v peaks above height threshold: " << vcands.size() << endl;
  cout << "Rejected candidate v peaks as too close to stronger accepted peak: " << nRejectedTooClose << endl;
  cout << "Accepted v peaks: " << vpeaks.size() << endl;
  for (size_t i=0; i<vpeaks.size(); ++i) {
    cout << "  vband " << i << ": peak v=" << vpeaks[i].v
         << ", smoothed height=" << vpeaks[i].height
         << ", assigned events=" << bandCounts[i] << endl;
  }
  cout << "Boundaries:";
  for (double b : vbounds) cout << " " << b;
  cout << endl;

  TString outbase = Form("plots/yfp_ypfp_vband_autoPeakSep_pages_run%d_delta_%g_to_%g", nrun, deltaMin, deltaMax);
  outbase.ReplaceAll("-", "m");
  outbase.ReplaceAll(".", "p");
  TString outpdf = outbase + ".pdf";
  TString outroot = outbase + ".root";
  gSystem->mkdir("plots", kTRUE);

  TCanvas* c1 = new TCanvas("c_vband_assign", "Pooled local PC1 + v-band assignment", 1200, 950);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetLogz();
  hOrig->Draw("COLZ");

  double L = 0.55 * std::sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
  TLine* pc1 = new TLine(mx - L*dxRaw, my - L*dyRaw, mx + L*dxRaw, my + L*dyRaw);
  pc1->SetLineColor(kOrange+7);
  pc1->SetLineWidth(4);
  pc1->Draw("same");

  TLine* pc2 = new TLine(mx - L*dxRaw2, my - L*dyRaw2, mx + L*dxRaw2, my + L*dyRaw2);
  pc2->SetLineColor(kMagenta+2);
  pc2->SetLineWidth(2);
  pc2->SetLineStyle(2);
  pc2->Draw("same");

  TMarker* cen = new TMarker(mx, my, 29);
  cen->SetMarkerColor(kBlack);
  cen->SetMarkerSize(1.5);
  cen->Draw("same");

  for (size_t k=0; k<peaks.size(); ++k) {
    const double xraw = mx + sx*peaks[k].xz;
    const double yraw = my + sy*peaks[k].yz;
    TMarker* m = new TMarker(xraw, yraw, 20);
    m->SetMarkerColor(kBlack);
    m->SetMarkerSize(0.9);
    m->Draw("same");
  }

  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.032);
  lat.DrawLatex(0.12,0.92,Form("Pooled PC1 %.2f deg z-space; v peaks=%zu", thetaDegZ, vpeaks.size()));

  c1->cd(2);
  gPad->SetLogz();
  hRot->Draw("COLZ");
  for (const auto& pk : vpeaks) {
    TLine* lpk = new TLine(umin, pk.v, umax, pk.v);
    lpk->SetLineColor(kRed+1);
    lpk->SetLineWidth(2);
    lpk->Draw("same");
  }
  for (double b : vbounds) {
    TLine* lb = new TLine(umin, b, umax, b);
    lb->SetLineColor(kRed+1);
    lb->SetLineStyle(2);
    lb->SetLineWidth(2);
    lb->Draw("same");
  }

  c1->cd(3);
  hV->Draw("HIST");
  hVSmooth->SetLineColor(kRed+1);
  hVSmooth->SetLineWidth(2);
  hVSmooth->Draw("HIST SAME");
  for (const auto& pk : vpeaks) {
    TLine* lpk = new TLine(pk.v, 0, pk.v, hV->GetMaximum()*1.05);
    lpk->SetLineColor(kRed+1);
    lpk->SetLineWidth(2);
    lpk->Draw("same");
  }
  for (double b : vbounds) {
    TLine* lb = new TLine(b, 0, b, hV->GetMaximum()*1.05);
    lb->SetLineColor(kRed+1);
    lb->SetLineStyle(2);
    lb->Draw("same");
  }

  c1->cd(4);
  TH2D* hFrame = new TH2D("hFrameBands", Form("v-band assignments, Run %d;ypfp;yfp", nrun), 10, xmin, xmax, 10, ymin, ymax);
  hFrame->SetMinimum(0); hFrame->SetMaximum(1);
  hFrame->Draw("AXIS");
  std::vector<int> colors = {kRed+1, kOrange+7, kSpring+5, kGreen+2, kAzure+7, kBlue+1, kMagenta+1, kViolet+7, kCyan+2, kGray+2};
  std::vector<TGraph*> graphs(vpeaks.size(), nullptr);
  for (size_t b=0; b<vpeaks.size(); ++b) graphs[b] = new TGraph();
  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)graphs.size()) continue;
    int p = graphs[b]->GetN();
    graphs[b]->SetPoint(p, xvals[i], yvals[i]);
  }
  for (size_t b=0; b<graphs.size(); ++b) {
    graphs[b]->SetMarkerStyle(20);
    graphs[b]->SetMarkerSize(0.18);
    graphs[b]->SetMarkerColor(colors[b % colors.size()]);
    graphs[b]->Draw("P SAME");
  }

  // Build per-band diagnostics and write multi-page PDF.
  std::vector<double> bandLow(vpeaks.size(), vmin), bandHigh(vpeaks.size(), vmax);
  for (size_t b=0; b<vpeaks.size(); ++b) {
    bandLow[b]  = (b==0) ? vmin : vbounds[b-1];
    bandHigh[b] = (b+1==vpeaks.size()) ? vmax : vbounds[b];
  }

  c1->SaveAs(outpdf + "[");
  c1->SaveAs(outpdf);

  std::vector<TH2D*> hBandOrig(vpeaks.size(), nullptr);
  std::vector<TH2D*> hBandRot(vpeaks.size(), nullptr);
  std::vector<TH1D*> hBandU(vpeaks.size(), nullptr);
  std::vector<TH1D*> hBandV(vpeaks.size(), nullptr);

  for (size_t b=0; b<vpeaks.size(); ++b) {
    hBandOrig[b] = new TH2D(Form("hBandOrig_%zu", b),
      Form("Band %zu only in original coords, Run %d;ypfp;yfp", b, nrun),
      220, xmin, xmax, 220, ymin, ymax);
    hBandRot[b] = new TH2D(Form("hBandRot_%zu", b),
      Form("Band %zu only in rotated coords, Run %d;u;v", b, nrun),
      220, umin, umax, 220, vmin, vmax);
    hBandU[b] = new TH1D(Form("hBandU_%zu", b),
      Form("Band %zu u distribution, Run %d;u;counts", b, nrun),
      220, umin, umax);
    hBandV[b] = new TH1D(Form("hBandV_%zu", b),
      Form("Band %zu v distribution, Run %d;v;counts", b, nrun),
      220, vmin, vmax);
  }

  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)vpeaks.size()) continue;
    hBandOrig[b]->Fill(xvals[i], yvals[i]);
    hBandRot[b]->Fill(uvals[i], vvals[i]);
    hBandU[b]->Fill(uvals[i]);
    hBandV[b]->Fill(vvals[i]);
  }

  for (size_t b=0; b<vpeaks.size(); ++b) {
    TCanvas* cb = new TCanvas(Form("c_band_%zu", b), Form("Band %zu diagnostics", b), 1200, 950);
    cb->Divide(2,2);

    cb->cd(1);
    gPad->SetLogz();
    hOrig->Draw("COLZ");
    TGraph* gall = new TGraph();
    TGraph* gsel = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int p = gall->GetN();
      gall->SetPoint(p, xvals[i], yvals[i]);
      if (bandIndex[i] == (int)b) {
        int q = gsel->GetN();
        gsel->SetPoint(q, xvals[i], yvals[i]);
      }
    }
    gall->SetMarkerStyle(20); gall->SetMarkerSize(0.10); gall->SetMarkerColor(kGray+1);
    gall->Draw("P SAME");
    gsel->SetMarkerStyle(20); gsel->SetMarkerSize(0.22); gsel->SetMarkerColor(colors[b % colors.size()]);
    gsel->Draw("P SAME");
    pc1->Draw("same");
    pc2->Draw("same");
    lat.DrawLatex(0.12,0.92,Form("Band %zu highlighted in original coords", b));
    lat.DrawLatex(0.12,0.88,Form("peak v=%.3f, range [%.3f, %.3f], N=%d", vpeaks[b].v, bandLow[b], bandHigh[b], bandCounts[b]));

    cb->cd(2);
    gPad->SetLogz();
    hRot->Draw("COLZ");
    TGraph* galluv = new TGraph();
    TGraph* gseluv = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int p = galluv->GetN();
      galluv->SetPoint(p, uvals[i], vvals[i]);
      if (bandIndex[i] == (int)b) {
        int q = gseluv->GetN();
        gseluv->SetPoint(q, uvals[i], vvals[i]);
      }
    }
    galluv->SetMarkerStyle(20); galluv->SetMarkerSize(0.10); galluv->SetMarkerColor(kGray+1);
    galluv->Draw("P SAME");
    gseluv->SetMarkerStyle(20); gseluv->SetMarkerSize(0.22); gseluv->SetMarkerColor(colors[b % colors.size()]);
    gseluv->Draw("P SAME");
    for (const auto& pk : vpeaks) {
      TLine* lpk = new TLine(umin, pk.v, umax, pk.v);
      lpk->SetLineColor(kRed+1);
      lpk->SetLineWidth(pk.v == vpeaks[b].v ? 3 : 1);
      lpk->SetLineStyle(pk.v == vpeaks[b].v ? 1 : 2);
      lpk->Draw("same");
    }
    TLine* llo = new TLine(umin, bandLow[b], umax, bandLow[b]);
    llo->SetLineColor(colors[b % colors.size()]); llo->SetLineWidth(3); llo->SetLineStyle(2); llo->Draw("same");
    TLine* lhi = new TLine(umin, bandHigh[b], umax, bandHigh[b]);
    lhi->SetLineColor(colors[b % colors.size()]); lhi->SetLineWidth(3); lhi->SetLineStyle(2); lhi->Draw("same");

    cb->cd(3);
    hV->Draw("HIST");
    hVSmooth->SetLineColor(kRed+1);
    hVSmooth->SetLineWidth(2);
    hVSmooth->Draw("HIST SAME");
    double ymaxV = hV->GetMaximum()*1.08;
    for (const auto& pk : vpeaks) {
      TLine* lpk = new TLine(pk.v, 0, pk.v, ymaxV);
      lpk->SetLineColor(kRed+1);
      lpk->SetLineWidth(pk.v == vpeaks[b].v ? 3 : 1);
      lpk->SetLineStyle(pk.v == vpeaks[b].v ? 1 : 2);
      lpk->Draw("same");
    }
    TLine* bv1 = new TLine(bandLow[b], 0, bandLow[b], ymaxV);
    TLine* bv2 = new TLine(bandHigh[b], 0, bandHigh[b], ymaxV);
    bv1->SetLineColor(colors[b % colors.size()]); bv2->SetLineColor(colors[b % colors.size()]);
    bv1->SetLineWidth(3); bv2->SetLineWidth(3);
    bv1->SetLineStyle(2); bv2->SetLineStyle(2);
    bv1->Draw("same"); bv2->Draw("same");
    hBandV[b]->SetLineColor(colors[b % colors.size()]);
    hBandV[b]->SetLineWidth(2);
    hBandV[b]->Draw("HIST SAME");

    cb->cd(4);
    gPad->SetLogz();
    hBandOrig[b]->Draw("COLZ");
    pc1->Draw("same");
    pc2->Draw("same");
    TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.035);
    lat2.DrawLatex(0.12,0.92,Form("Band %zu only", b));
    lat2.DrawLatex(0.12,0.87,Form("N=%d events", bandCounts[b]));
    lat2.DrawLatex(0.12,0.82,Form("v peak=%.3f", vpeaks[b].v));
    lat2.DrawLatex(0.12,0.77,Form("v range [%.3f, %.3f]", bandLow[b], bandHigh[b]));

    cb->SaveAs(outpdf);
  }

  c1->SaveAs(outpdf + "]");

  TFile fout(outroot, "RECREATE");
  hOrig->Write();
  hZ->Write();
  hRot->Write();
  hU->Write();
  hV->Write();
  hVSmooth->Write();
  for (size_t b=0; b<graphs.size(); ++b) graphs[b]->Write(Form("g_vband_%zu", b));
  for (size_t b=0; b<vpeaks.size(); ++b) {
    hBandOrig[b]->Write();
    hBandRot[b]->Write();
    hBandU[b]->Write();
    hBandV[b]->Write();
  }
  fout.WriteObject(c1, "c_vband_assign_overview");
  fout.Close();

  cout << "Wrote: " << outpdf << endl;
  cout << "Wrote: " << outroot << endl;
}
