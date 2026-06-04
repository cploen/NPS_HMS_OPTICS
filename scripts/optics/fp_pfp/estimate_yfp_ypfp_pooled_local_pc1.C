// estimate_yfp_ypfp_pooled_local_pc1.C
// Diagnostic: estimate the shared *local island* long-axis direction in HMS ypfp-yfp space.
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
//   hcana -l -q 'estimate_yfp_ypfp_pooled_local_pc1.C(1544,-5,0,"auto_ycut")'
//
// Useful tuning, if needed:
//   hcana -l -q 'estimate_yfp_ypfp_pooled_local_pc1.C(1544,-5,0,"auto_ycut",7,0.55,0.45)'
// where the last numbers are: nPeaksTarget, peakSepZ, coreRadiusZ.

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

static std::vector<TString> SplitCSV_localpc1(const TString& line) {
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

static bool ReadOpticsRunInfo_localpc1(int nrun, OpticsRunInfo& info,
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

    auto tok = SplitCSV_localpc1(line);
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

static TCutG* GetFirstCutG_localpc1(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj && obj->InheritsFrom(TCutG::Class())) return (TCutG*)obj;
  }
  return nullptr;
}

static TCutG* LoadYtarCut_localpc1(int nrun, const TString& ytarTag, int foilIndex=0) {
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

  TCutG* first = GetFirstCutG_localpc1(f);
  if (first) cout << "Loaded first available TCutG from " << fname << ": " << first->GetName() << endl;
  else       cout << "WARNING: no TCutG found in " << fname << endl;
  return first;
}

static void Eigen2x2Symmetric_localpc1(double sxx, double syy, double sxy,
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

void estimate_yfp_ypfp_pooled_local_pc1(Int_t nrun=1544,
                                        Double_t deltaMin=-5.0,
                                        Double_t deltaMax=0.0,
                                        TString ytarTag="auto_ycut",
                                        Int_t nPeaksTarget=7,
                                        Double_t peakSepZ=0.55,
                                        Double_t coreRadiusZ=0.45,
                                        Bool_t useYtarCut=true,
                                        Int_t foilIndex=0,
                                        Long64_t maxEvents=-1) {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  OpticsRunInfo info;
  ReadOpticsRunInfo_localpc1(nrun, info);

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
  if (useYtarCut) ytarCut = LoadYtarCut_localpc1(nrun, ytarTag, foilIndex);

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
  cout << "\n=== Pooled local PC1 diagnostic ===" << endl;
  cout << "Run: " << nrun << endl;
  cout << "Input: " << inroot << endl;
  cout << "Metadata opticsID: " << info.opticsID << ", NumFoil: " << info.numFoil << endl;
  cout << "Delta slice: [" << deltaMin << ", " << deltaMax << ") %" << endl;
  cout << "Requested density peaks: " << nPeaksTarget << endl;
  cout << "Peak suppression radius in z-space: " << peakSepZ << endl;
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
    if ((int)peaks.size() >= nPeaksTarget) break;
  }

  if ((int)peaks.size() < nPeaksTarget) {
    cout << "WARNING: found only " << peaks.size() << " peaks." << endl;
  }

  cout << "\nSelected density peaks in z-space:" << endl;
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
    Eigen2x2Symmetric_localpc1(cxx, cyy, cxy, lamA, lamB, cc, ss, th);
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
  Eigen2x2Symmetric_localpc1(pool_xx, pool_yy, pool_xy, lam1, lam2, c, s, theta);

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

  TString outbase = Form("plots/yfp_ypfp_pooled_local_pc1_run%d_delta_%g_to_%g", nrun, deltaMin, deltaMax);
  outbase.ReplaceAll("-", "m");
  outbase.ReplaceAll(".", "p");
  TString outpdf = outbase + ".pdf";
  TString outroot = outbase + ".root";
  gSystem->mkdir("plots", kTRUE);

  TCanvas* c1 = new TCanvas("c_pooled_local_pc1", "Pooled local PC1 diagnostic", 1200, 950);
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
  lat.DrawLatex(0.12,0.92,Form("Pooled local PC1: %.2f deg z-space; raw overlay %.2f deg; cores=%d", thetaDegZ, thetaDegRaw, nValidCores));

  c1->cd(2);
  gPad->SetLogz();
  hRot->Draw("COLZ");

  c1->cd(3);
  hU->Draw("HIST");

  c1->cd(4);
  hV->Draw("HIST");

  c1->SaveAs(outpdf);

  TFile fout(outroot, "RECREATE");
  hOrig->Write();
  hZ->Write();
  hRot->Write();
  hU->Write();
  hV->Write();
  fout.WriteObject(c1, "c_pooled_local_pc1");
  fout.Close();

  cout << "Wrote: " << outpdf << endl;
  cout << "Wrote: " << outroot << endl;
}
