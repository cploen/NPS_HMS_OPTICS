// estimate_yfp_ypfp_scaled_pc1.C
// Diagnostic: estimate scaled/global PC1 shared long-axis direction in the HMS ypfp-yfp event cloud.
// Default test: run 1544, delta in [-5,0] %, single foil at ztar ~= 0.
//
// Run from repo top directory, e.g.
//   hcana -l -q 'estimate_yfp_ypfp_scaled_pc1.C(1544,-5,0,"auto_ycut")'
// or
//   root -l -q 'estimate_yfp_ypfp_scaled_pc1.C(1544,-5,0,"auto_ycut")'

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
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

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

static std::vector<TString> SplitCSV(const TString& line) {
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

static bool ReadOpticsRunInfo(int nrun, OpticsRunInfo& info,
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

    auto tok = SplitCSV(line);
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

static TCutG* GetFirstCutG(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj && obj->InheritsFrom(TCutG::Class())) return (TCutG*)obj;
  }
  return nullptr;
}

static TCutG* LoadYtarCut(int nrun, const TString& ytarTag, int foilIndex=0) {
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

  TCutG* first = GetFirstCutG(f);
  if (first) cout << "Loaded first available TCutG from " << fname << ": " << first->GetName() << endl;
  else       cout << "WARNING: no TCutG found in " << fname << endl;
  return first;
}

void estimate_yfp_ypfp_scaled_pc1(Int_t nrun=1544,
                           Double_t deltaMin=-5.0,
                           Double_t deltaMax=0.0,
                           TString ytarTag="auto_ycut",
                           Bool_t useYtarCut=true,
                           Int_t foilIndex=0,
                           Long64_t maxEvents=-1) {

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  OpticsRunInfo info;
  ReadOpticsRunInfo(nrun, info);

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
  if (useYtarCut) ytarCut = LoadYtarCut(nrun, ytarTag, foilIndex);

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

  // Match the visual ypfp-yfp plot convention used in the existing optics macros:
  // x = ypfp, y = yfp.
  std::vector<double> xvals; // ypfp
  std::vector<double> yvals; // yfp

  Long64_t nentries = T->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

  Long64_t nPassBasic=0, nPassYtar=0, nPassDelta=0;
  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    // Baseline electron / good-event cuts used in this optics workflow.
    if (!(sumnpe > 6.0 && etracknorm > 0.65)) continue;
    nPassBasic++;

    if (useYtarCut && ytarCut) {
      // Existing ytar-delta cuts are defined as x=ytar, y=delta.
      if (!ytarCut->IsInside(ytar, delta)) continue;
    }
    nPassYtar++;

    if (!(delta >= deltaMin && delta < deltaMax)) continue;
    nPassDelta++;

    xvals.push_back(ypfp);
    yvals.push_back(yfp);
  }

  const Long64_t N = xvals.size();
  cout << "\n=== Scaled PC1 diagnostic ===" << endl;
  cout << "Run: " << nrun << endl;
  cout << "Input: " << inroot << endl;
  cout << "Metadata opticsID: " << info.opticsID << ", NumFoil: " << info.numFoil << endl;
  cout << "Delta slice: [" << deltaMin << ", " << deltaMax << ") %" << endl;
  cout << "Events after PID/basic: " << nPassBasic << endl;
  cout << "Events after ytar cut:  " << nPassYtar << endl;
  cout << "Events in delta slice: " << nPassDelta << endl;
  cout << "Events used for PCA:   " << N << endl;

  if (N < 10) {
    cout << "ERROR: too few events for PCA." << endl;
    return;
  }

  double mx=0, my=0;
  for (Long64_t i=0; i<N; ++i) { mx += xvals[i]; my += yvals[i]; }
  mx /= (double)N;
  my /= (double)N;

  // Raw RMS values. These are used to standardize the coordinates before PCA:
  //   xz = (ypfp - <ypfp>) / sigma_ypfp
  //   yz = (yfp  - <yfp>)  / sigma_yfp
  // This prevents the numerically large yfp coordinate from dominating the PCA.
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
    cout << "ERROR: zero RMS in one coordinate; cannot perform scaled PCA." << endl;
    return;
  }

  double sxx=0, syy=0, sxy=0;
  for (Long64_t i=0; i<N; ++i) {
    const double xz = (xvals[i] - mx) / sx;
    const double yz = (yvals[i] - my) / sy;
    sxx += xz*xz;
    syy += yz*yz;
    sxy += xz*yz;
  }
  sxx /= (double)(N-1);
  syy /= (double)(N-1);
  sxy /= (double)(N-1);

  // 2D symmetric eigensystem by closed form, now in standardized coordinates.
  const double tr = sxx + syy;
  const double detTerm = std::sqrt((sxx - syy)*(sxx - syy) + 4.0*sxy*sxy);
  const double lam1 = 0.5*(tr + detTerm);
  const double lam2 = 0.5*(tr - detTerm);

  // PC1 angle relative to standardized x=ypfp axis.
  const double theta = 0.5 * std::atan2(2.0*sxy, sxx - syy);
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const double thetaDegScaled = theta * 180.0 / TMath::Pi();

  // Direction of the scaled-PC1 line when mapped back into raw plot coordinates.
  // If z-space direction is (c,s), raw-space direction is (sx*c, sy*s).
  double dxRaw = sx * c;
  double dyRaw = sy * s;
  const double normRaw = std::sqrt(dxRaw*dxRaw + dyRaw*dyRaw);
  dxRaw /= normRaw;
  dyRaw /= normRaw;
  const double thetaDegRaw = std::atan2(dyRaw, dxRaw) * 180.0 / TMath::Pi();

  cout << "Mean x=<ypfp>: " << mx << endl;
  cout << "Mean y=<yfp>:  " << my << endl;
  cout << "Raw covariance matrix in plot coords [x=ypfp,y=yfp]:" << endl;
  cout << "  [ " << raw_sxx << "   " << raw_sxy << " ]" << endl;
  cout << "  [ " << raw_sxy << "   " << raw_syy << " ]" << endl;
  cout << "Raw RMS: sigma_ypfp=" << sx << ", sigma_yfp=" << sy << endl;
  cout << "Scaled covariance/correlation matrix [x=(ypfp-mean)/sigma, y=(yfp-mean)/sigma]:" << endl;
  cout << "  [ " << sxx << "   " << sxy << " ]" << endl;
  cout << "  [ " << sxy << "   " << syy << " ]" << endl;
  cout << "lambda1/lambda2 scaled: " << lam1 << " / " << lam2 << endl;
  cout << "Scaled PC1 angle from +standardized ypfp axis: " << thetaDegScaled << " deg" << endl;
  cout << "Scaled PC1 vector in z-space: (" << c << ", " << s << ")" << endl;
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
    const double xz = (xvals[i] - mx) / sx;
    const double yz = (yvals[i] - my) / sy;
    const double u =  xz*c + yz*s;
    const double v = -xz*s + yz*c;
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

  TH2D* hRot = new TH2D("hRot", Form("Scaled-PCA rotated coords, Run %d;u = scaled PC1 coordinate;v = scaled PC2 coordinate", nrun),
                        220, umin, umax, 220, vmin, vmax);
  TH1D* hU = new TH1D("hU", Form("Run %d;u = along scaled PC1;counts", nrun), 220, umin, umax);
  TH1D* hV = new TH1D("hV", Form("Run %d;v = across scaled PC1;counts", nrun), 220, vmin, vmax);

  for (Long64_t i=0; i<N; ++i) {
    hRot->Fill(uvals[i], vvals[i]);
    hU->Fill(uvals[i]);
    hV->Fill(vvals[i]);
  }

  TString outbase = Form("plots/yfp_ypfp_scaled_pc1_run%d_delta_%g_to_%g", nrun, deltaMin, deltaMax);
  outbase.ReplaceAll("-", "m");
  outbase.ReplaceAll(".", "p");
  TString outpdf = outbase + ".pdf";
  TString outroot = outbase + ".root";
  gSystem->mkdir("plots", kTRUE);

  TCanvas* c1 = new TCanvas("c_pc1", "PC1 diagnostic", 1200, 950);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetLogz();
  hOrig->Draw("COLZ");

  // Draw PC1 and PC2 through centroid. Use span based on plot range.
  double L = 0.55 * std::sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
  TLine* pc1 = new TLine(mx - L*dxRaw, my - L*dyRaw, mx + L*dxRaw, my + L*dyRaw);
  pc1->SetLineColor(kOrange+7);
  pc1->SetLineWidth(4);
  pc1->Draw("same");

  // Raw-space direction corresponding to scaled-space PC2=(-s,c).
  double dxRaw2 = -sx * s;
  double dyRaw2 =  sy * c;
  const double normRaw2 = std::sqrt(dxRaw2*dxRaw2 + dyRaw2*dyRaw2);
  dxRaw2 /= normRaw2;
  dyRaw2 /= normRaw2;
  TLine* pc2 = new TLine(mx - L*dxRaw2, my - L*dyRaw2, mx + L*dxRaw2, my + L*dyRaw2);
  pc2->SetLineColor(kMagenta+2);
  pc2->SetLineWidth(2);
  pc2->SetLineStyle(2);
  pc2->Draw("same");
  TMarker* cen = new TMarker(mx, my, 29);
  cen->SetMarkerColor(kBlack);
  cen->SetMarkerSize(1.8);
  cen->Draw("same");
  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.12,0.92,Form("Scaled PC1: %.2f deg in z-space; raw overlay %.2f deg", thetaDegScaled, thetaDegRaw));

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
  hRot->Write();
  hU->Write();
  hV->Write();
  fout.WriteObject(c1, "c_pc1");
  fout.Close();

  cout << "Wrote: " << outpdf << endl;
  cout << "Wrote: " << outroot << endl;
}
