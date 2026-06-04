// dbscan_fp_pfp_islands.C
// Diagnostic DBSCAN clustering for HMS optics FP/PFP islands.
//
// Goal:
//   Given a replay ROOT file, a foil TCutG in (delta,ytar), and a delta slice,
//   test whether DBSCAN can separate the FP/PFP islands without angle scanning.
//
// Example:
//   hcana -l -q 'dbscan_fp_pfp_islands.C("/path/to/replay.root","ytar_delta_1544_multifoil_test1_multifoil_cut.root","cut_foil0",-1.0,1.0,"Y",0.22,12,50000,"dbscan_test")'
//
// Plane options:
//   "Y" -> cluster in (H.dc.y_fp,  H.dc.yp_fp)
//   "X" -> cluster in (H.dc.x_fp,  H.dc.xp_fp)
//
// Notes:
//   - Coordinates are standardized before DBSCAN, so eps is dimensionless.
//   - This is diagnostic code: tune eps/minPts and inspect the plots.
//   - Output includes a ROOT tree with event index, cluster id, FP/PFP coords, sieve coords.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TGraph.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TRandom3.h>

#include <vector>
#include <queue>
#include <map>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

struct DBPoint {
  Long64_t entry = -1;
  double raw1 = 0.0;       // yfp or xfp
  double raw2 = 0.0;       // ypfp or xpfp
  double z1 = 0.0;         // standardized raw1
  double z2 = 0.0;         // standardized raw2
  double delta = 0.0;
  double ytar = 0.0;
  double xsieve = 0.0;
  double ysieve = 0.0;
  int cluster = -99;       // -99 unvisited, -1 noise, >=0 cluster id
};

static double SafeStd(const std::vector<double>& v, double mean) {
  if (v.size() < 2) return 1.0;
  double s2 = 0.0;
  for (double x : v) s2 += (x - mean) * (x - mean);
  double s = std::sqrt(s2 / double(v.size() - 1));
  if (!std::isfinite(s) || s <= 1e-12) return 1.0;
  return s;
}

static std::vector<int> Palette() {
  return {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kCyan+2, kOrange+7,
          kViolet+1, kSpring+5, kAzure+6, kPink+7, kTeal+3, kYellow+2,
          kRed-4, kBlue-4, kGreen-6, kMagenta-6, kOrange-3, kCyan-6};
}

static TString Sanitize(TString s) {
  s.ReplaceAll("/","_"); s.ReplaceAll(" ","_"); s.ReplaceAll(".","p");
  s.ReplaceAll("-","m"); s.ReplaceAll("+","p");
  return s;
}

// Cell-grid accelerated radius query in standardized coordinates.
class NeighborGrid {
public:
  NeighborGrid(const std::vector<DBPoint>& pts, double eps_) : p(pts), eps(eps_) {
    cellSize = eps;
    if (cellSize <= 0) cellSize = 0.2;
    for (int i = 0; i < (int)p.size(); ++i) {
      auto key = CellKey(p[i].z1, p[i].z2);
      grid[key].push_back(i);
    }
  }

  std::vector<int> RegionQuery(int idx) const {
    std::vector<int> out;
    auto c = CellKey(p[idx].z1, p[idx].z2);
    const double eps2 = eps * eps;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        auto it = grid.find(std::make_pair(c.first + dx, c.second + dy));
        if (it == grid.end()) continue;
        for (int j : it->second) {
          double d1 = p[j].z1 - p[idx].z1;
          double d2 = p[j].z2 - p[idx].z2;
          if (d1*d1 + d2*d2 <= eps2) out.push_back(j);
        }
      }
    }
    return out;
  }

private:
  const std::vector<DBPoint>& p;
  double eps = 0.2;
  double cellSize = 0.2;
  std::map<std::pair<int,int>, std::vector<int>> grid;

  std::pair<int,int> CellKey(double x, double y) const {
    return std::make_pair((int)std::floor(x / cellSize), (int)std::floor(y / cellSize));
  }
};

static int RunDBSCAN(std::vector<DBPoint>& pts, double eps, int minPts) {
  NeighborGrid ng(pts, eps);
  int clusterId = 0;

  for (int i = 0; i < (int)pts.size(); ++i) {
    if (pts[i].cluster != -99) continue;

    std::vector<int> neigh = ng.RegionQuery(i);
    if ((int)neigh.size() < minPts) {
      pts[i].cluster = -1; // provisional noise; may be absorbed later
      continue;
    }

    pts[i].cluster = clusterId;
    std::queue<int> q;
    for (int j : neigh) if (j != i) q.push(j);

    while (!q.empty()) {
      int j = q.front(); q.pop();

      if (pts[j].cluster == -1) pts[j].cluster = clusterId;
      if (pts[j].cluster != -99) continue;

      pts[j].cluster = clusterId;
      std::vector<int> neigh2 = ng.RegionQuery(j);
      if ((int)neigh2.size() >= minPts) {
        for (int k : neigh2) {
          if (pts[k].cluster == -99 || pts[k].cluster == -1) q.push(k);
        }
      }
    }

    clusterId++;
  }

  return clusterId;
}

void dbscan_fp_pfp_islands(
  const char* inputRoot = "replay.root",
  const char* cutRoot   = "ytar_delta_cut.root",
  const char* foilCutName = "delta_vs_ytar_cut_foil0",
  double deltaMin = -1.0,
  double deltaMax =  1.0,
  const char* plane = "Y",
  double eps = 0.22,
  int minPts = 12,
  int maxEvents = 50000,
  const char* outTag = "dbscan_fp_pfp"
) {
  gStyle->SetOptStat(0);

  const bool useY = (TString(plane).BeginsWith("Y", TString::kIgnoreCase));
  const TString coord1Name = useY ? "H.dc.y_fp"  : "H.dc.x_fp";
  const TString coord2Name = useY ? "H.dc.yp_fp" : "H.dc.xp_fp";
  const TString planeTag = useY ? "Y" : "X";

  TFile* fin = TFile::Open(inputRoot, "READ");
  if (!fin || fin->IsZombie()) { std::cerr << "ERROR: cannot open " << inputRoot << "\n"; return; }
  TTree* T = (TTree*)fin->Get("T");
  if (!T) { std::cerr << "ERROR: no tree named T in " << inputRoot << "\n"; return; }

  TCutG* foilCut = nullptr;
  TFile* fcut = TFile::Open(cutRoot, "READ");
  if (fcut && !fcut->IsZombie()) {
    foilCut = (TCutG*)fcut->Get(foilCutName);
    if (!foilCut) std::cerr << "WARNING: no TCutG named " << foilCutName << " in " << cutRoot << "; using delta slice only.\n";
  } else {
    std::cerr << "WARNING: cannot open cut file " << cutRoot << "; using delta slice only.\n";
  }

  double cer=0, cal=0, ytar=0, delta=0, yfp=0, ypfp=0, xfp=0, xpfp=0, xs=0, ys=0;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("H.cer.npeSum",1);           T->SetBranchAddress("H.cer.npeSum", &cer);
  T->SetBranchStatus("H.cal.etottracknorm",1);    T->SetBranchAddress("H.cal.etottracknorm", &cal);
  T->SetBranchStatus("H.gtr.y",1);                T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchStatus("H.gtr.dp",1);               T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchStatus("H.dc.y_fp",1);              T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchStatus("H.dc.yp_fp",1);             T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchStatus("H.dc.x_fp",1);              T->SetBranchAddress("H.dc.x_fp", &xfp);
  T->SetBranchStatus("H.dc.xp_fp",1);             T->SetBranchAddress("H.dc.xp_fp", &xpfp);
  T->SetBranchStatus("H.extcor.xsieve",1);        T->SetBranchAddress("H.extcor.xsieve", &xs);
  T->SetBranchStatus("H.extcor.ysieve",1);        T->SetBranchAddress("H.extcor.ysieve", &ys);

  std::vector<DBPoint> pts;
  std::vector<double> v1, v2;

  Long64_t nEnt = T->GetEntries();
  Long64_t stride = 1;
  if (maxEvents > 0 && nEnt > maxEvents*3LL) stride = std::max<Long64_t>(1, nEnt / maxEvents);

  for (Long64_t i = 0; i < nEnt; i += stride) {
    T->GetEntry(i);

    // Baseline electron-ish cuts. Loosen/remove here if needed.
    if (cer <= 2.0) continue;
    if (cal <= 0.6) continue;
    if (delta < deltaMin || delta > deltaMax) continue;
    if (foilCut && !foilCut->IsInside(ytar,delta)) continue; // assumes TCutG axes are delta vs ytar

    DBPoint p;
    p.entry = i;
    p.raw1 = useY ? yfp : xfp;
    p.raw2 = useY ? ypfp : xpfp;
    p.delta = delta;
    p.ytar = ytar;
    p.xsieve = xs;
    p.ysieve = ys;
    if (!std::isfinite(p.raw1) || !std::isfinite(p.raw2)) continue;
    pts.push_back(p);
    v1.push_back(p.raw1);
    v2.push_back(p.raw2);
    if (maxEvents > 0 && (int)pts.size() >= maxEvents) break;
  }

  if (pts.size() < 10) {
    std::cerr << "ERROR: only selected " << pts.size() << " events. Check cuts/delta slice/cut axes.\n";
    return;
  }

  double mean1 = 0, mean2 = 0;
  for (double x : v1) mean1 += x;
  for (double x : v2) mean2 += x;
  mean1 /= v1.size(); mean2 /= v2.size();
  double sig1 = SafeStd(v1, mean1);
  double sig2 = SafeStd(v2, mean2);

  for (auto& p : pts) {
    p.z1 = (p.raw1 - mean1) / sig1;
    p.z2 = (p.raw2 - mean2) / sig2;
  }

  int nClusters = RunDBSCAN(pts, eps, minPts);
  std::vector<int> counts(nClusters, 0);
  int nNoise = 0;
  for (auto& p : pts) {
    if (p.cluster >= 0 && p.cluster < nClusters) counts[p.cluster]++;
    else if (p.cluster == -1) nNoise++;
  }

  std::cout << "\nDBSCAN FP/PFP diagnostic\n";
  std::cout << "  inputRoot   = " << inputRoot << "\n";
  std::cout << "  foilCut     = " << foilCutName << " from " << cutRoot << "\n";
  std::cout << "  plane       = " << planeTag << " using (" << coord1Name << ", " << coord2Name << ")\n";
  std::cout << "  delta       = [" << deltaMin << ", " << deltaMax << "]\n";
  std::cout << "  selected    = " << pts.size() << " events\n";
  std::cout << "  eps/minPts  = " << eps << " / " << minPts << "\n";
  std::cout << "  clusters    = " << nClusters << ", noise = " << nNoise << "\n";
  for (int c = 0; c < nClusters; ++c) std::cout << "    cluster " << std::setw(2) << c << ": " << counts[c] << "\n";

  TString tag = TString(outTag) + "_" + planeTag + "_" + Sanitize(foilCutName) + Form("_dm%.2f_dp%.2f", deltaMin, deltaMax);
  TString pdfName = tag + ".pdf";
  TString rootName = tag + ".root";

  double r1min=*std::min_element(v1.begin(),v1.end()), r1max=*std::max_element(v1.begin(),v1.end());
  double r2min=*std::min_element(v2.begin(),v2.end()), r2max=*std::max_element(v2.begin(),v2.end());
  double pad1 = 0.08*(r1max-r1min), pad2 = 0.08*(r2max-r2min);
  r1min-=pad1; r1max+=pad1; r2min-=pad2; r2max+=pad2;

  TH2D* hRaw = new TH2D("hRaw", Form("Selected events; %s; %s", coord1Name.Data(), coord2Name.Data()), 240,r1min,r1max,240,r2min,r2max);
  TH2D* hNoise = new TH2D("hNoise", Form("DBSCAN noise; %s; %s", coord1Name.Data(), coord2Name.Data()), 240,r1min,r1max,240,r2min,r2max);
  for (auto& p : pts) {
    hRaw->Fill(p.raw1, p.raw2);
    if (p.cluster < 0) hNoise->Fill(p.raw1, p.raw2);
  }

  std::vector<int> pal = Palette();
  std::vector<TGraph*> gFP(nClusters, nullptr), gSieve(nClusters, nullptr);
  std::vector<std::vector<double>> gx(nClusters), gy(nClusters), gsx(nClusters), gsy(nClusters);
  std::vector<double> nx, ny, nsx, nsy;

  for (auto& p : pts) {
    if (p.cluster >= 0) {
      gx[p.cluster].push_back(p.raw1); gy[p.cluster].push_back(p.raw2);
      gsx[p.cluster].push_back(p.xsieve); gsy[p.cluster].push_back(p.ysieve);
    } else {
      nx.push_back(p.raw1); ny.push_back(p.raw2); nsx.push_back(p.xsieve); nsy.push_back(p.ysieve);
    }
  }
  for (int c = 0; c < nClusters; ++c) {
    gFP[c] = new TGraph(gx[c].size(), gx[c].data(), gy[c].data());
    gFP[c]->SetMarkerStyle(20); gFP[c]->SetMarkerSize(0.35); gFP[c]->SetMarkerColor(pal[c % pal.size()]);
    gSieve[c] = new TGraph(gsx[c].size(), gsx[c].data(), gsy[c].data());
    gSieve[c]->SetMarkerStyle(20); gSieve[c]->SetMarkerSize(0.35); gSieve[c]->SetMarkerColor(pal[c % pal.size()]);
  }
  TGraph* gNoiseFP = new TGraph(nx.size(), nx.data(), ny.data());
  gNoiseFP->SetMarkerStyle(6); gNoiseFP->SetMarkerColor(kGray+1); gNoiseFP->SetMarkerSize(0.3);
  TGraph* gNoiseSieve = new TGraph(nsx.size(), nsx.data(), nsy.data());
  gNoiseSieve->SetMarkerStyle(6); gNoiseSieve->SetMarkerColor(kGray+1); gNoiseSieve->SetMarkerSize(0.3);

  TCanvas* c = new TCanvas("c", "dbscan fp/pfp", 1100, 850);
  c->Print(pdfName + "[");

  hRaw->Draw("COLZ");
  c->SetLogz(1);
  c->Print(pdfName);
  c->SetLogz(0);

  hRaw->Draw("AXIS");
  hRaw->SetTitle(Form("DBSCAN clusters: plane %s, eps=%.3f, minPts=%d, clusters=%d, noise=%d; %s; %s",
                      planeTag.Data(), eps, minPts, nClusters, nNoise, coord1Name.Data(), coord2Name.Data()));
  gNoiseFP->Draw("P SAME");
  for (int cc = 0; cc < nClusters; ++cc) gFP[cc]->Draw("P SAME");
  c->Print(pdfName);

  TH2D* hSieveAxis = new TH2D("hSieveAxis", "Same DBSCAN labels projected to sieve plane; H.extcor.xsieve; H.extcor.ysieve", 100,-15,15,100,-15,15);
  hSieveAxis->Draw("AXIS");
  gNoiseSieve->Draw("P SAME");
  for (int cc = 0; cc < nClusters; ++cc) gSieve[cc]->Draw("P SAME");
  c->Print(pdfName);

  TH1D* hCounts = new TH1D("hCounts", "DBSCAN cluster sizes; cluster id; events", std::max(1,nClusters), -0.5, nClusters-0.5);
  for (int cc = 0; cc < nClusters; ++cc) hCounts->SetBinContent(cc+1, counts[cc]);
  hCounts->Draw("HIST TEXT0");
  c->Print(pdfName);

  c->Print(pdfName + "]");

  TFile* fout = TFile::Open(rootName, "RECREATE");
  TTree* tout = new TTree("Tdbscan", "DBSCAN diagnostic selected events");
  Long64_t oentry; int ocluster; double oraw1, oraw2, oz1, oz2, odelta, oytar, oxs, oys;
  tout->Branch("entry", &oentry, "entry/L");
  tout->Branch("cluster", &ocluster, "cluster/I");
  tout->Branch("fp", &oraw1, "fp/D");
  tout->Branch("pfp", &oraw2, "pfp/D");
  tout->Branch("z_fp", &oz1, "z_fp/D");
  tout->Branch("z_pfp", &oz2, "z_pfp/D");
  tout->Branch("delta", &odelta, "delta/D");
  tout->Branch("ytar", &oytar, "ytar/D");
  tout->Branch("xsieve", &oxs, "xsieve/D");
  tout->Branch("ysieve", &oys, "ysieve/D");
  for (auto& p : pts) {
    oentry=p.entry; ocluster=p.cluster; oraw1=p.raw1; oraw2=p.raw2; oz1=p.z1; oz2=p.z2;
    odelta=p.delta; oytar=p.ytar; oxs=p.xsieve; oys=p.ysieve;
    tout->Fill();
  }
  TObjString meta(Form("plane=%s eps=%.6g minPts=%d deltaMin=%.6g deltaMax=%.6g input=%s cutFile=%s cutName=%s mean1=%.8g sig1=%.8g mean2=%.8g sig2=%.8g",
                       planeTag.Data(), eps, minPts, deltaMin, deltaMax, inputRoot, cutRoot, foilCutName, mean1, sig1, mean2, sig2));
  meta.Write("metadata");
  tout->Write();
  hRaw->Write();
  hNoise->Write();
  hCounts->Write();
  fout->Close();

  std::cout << "\nWrote:\n  " << pdfName << "\n  " << rootName << "\n\n";
}

