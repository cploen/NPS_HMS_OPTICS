// make_elastic_overlays.C
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <vector>
#include <string>
#include <utility>   // for std::pair
#include <sstream>   // for std::ostringstream
#include <cmath>
#include <iostream>
#include <TStyle.h>

// ——————————————————————————————————————
// 1) Top‐level struct & settings
// ——————————————————————————————————————
struct KinSetting {
  Int_t                                      code;       // e.g. 6117, 6667, etc.
  std::string                                label;      // e.g. "6p117", "6p667"
  std::vector<std::vector<Int_t>>            runGroups;  // each sub‐vector is one logical dataset
  std::vector<double>                        angles;     // matching one entry per group
  std::pair<double,double>                   xBcut;      // {xB_min, xB_max}
  double				     Mp_peak;    // from simc simulation
};

static const std::vector<KinSetting> allSettings = {
  // 6.667 GeV, three singleton runs
  { 
    6667, "6p667",
    { {1534}, {1535}, {1536} },
    { 19.145, 17.995, 16.850 },
    { 0.983, 1.012 },
    { 0.944992}
  },
  // 5.878 GeV, three singleton runs
  { 
    5878, "5p878",
    { {1714}, {1715}, {1716} },
    { 22.835, 21.685, 20.545 },
    { 0.99, 1.01 },
    {0.946603}
  },
  // 6.117 GeV, five fragments → three groups
  { 
    6117, "6p117",
    { {1250}, {1251,1252}, {1253,1259} },
    { 19.260, 20.700, 22.120 },
    { 0.988, 1.011 },
    {0.944904}
  },
  // 5.639 GeV, three singleton runs
  { 
    5639, "5p639",
    { {1267}, {1268}, {1269} },
    { 21.270, 22.700, 24.300 },
    { 0.985, 1.015 },
    std::numeric_limits<double>::quiet_NaN()
  }
};

// ——————————————————————————————————————
// 2) The single macro you .x in ROOT
// ——————————————————————————————————————
void make_elastic_overlays(Int_t pRun=6667) {
  gStyle->SetOptStat(0);

  // — find the matching KinSetting —
  const KinSetting* ks = nullptr;
  for (const auto &s : allSettings) {
    if (s.code == pRun) { ks = &s; break; }
  }
  if (!ks) {
    std::cerr << "✖ Invalid pRun " << pRun << std::endl;
    return;
  }

  const int N = ks->runGroups.size();

  // — build human‐readable group names like "1251,1252" —
  std::vector<std::string> groupNames(N);
  for (int i = 0; i < N; ++i) {
    std::ostringstream oss;
    for (size_t j = 0; j < ks->runGroups[i].size(); ++j) {
      if (j) oss << ",";
      oss << ks->runGroups[i][j];
    }
    groupNames[i] = oss.str();
  }

  // — placeholders for your future M_p peaks —
//  std::vector<double> Mp_peak(N, std::nan(""));

  // — book your histograms —
  std::vector<TH1F*> hKinW(N), hDelta(N), hW2(N), hYptar(N), hXb_uncut(N);
  std::vector<TH2F*> h2_xptarW(N), h2_yptarW(N), h2_ytarW(N), h2_xfpW(N), h2_xpfpW(N) ;

for (int i = 0; i < N; ++i) {
  double ang = ks->angles[i];
  TString tag = Form("%s_grp%d", ks->label.c_str(), i);
  const char* gn = groupNames[i].c_str();

  hXb_uncut[i] = new TH1F(
    Form("hXbUncut_%s", tag.Data()),
    Form("Runs %s @%.3f#circ; x_{B}; Counts", gn, ang),
    100, 0.8, 1.2
  );

    hKinW[i] = new TH1F(
      Form("hKinW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W [GeV]; Counts", gn, ang),
      100, 0.8, 1.1
    );
    hDelta[i] = new TH1F(
      Form("hDelta_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; delta p/p [%%]; Counts", gn, ang),
      100, -15, 15
    );
    hW2[i] = new TH1F(
      Form("hW2_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W^{2} [GeV]^{2}; Counts", gn, ang),
      100, 0.5, 1.2
    );
    hYptar[i] = new TH1F(
      Form("hYptar_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; yptar [#circ]; Counts", gn, ang),
      100, -0.15, 0.15
    );

    h2_xptarW[i] = new TH2F(
      Form("h2_xptarW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W [GeV]; xptar [#circ]", gn, ang),
      100, 0.8, 1.2, 100, -0.1, 0.1
    );
    h2_yptarW[i] = new TH2F(
      Form("h2_yptarW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ;  W [GeV]; yptar [#circ]", gn, ang),
      100, 0.8, 1.2, 100, -0.15, 0.15
    );
    h2_ytarW[i] = new TH2F(
      Form("h2_ytarW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W [GeV]; y_{tar} [#circ]", gn, ang),
      100, 0.8, 1.2, 100, -3.0, 3.0
    );

    h2_xfpW[i] = new TH2F(
      Form("h2_xfpW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W [GeV]; x_{fp} [mm]", gn, ang),
      100, 0.8, 1.2, //W
      100, -60, 60   // xfp
    );

    h2_xpfpW[i] = new TH2F(
      Form("h2_xpfpW_%s", tag.Data()),
      Form("Runs %s @%.3f#circ; W [GeV]; x'_{fp} [mrad]", gn, ang),
      100, 0.8, 1.2,
      100, -0.06, 0.06   // or whatever your focal-plane angular range is
    );
  }


  // — loop over each group, build a TChain, set branches, fill histos —
  for (int i = 0; i < N; ++i) {
    TChain chain("T");
    for (auto run : ks->runGroups[i]) {
      chain.Add(Form("../ROOTfiles/OPTICS/nps_hms_optics_%d_1_-1.root", run));
    }

    // Branches
    Double_t kinW, xbcalc, yptar, delta, npeSum, etracknorm, xptar, xtar, ytar, xfp, xpfp;
    chain.SetBranchAddress("H.kin.W",        &kinW);
    chain.SetBranchAddress("H.kin.x_bj",     &xbcalc);
    chain.SetBranchAddress("H.gtr.ph",       &yptar);
    chain.SetBranchAddress("H.gtr.dp",       &delta);
    chain.SetBranchAddress("H.cal.etracknorm",&etracknorm);
    chain.SetBranchAddress("H.cer.npeSum",   &npeSum);
    chain.SetBranchAddress("H.gtr.th",	     &xptar);
    chain.SetBranchAddress("H.gtr.x", 	     &xtar);
    chain.SetBranchAddress("H.gtr.y", 	     &ytar);
    chain.SetBranchAddress("H.dc.x_fp",  &xfp);
    chain.SetBranchAddress("H.dc.xp_fp", &xpfp);

    Long64_t nentries = chain.GetEntries();
    double xbmin = ks->xBcut.first, xbmax = ks->xBcut.second;

    for (Long64_t ev = 0; ev < nentries; ++ev) {
      chain.GetEntry(ev);
      
      hXb_uncut[i]->Fill(xbcalc);

      if (npeSum > 6 && etracknorm > 0.65 
          && xbcalc > xbmin && xbcalc < xbmax) {
        hKinW[i]  ->Fill(kinW);
        hDelta[i] ->Fill(delta);
        hW2[i]    ->Fill(kinW*kinW);
        hYptar[i] ->Fill(yptar);

        h2_xptarW[i]->Fill(kinW,xptar);
        h2_yptarW[i]->Fill(kinW,yptar);
        h2_ytarW[i] ->Fill(kinW,ytar);
        h2_xfpW[i] ->Fill(kinW,xfp);
        h2_xpfpW[i]->Fill(kinW,xpfp);
      }
    }
  }

TH2F* h2_xptarW_all = (TH2F*)h2_xptarW[0]->Clone("h2_xptarW_all");
for (int i = 1; i < N; ++i) h2_xptarW_all->Add(h2_xptarW[i]);
h2_xptarW_all->SetTitle("Combined; W [GeV]; xptar [#circ]");

TH2F* h2_yptarW_all = (TH2F*)h2_yptarW[0]->Clone("h2_yptarW_all");
for (int i = 1; i < N; ++i) h2_yptarW_all->Add(h2_yptarW[i]);
h2_yptarW_all->SetTitle("Combined; W [GeV]; yptar [#circ]");

TH2F* h2_ytarW_all  = (TH2F*)h2_ytarW[0]->Clone("h2_ytarW_all");
for (int i = 1; i < N; ++i) h2_ytarW_all->Add(h2_ytarW[i]);
h2_ytarW_all->SetTitle("Combined; W [GeV]; y_{tar} [#circ]");

TH2F* all_xfpW   = (TH2F*)h2_xfpW[0]  ->Clone("all_xfpW");
TH2F* all_xpfpW  = (TH2F*)h2_xpfpW[0] ->Clone("all_xpfpW");
for(int i=1;i<N;++i) {
  all_xfpW  ->Add(h2_xfpW[i]);
  all_xpfpW ->Add(h2_xpfpW[i]);
}
all_xfpW ->SetTitle("Combined; W [GeV]; x_{fp} [mm]");
all_xpfpW->SetTitle("Combined; W [GeV]; x'_{fp} [mrad]");


// 3) After filling, draw the overlay with cut markers:
TCanvas cXb("XbUncut","Uncut x_{B} distribution",800,600);
cXb.SetTopMargin(0.12);
gStyle->SetOptStat(0);

// find global maximum
double maxY=0;
for (auto& h : hXb_uncut) {
  h->SetStats(kFALSE);
  maxY = std::max(maxY, h->GetMaximum());
}
hXb_uncut[0]->SetMaximum(maxY * 1.1);

TLegend leg(0.15,0.75,0.4,0.9);
leg.SetTextSize(0.03);

// overlay the runs
for (int i = 0; i < N; ++i) {
  hXb_uncut[i]->SetLineColor(i+1);
  hXb_uncut[i]->SetLineWidth(2);
  if (i == 0) {
    hXb_uncut[i]->Draw("hist");
    hXb_uncut[i]->GetXaxis()->SetTitle("x_{B}");
    hXb_uncut[i]->GetYaxis()->SetTitle("Counts");
  } else {
    hXb_uncut[i]->Draw("hist same");
  }
  leg.AddEntry(hXb_uncut[i],
               Form("Runs %s @%.3f#circ",groupNames[i].c_str(), ks->angles[i]),
               "l");
}
// draw the cut window
double xbmin = ks->xBcut.first,
       xbmax = ks->xBcut.second;
TLine lmin(xbmin, 0, xbmin, maxY*1.1),
      lmax(xbmax, 0, xbmax, maxY*1.1);
lmin.SetLineColor(kBlack); lmin.SetLineStyle(2);
lmax.SetLineColor(kBlack); lmax.SetLineStyle(2);
lmin.Draw(); lmax.Draw();
leg.AddEntry(&lmin, Form("x_{B} cut [%.3f,%.3f]", xbmin, xbmax), "l");

leg.Draw();
cXb.SaveAs(Form("%s_xB_uncut.pdf", ks->label.c_str()));

std::cout << "=== Fill summary for pRun=" << pRun << " ===\n";
for (int i = 0; i < N; ++i) {
  std::cout 
    << " Group " << i << " (" << groupNames[i] << " @"
    << ks->angles[i] << "°): "
    << "W entries="    << hKinW[i]->GetEntries() 
    << ", Δ entries="  << hDelta[i]->GetEntries() 
    << ", W2 entries=" << hW2[i]->GetEntries() 
    << ", yptar entries=" << hYptar[i]->GetEntries()
    << std::endl;
}

TCanvas c1("all_xptarW","x′ vs W (all runs)",800,600);
c1.SetTopMargin(0.12);
gStyle->SetOptStat(0);
h2_xptarW_all->Draw("COLZ");
c1.SaveAs(Form("%s_xptar_vs_W_all.pdf", ks->label.c_str()));

TCanvas c2("all_yptarW","x′ vs W (all runs)",800,600);
c2.SetTopMargin(0.12);
gStyle->SetOptStat(0);
h2_yptarW_all->Draw("COLZ");
c2.SaveAs(Form("%s_yptar_vs_W_all.pdf", ks->label.c_str()));

TCanvas c3("all_ytarW","x′ vs W (all runs)",800,600);
c3.SetTopMargin(0.12);
gStyle->SetOptStat(0);
h2_ytarW_all->Draw("COLZ");
c3.SaveAs(Form("%s_ytar_vs_W_all.pdf", ks->label.c_str()));

TCanvas c4("all_xfpW","x_{fp} vs W",800,600);
c4.SetTopMargin(0.12); gStyle->SetOptStat(0);
all_xfpW->Draw("COLZ");
c4.SaveAs(Form("%s_xfp_vs_W_all.pdf", ks->label.c_str()));

TCanvas c5("all_xpfpW","x'_{fp} vs W",800,600);
c5.SetTopMargin(0.12); gStyle->SetOptStat(0);
all_xpfpW->Draw("COLZ");
c5.SaveAs(Form("%s_xpfp_vs_W_all.pdf", ks->label.c_str()));

  // — write histograms to a ROOT file —
  TFile outfile(Form("%s_histos.root", ks->label.c_str()), "RECREATE");
  for (int i = 0; i < N; ++i) {
    hKinW[i]->Write();
    hDelta[i]->Write();
    hW2[i]   ->Write();
    hXb_uncut[i]->Write();
    hYptar[i]->Write();
    h2_xptarW[i]->Write();
    h2_yptarW[i]->Write();
    h2_ytarW[i] ->Write();
 }
  outfile.Close();

  // — overlay helper —
  auto makeOverlay = [&](auto &hists, const char *name, 
                         const char *xTitle, bool drawPeak){
    TCanvas c(name, name, 800,600);
    c.SetTopMargin(0.12);            // give a bit more headroom
    gStyle->SetOptStat(0);           // turn off stats box

    // 1) find the largest bin height across all histos
    double maxY = 0;
    for (auto &h : hists) {
      maxY = std::max(maxY, h->GetMaximum());
      h->SetStats(kFALSE);           // also disable per-histo stats
    }

    // 2) extend the y-axis on the _first_ histogram
    hists[0]->SetMaximum(maxY * 1.1);

    TLegend leg(0.15,0.75,0.4,0.9);
    leg.SetTextSize(0.03);

    for (int i = 0; i < N; ++i) {
      hists[i]->SetLineColor(i+1);
      hists[i]->SetLineWidth(2);
      if (i == 0) {
        hists[i]->Draw("hist");
        hists[i]->GetXaxis()->SetTitle(xTitle);
        hists[i]->GetYaxis()->SetTitle("Counts");
      } else {
        hists[i]->Draw("hist same");
      }
      leg.AddEntry(hists[i],
                   Form("Runs %s @%.3f #circ", 
                        groupNames[i].c_str(), ks->angles[i]),
                   "l");
      // optional Mp line
      if (drawPeak && std::isfinite(Mp_peak[i])) {
        TLine line(Mp_peak[i], gPad->GetUymin(),
                   Mp_peak[i], gPad->GetUymax());
        line.SetLineStyle(2);
        line.SetLineColor(i+1);
        line.Draw();
        leg.AddEntry(&line, 
                     Form("M_p (runs %s)", groupNames[i].c_str()),
                     "l");
      }
    }
    leg.Draw();
    c.SaveAs(Form("%s_%s.pdf", ks->label.c_str(), name));
  };

  // — produce all your overlays —
  makeOverlay(hKinW,  "KinW",  "W [GeV]",        true);
  makeOverlay(hDelta, "Delta", "δp/p [%]",       false);
  makeOverlay(hW2,    "W2",    "W^{2} [GeV^{2}]", true);
  makeOverlay(hYptar, "Yptar", "yptar [deg]",        false);

  // 2‐D overlays (no Mp line)
  {
  TCanvas c("xptar_vs_W","xptar vs W", 1200, 400);
  c.Divide(N,1, 0.01, 0.01);  // N pads side by side, tiny margins

  for(int i=0;i<N;++i){
    c.cd(i+1);
    h2_xptarW[i]->Draw("colz");
    TLatex t; 
    t.DrawLatexNDC(0.15,0.85, Form("Runs %s @%.3f#circ",groupNames[i].c_str(),ks->angles[i]));
  }
  c.SaveAs(Form("%s_xptar_vs_W.pdf", ks->label.c_str()));
  }
  {
  TCanvas c("yptar_vs_W","y′ vs W", 1200, 400);
  c.Divide(N,1, 0.01, 0.01);  // N pads side by side, tiny margins

  for(int i=0;i<N;++i){
    c.cd(i+1);
    h2_yptarW[i]->Draw("colz");
    TLatex t; 
    t.DrawLatexNDC(0.15,0.85, Form("Runs %s @%.3f#circ",groupNames[i].c_str(),ks->angles[i]));
  }
  c.SaveAs(Form("%s_yptar_vs_W.pdf", ks->label.c_str()));
  }
  {
 TCanvas c("ytar_vs_W_panels", "y_{tar} vs W (panels)", 1200, 400);
  c.Divide(N, 1, 0.01, 0.01);  // N panels horizontally

  for (int i = 0; i < N; ++i) {
    c.cd(i+1);
    h2_ytarW[i]->SetStats(kFALSE);
    h2_ytarW[i]->Draw("colz");
    TLatex label;
    label.DrawLatexNDC(0.15, 0.85,
      Form("Runs %s @%.3f#circ", groupNames[i].c_str(), ks->angles[i])
    );
  }
  c.SaveAs(Form("%s_ytar_vs_W_panels.pdf", ks->label.c_str()));
  }



}
