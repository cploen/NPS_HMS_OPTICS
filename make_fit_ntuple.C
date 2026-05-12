#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void make_fit_ntuple(Int_t nrun=1814,Int_t FileID=-2){
 Bool_t CutYtarFlag=kTRUE;
 Bool_t CutYpFpYFpFlag=kTRUE;
 Bool_t CutXpFpXFpFlag=kTRUE;

 gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);

 //  Get info for that optics run
 TString OpticsFile = Form("DATfiles/list_of_optics_run.dat");
   ifstream file_optics(OpticsFile.Data());
 TString opticsline;
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Double_t ymis =0.0;
  Int_t NumFoil=0;
  TString temp;

  vector <Double_t> ztar_foil;
  Int_t ndelcut;
  vector<Double_t > delcut;

  if (file_optics.is_open()) {
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      ndelcut = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ymis = temp.Atof();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
        temp.ReadToDelim(file_optics);
	ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut;nd++) {
        temp.ReadToDelim(file_optics,',');
	delcut.push_back(temp.Atof());
      }
        temp.ReadToDelim(file_optics);
	delcut.push_back(temp.Atof());
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }

  Double_t     y_mis;
  Double_t     x_mis;
 
  if (TMath::Abs(CentAngle)<40) {y_mis = 0.1*(0.52-0.012*TMath::Abs(CentAngle)+0.002*TMath::Abs(CentAngle)*TMath::Abs(CentAngle));} // cm
  else{y_mis = 0.1*(0.52-0.012*40. + 0.002*40.*40.);} // cm
 
  if (TMath::Abs(CentAngle)<50) {x_mis = 0.1*(2.37-0.086*TMath::Abs(CentAngle)+0.0012*TMath::Abs(CentAngle)*TMath::Abs(CentAngle));}
  else{x_mis = 0.1*(2.37-0.086*50.+0.0012*50.*50.);}
 
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << " y_mis = " << y_mis<< " and x_mis = " << x_mis<< endl;
  if (NumFoil==0) return;
  for (Int_t nf=0;nf<NumFoil;nf++) {
    cout << nf << " foil = " << ztar_foil[nf] << endl;
  }
  
  vector <Double_t> ys_cent;
  vector <Double_t> xs_cent;
  for (Int_t nys=0;nys<9;nys++) {
    Double_t ypos=(nys-4)*0.6*2.54;
    ys_cent.push_back(ypos);
    Double_t xpos=(nys-4)*2.54;
    xs_cent.push_back(xpos);
  }
 
   TString inputroot;
   TString outputroot;
//   inputroot=Form("ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_%s_1_%d.root", OpticsID.Data(),FileID); 
  
   inputroot=Form("ROOTfiles/OPTICS/nps_hms_optics_%s_1_%d.root", OpticsID.Data(),FileID); 

 outputroot= Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
  
	TH1F *hxbpm_tar = new TH1F("hxbpm_tar",Form("Run %d ; Xbpm_tar ; Counts",nrun),100,-2.,2.);
	TH1F *hybpm_tar = new TH1F("hybpm_tar",Form("Run %d ; Ybpm_tar ; Counts",nrun),100,-2.,2.);
 

  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile=Form("cuts/ytar_delta_%s_%d_cut.root",OpticsID.Data(),FileID);
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
 

   cout << " Cut file = " << YtarDeltaCutFile << endl;
   for (Int_t nc=0;nc<NumFoil;nc++) {
      fYtarDeltaCut->cd();
      TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
      if (tempcut) {
      Int_t npt = tempcut->GetN();
      cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
      ytar_delta_cut.push_back(tempcut);
      } else {
      cout << " No hYtarDelta_cut = " << nc << endl;
      }
   }
  }
 
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_yfp_cut;
  vector<vector<vector<Int_t> > > ypfp_yfp_cut_flag;
  ypfp_yfp_cut.resize(NumFoil);
  ypfp_yfp_cut_flag.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
          ypfp_yfp_cut[nf].resize(ndelcut);
          ypfp_yfp_cut_flag[nf].resize(ndelcut);
   }
  if (CutYpFpYFpFlag) {
   outCutFile=Form("cuts/YpFpYFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<9;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
      Int_t npt = tempg->GetN();
      //cout << "hYpFpYFp_cut = " << nf << " " << nd << " " << nc << " npts = " << npt << endl;
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	  ypfp_yfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  ypfp_yfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }
//
  TString xpfp_xfp_outCutFile;
  TFile *xpfp_xfp_fcut;
  vector<vector<vector<TCutG*> > > xpfp_xfp_cut;
  vector<vector<vector<Int_t> > > xpfp_xfp_cut_flag;
  xpfp_xfp_cut.resize(NumFoil);
  xpfp_xfp_cut_flag.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
          xpfp_xfp_cut[nf].resize(ndelcut);
          xpfp_xfp_cut_flag[nf].resize(ndelcut);
   }
  if (CutXpFpXFpFlag) {
    xpfp_xfp_outCutFile=Form("cuts/XpFpXFp_%s_%d_cut.root",OpticsID.Data(),FileID);
    xpfp_xfp_fcut = new TFile(xpfp_xfp_outCutFile);
    cout << "xpfp_xfp_ Cut file = " << xpfp_xfp_outCutFile << endl;
//added cout lines below
cout << "number of delta regions = " << ndelcut << endl;

    xpfp_xfp_fcut->cd();
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<9;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hXpFpXFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
	  xpfp_xfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hXpFpXFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  xpfp_xfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }

TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);
 Double_t  xbpm_tar;
   tsimc->SetBranchAddress("H.rb.raster.fr_xbpm_tar",&xbpm_tar);
 Double_t  ybpm_tar;
   tsimc->SetBranchAddress("H.rb.raster.fr_ybpm_tar",&ybpm_tar);
 Double_t frx;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&frx);
 Double_t fry;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&fry);

 Double_t xptarT,ytarT,yptarT,ysieveT,xsieveT,ztarT,ztar,xtarT;
   
   TFile hroot(outputroot,"recreate");
   TTree *otree = new TTree("TFit","FitTree");
   otree->Branch("ys",&ysieve);
   otree->Branch("ysT",&ysieveT);
   otree->Branch("xs",&xsieve);
   otree->Branch("xsT",&xsieveT);
   otree->Branch("ztarT",&ztarT);
   otree->Branch("ztar",&ztar);
   otree->Branch("xtar",&xtar);
   otree->Branch("xtarT",&xtarT);
   otree->Branch("xptar",&xptar);
   otree->Branch("yptar",&yptar);
   otree->Branch("ytar",&ytar);
   otree->Branch("xptarT",&xptarT);
   otree->Branch("yptarT",&yptarT);
   otree->Branch("ytarT",&ytarT);
   otree->Branch("delta",&delta);
   otree->Branch("xpfp",&xpfp);
   otree->Branch("ypfp",&ypfp);
   otree->Branch("xfp",&xfp);
   otree->Branch("yfp",&yfp);
   otree->Branch("reactxcalc",&reactx);
   otree->Branch("reactycalc",&reacty);

 Double_t zdis_sieve = 168.;
 
// loop over entries
 Long64_t nentries = tsimc->GetEntries();
 cout << " start loop for get beam positions " << nentries << endl;
 for (int i = 0; i < nentries; i++) {
   tsimc->GetEntry(i);
   if (i%50000==0) cout << " Entry = " << i << endl;
   if (etracknorm>.8 && sumnpe > 6. && delta>-10 && delta<10) {
     hxbpm_tar->Fill(xbpm_tar);
     hybpm_tar->Fill(ybpm_tar);		  
   }}

//Double_t xbeam = -hxbpm_tar->GetMean(); // OLD WAY: constant beam in Hall coordinates
//Double_t ybeam = hybpm_tar->GetMean();

//loop over entries
cout << " start loop " << nentries << endl;
CentAngle=CentAngle*3.14159/180.;
for (int i = 0; i < nentries; i++) {
  tsimc->GetEntry(i);
  Double_t xbeam = -xbpm_tar;
  Double_t ybeam =  ybpm_tar;

  if (i%50000==0) cout << " Entry = " << i << endl;
  if (etracknorm>.8 && sumnpe > 6. && delta>-10 && delta<10) {
    Int_t nf_found=-1, nd_found=-1,ny_found=-1,nx_found=-1;
    for  (UInt_t nf=0;nf<ytar_delta_cut.size();nf++) {
      if (ytar_delta_cut[nf]->IsInside(ytar,delta)) nf_found=nf;
    } 
    for  (UInt_t nd=0;nd<ndelcut;nd++) {
     if ( delta >=delcut[nd] && delta <delcut[nd+1])  nd_found=nd;
    }
    if (nf_found!=-1 && nd_found!=-1) {
      for  (UInt_t ny=0;ny<9;ny++) {
        if (ypfp_yfp_cut[nf_found][nd_found][ny] && ypfp_yfp_cut[nf_found][nd_found][ny]->IsInside(ypfp,yfp)) ny_found=ny;
	}
        for  (UInt_t nx=0;nx<9;nx++) {
	  if (xpfp_xfp_cut[nf_found][nd_found][nx] && xpfp_xfp_cut[nf_found][nd_found][nx]->IsInside(xpfp,xfp)) nx_found=nx;
	  }
    }
    if (nf_found !=-1 && nd_found!=-1 && ny_found!=-1 && nx_found!=-1) {
      Double_t ytar_cent = ztar_foil[nf_found]*TMath::Sin(CentAngle) + xbeam*cos(CentAngle)-y_mis;
      yptarT = (ys_cent[ny_found]-ytar_cent)/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));
      ytarT = +ztar_foil[nf_found]*(TMath::Sin(CentAngle)-yptarT*TMath::Cos(CentAngle)) +xbeam*(TMath::Cos(CentAngle)+yptarT*TMath::Sin(CentAngle))- y_mis;
      xptarT = (xs_cent[nx_found])/(zdis_sieve-ztar_foil[nf_found]*TMath::Cos(CentAngle));

      xtarT = -reacty - x_mis - xptarT*ztar_foil[nf_found]*TMath::Cos(CentAngle);
      ysieveT=ys_cent[ny_found];
      xsieveT=xs_cent[nx_found];
      ztarT=ztar_foil[nf_found];
      ztar=reactz;
      otree->Fill();
     }
   }
}

otree->Write();
}
