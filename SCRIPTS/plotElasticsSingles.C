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
#include <string>
#include <TLegend.h>
using namespace std;

void plotElasticsSingles(Int_t pRun=6117){
//4 digit momentum, no decimal, as arg:
//6667, 6117, 5878, 5369
//
//Readme: set runs, kinematics, and output file type in beginning of script

const int nSettings=3;
//double centAng;
//double centMom;
std::string MatrixEl="ang";
	cout << " Plotting runs for HMS_P = 5.878 GeV " << endl;
	//"KinC-x60-3" Runs 1714, 1715, 1716
	TFile *f = TFile::Open("./Output/elastic/hist/Singles_5p878ang.root");
	TFile *fout = new TFile("Output/elastic/hist/Overlays_5p878ang.root","RECREATE");
	//std::string kinSetting="KinC-x60-3"; 
	 std::string kinSetting="5p878";
 
	int nRun[nSettings]={1716,1715,1714};
	double centAng[nSettings]={20.545,21.685,22.835};//deg
	double centMom=-5.878;

/*
if (pRun == 6667) {
	cout << " Plotting runs for HMS_P = 6.667 GeV " << endl;
	//"KinC-x50-2" Runs 1534,1535,1536
	TFile *f = TFile::Open("./Output/elastic/hist/Singles_6p667ang.root");
	TFile *fout = new TFile("Output/elastic/hist/Overlays_6p667ang.root","RECREATE");
//	 std::string kinSetting="KinC-x50-2"; 
	 std::string kinSetting="6p667"; 

	int nRun[nSettings]={1534,1535,1536};
	double centAng[nSettings]={19.145,17.995,16.850};//deg
	double centMom=-6.667;
} else if (pRun == 6117){
	cout << " Plotting runs for HMS_P = 6.117 GeV " << endl;
	//KinC-x36-3 Runs 1250,1251,1259
	TFile *f = TFile::Open("./Output/elastic/hist/Singles_6p117ang.root");
	TFile *fout = new TFile("Output/elastic/hist/Overlays_6p117ang.root","RECREATE");
	//std::string kinSetting="KinC-x36-3"; 
	 std::string kinSetting="6p117"; 

	int nRun[nSettings]={1250,1251,1259}; 
	double centAng[nSettings]={19.26,20.69,22.12};//deg
	double centMom=-6.117;
} else if (pRun == 5878){
	cout << " Plotting runs for HMS_P = 5.878 GeV " << endl;
	//"KinC-x60-3" Runs 1714, 1715, 1716
	TFile *f = TFile::Open("./Output/elastic/hist/Singles_5p878ang.root");
	TFile *fout = new TFile("Output/elastic/hist/Overlays_5p878ang.root","RECREATE");
	//std::string kinSetting="KinC-x60-3"; 
	 std::string kinSetting="5p878";
 
	int nRun[nSettings]={1716,1715,1714};
	double centAng[nSettings]={20.545,21.685,22.835};//deg
	double centMom=-5.878;
} else if (pRun == 5639){
	cout << " Plotting runs for HMS_P = 5.639 GeV " << endl;
	//"KinC-x25-3'" Runs 1267,1268,1269
	TFile *f = TFile::Open("./Output/elastic/hist/Singles_5p639ang.root");
	TFile *fout = new TFile("Output/elastic/hist/Overlays_5p639ang.root","RECREATE");
	//std::string kinSetting="KinC-x25-3'"; 
	std::string kinSetting="5p639"; 

	int nRun[nSettings]={1267,1268,1269};
	double centAng[nSettings]={21.27,22.7,24.30};//deg
	double centMom=-5.639;
}else {
	cout << "Invalid momentum setting" << endl;
	return;
}
*/
std::string pdf_file_name=Form("Output/elastic/plots/Singles2D_");
std::string file_format=".pdf";

//make the output file
TCanvas *canvas = new TCanvas("canvas", "plots",800,800);
canvas->SetFillColor(0);
canvas->SetBorderMode(0);
canvas->SetBorderSize(0);
canvas->SetFrameFillColor(0);
canvas->SetFrameBorderMode(0);
gROOT->SetBatch(true);
gStyle->SetOptStat(0);
canvas->SetGridx();
canvas->SetGridy();

//save plots
canvas->Update();

//2D Histos Singles
TH2F *hWdp0 = (TH2F*)f->Get("h2_kinWdp_0");
  hWdp0->SetName("hWdp0");
  hWdp0->SetTitle(Form("h_kinWdp Run %d at %f GeV/c and %f deg",nRun[0],centMom,centAng[0] )); 
  hWdp0->SetLineColor(kRed);
  hWdp0->SetLineWidth(2);
  hWdp0->Scale(1./hWdp0->GetMaximum());
  hWdp0->Draw("histo");

/*
TH2F *hWdp1 = (TH2F*)f->Get("h_kinWdp_1");
//  hWdp1->SetName("hWdp1");
  hWdp1->SetTitle(Form("h_kinWdp Run %d at %f GeV/c and %f deg",nRun[1],centMom,centAng[1] )); 
  hWdp1->SetLineColor(kGreen);
  hWdp1->SetLineWidth(2);
  //hW22->Scale(1./hW22->GetMaximum());
  hWdp1->Draw("histo");

  TLegend *legL = new TLegend(0.125,0.875,0.4,0.75);//Top left, quite good
//TLegend *legL = new TLegend(0.15,0.85,0.4,0.7);//top left, not bad
//TLegend *legL = new TLegend(0.1,0.1,0.4,0.2);//bottom left
//TLegend *legR = new TLegend(0.5,0.6,0.8,0.8);//upper right
  legL->SetBorderSize(0); //no border
  legL->SetHeader("HMS Theta (deg)","C"); // option "C" allows to center the header
  legL->AddEntry(hWdp0,Form("Run %d CentAngle %f",nRun[0], centAng[0]),"l");
  legL->AddEntry(hWdp1,Form("Run %d CentAngle %f",nRun[1], centAng[1]),"l");
  //legL->AddEntry(hWdp2,Form("Run %d CentAngle %f",nRun[2], centAng[2]),"l");
  legL->Draw();
*/
//canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+"(").c_str());

canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+")").c_str());

fout->Write();
fout->Close();


}
