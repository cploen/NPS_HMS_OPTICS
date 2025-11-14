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

void plotElasticsOverlays(Int_t pRun=6117){
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
std::string pdf_file_name=Form("Output/elastic/plots/Overlays_");
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

//1D Histos
TH1F *hW20 = (TH1F*)f->Get("h_kinW2_0");
  hW20->SetName("hW20");
  hW20->SetTitle(Form("h_kinW2 at %f GeV/c",centMom)); 
  hW20->SetLineColor(kRed);
  hW20->SetLineWidth(2);
  //hW20->Scale(1./hW20->GetMaximum());
  hW20->Draw("histo");


TH1F *hW22 = (TH1F*)f->Get("h_kinW2_2");
  hW22->SetName("hW22");
  hW22->SetLineColor(kGreen);
  hW22->SetLineWidth(2);
  //hW22->Scale(1./hW22->GetMaximum());
  hW22->Draw("histosame");

TH1F *hW21 = (TH1F*)f->Get("h_kinW2_1");
  hW21->SetName("hW21");
  hW21->SetLineColor(kBlue);
  hW21->SetLineWidth(2);
  //hW21->Scale(1./hW21->GetMaximum());
  hW21->Draw("histosame");


  TLegend *legL = new TLegend(0.125,0.875,0.4,0.75);//Top left, quite good
//TLegend *legL = new TLegend(0.15,0.85,0.4,0.7);//top left, not bad
//TLegend *legL = new TLegend(0.1,0.1,0.4,0.2);//bottom left
//TLegend *legR = new TLegend(0.5,0.6,0.8,0.8);//upper right
  legL->SetBorderSize(0); //no border
  legL->SetHeader("HMS Theta (deg)","C"); // option "C" allows to center the header
  legL->AddEntry(hW20,Form("Run %d CentAngle %f",nRun[0], centAng[0]),"l");
  legL->AddEntry(hW21,Form("Run %d CentAngle %f",nRun[1], centAng[1]),"l");
  legL->AddEntry(hW22,Form("Run %d CentAngle %f",nRun[2], centAng[2]),"l");
  legL->Draw();

canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+"(").c_str());

TH1F *h0 = (TH1F*)f->Get("h_kinW_0");
  h0->SetName("h0");
 // h0->SetTitle("h_kinW %s",kinSetting.c_str());
  h0->SetTitle(Form("h_kinW at %f GeV/c",centMom));
  h0->SetLineColor(kRed);
  h0->SetLineWidth(2);
 // h0->Scale(1./h0->GetMaximum()); //normalize nevents
  h0->Draw("histo");

TH1F *h2 = (TH1F*)f->Get("h_kinW_2");
  h2->SetName("h2");
  h2->SetLineColor(kGreen);
  h2->SetLineWidth(2);
 // h2->Scale(1./h2->GetMaximum());
  h2->Draw("histosame");


TH1F *h1 = (TH1F*)f->Get("h_kinW_1");
  h1->SetName("h1");
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);
 // h2->Scale(1./h2->GetMaximum());
  h1->Draw("histosame");



  TLegend *legR = new TLegend(0.6,0.8,.89,.89);//This one's perfect! Very top RHS (x1,y1,x2,y2_)
//TLegend *legR = new TLegend(0.5,0.6,0.8,0.8);//top right
  legR->SetBorderSize(0); //no border
  legR->SetHeader("HMS Theta (deg)","C"); // option "C" allows to center the header
  legR->AddEntry(hW20,Form("Run %d CentAngle %f",nRun[0], centAng[0]),"l");
  legR->AddEntry(hW21,Form("Run %d CentAngle %f",nRun[1], centAng[1]),"l");
  legR->AddEntry(hW22,Form("Run %d CentAngle %f",nRun[2], centAng[2]),"l");
  legR->Draw();

canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+"(").c_str());
TH1F *hd0 = (TH1F*)f->Get("h_delta_0");
  hd0->SetName("hd0");
  hd0->SetTitle("h_delta");
  hd0->SetLineColor(kRed);
  hd0->SetLineWidth(2);
 // hd0->Scale(1./hd0->GetMaximum());
  hd0->Draw("histo");

TH1F *hd2 = (TH1F*)f->Get("h_delta_2");
  hd2->SetName("hd2");
  hd2->SetLineColor(kGreen);
  hd2->SetLineWidth(2);
 // hd2->Scale(1./hd2->GetMaximum());
  hd2->Draw("histosame");


TH1F *hd1 = (TH1F*)f->Get("h_delta_1");
  hd1->SetName("hd1");
  hd1->SetLineColor(kBlue);
  hd1->SetLineWidth(2);
 // hd1->Scale(1./hd1->GetMaximum());
  hd1->Draw("histosame");


// leg->Draw();
 canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+"(").c_str());
TH1F *hyptar0 = (TH1F*)f->Get("h_yptar_0");
  hyptar0->SetName("hyptar0");
 // h0->SetTitle("h_kinW %s",kinSetting.c_str());
  hyptar0->SetTitle("h_yptar");
  hyptar0->SetLineColor(kRed);
  hyptar0->SetLineWidth(2);
 // hyptar0->Scale(1./hyptar0->GetMaximum());
  hyptar0->Draw("histo");


TH1F *hyptar2 = (TH1F*)f->Get("h_yptar_2");
  hyptar2->SetName("hyptar2");
  hyptar2->SetLineColor(kGreen);
  hyptar2->SetLineWidth(2);
  //hyptar2->Scale(1./hyptar2->GetMaximum());
  hyptar2->Draw("histosame");

TH1F *hyptar1 = (TH1F*)f->Get("h_yptar_1");
  hyptar1->SetName("hyptar1");
  hyptar1->SetLineColor(kBlue);
  hyptar1->SetLineWidth(2);
 // hyptar1->Scale(1./hyptar1->GetMaximum());
  hyptar1->Draw("histosame");

//leg->Draw();
canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+"(").c_str());

// 2D Histos Overlays
TH2F *h20 = (TH2F*)f->Get("h2_dptheta_0");
h20->SetName("h20");
h20->SetTitle("Theta vs Delta");
h20->SetMarkerColor(kRed);
//h20->Scale(1./h20->GetMaximum());
h20->Draw("histo");

TH2F *h21 = (TH2F*)f->Get("h2_dptheta_1");
h21->SetName("h21");
h21->SetMarkerColor(kBlue);
//h21->Scale(1./h21->GetMaximum());
h21->Draw("histosame");

TH2F *h22 = (TH2F*)f->Get("h2_dptheta_2");
h22->SetName("h22");
h22->SetMarkerColor(kGreen);
//h22->Scale(1./h22->GetMaximum());
h22->Draw("histosame");

canvas->Print((pdf_file_name + kinSetting +  MatrixEl + file_format+"(").c_str());


TH2F *h2W2yp0 = (TH2F*)f->Get("h2_kinW2theta_0");
h2W2yp0->SetName("h2W2yp0");
h2W2yp0->SetTitle("Theta vs KinW2");
h2W2yp0->SetMarkerColor(kRed);
h2W2yp0->Draw("histo");

TH2F *h2W2yp1 = (TH2F*)f->Get("h2_kinW2theta_1");
h2W2yp1->SetName("h2W2yp1");
h2W2yp1->SetMarkerColor(kBlue);
h2W2yp1->Draw("histosame");

TH2F *h2W2yp2 = (TH2F*)f->Get("h2_kinW2theta_2");
h2W2yp2->SetName("h2W2yp2");
h2W2yp2->SetMarkerColor(kGreen);
h2W2yp2->Draw("histosame");

canvas->Print((pdf_file_name + kinSetting + MatrixEl + file_format+")").c_str());

fout->Write();
fout->Close();


}
