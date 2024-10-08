#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>



void fit_opt_matrix_v2(Int_t nSettings = 1) {
  //  Int_t nSettings = 1; //number of files
  Int_t FileID=-1;
  Int_t maxFoils=2;
  Int_t maxDel=10;

  vector<int> runTot;
  runTot.push_back(1264);
  runTot.push_back(1263);

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetPadLeftMargin(0.17);
  //
  string newcoeffsfilename="nps_hms_newfit_5pt6GeV_all.dat";
  string oldcoeffsfilename="DATfiles/hms_recon_coeff_opt2018.dat";

  int nfit=0,npar,nfit_max=45000,npar_final=0,max_order=6,norder;
  Int_t MaxPerBin=1000;
  Int_t MaxZtarPerBin=15000;

  //
  TH1F *hDelta = new TH1F("hDelta","Delta ",20,-10.,30.);
  TH1F *hDeltarecon = new TH1F("hDeltarecon","Delta Recon % ",40,-10,30);
  TH1F *hDeltadiff = new TH1F("hDeltadiff","Delta Diff % ",40,-1,1);
  TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon % ",40,-10,30);
  TH1F *hDeltanewdiff = new TH1F("hDeltanewdiff","Delta New Diff % ",40,-1,1);
  TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
  TH1F *hytarrecon = new TH1F("hytarrecon","ytar recon(cm)",70,-7.,7.);
  TH1F *hytarnew = new TH1F("hytarnew","ytar new(cm)",70,-7.,7.);
  TH1F *hytardiff = new TH1F("hytardiff","ytar diff(cm)",70,-5.,5.);
  TH1F *hytarnewdiff = new TH1F("hytarnewdiff","ytar new diff(cm)",70,-5.,5.);
  TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
  TH1F *hxptarrecon = new TH1F("hxptarrecon","xptar recon",100,-.12,.12);
  TH1F *hxptardiff = new TH1F("hxptardiff","xptar diff (mr)",100,-15,15);
  TH1F *hxptarnew = new TH1F("hxptarnew","xptar new recon",100,-.12,.12);
  TH1F *hxptarnewdiff = new TH1F("hxptarnewdiff","xptar new diff (mr)",100,-15,15);
  TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
  TH1F *hyptarrecon = new TH1F("hyptarrecon","yptar recon ",100,-.1,.1);
  TH1F *hyptardiff = new TH1F("hyptardiff","yptar diff (mr) ",100,-10,10);
  TH1F *hyptarnew = new TH1F("hyptarnew","yptar new recon ",100,-.1,.1);
  TH1F *hyptarnewdiff = new TH1F("hyptarnewdiff","yptar new diff (mr) ",100,-10,10);
  //
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
 ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
   if(!oldcoeffsfile.is_open()) {
     cout << " error opening reconstruction coefficient file: " << oldcoeffsfilename.c_str() << endl;
     return;
   } else {
     cout << "Open  Old coeff file = " << oldcoeffsfilename.c_str() << endl;
   }
   string line="!";
  int good = getline(oldcoeffsfile,line).good();
  int nskip=0;
  if ( line[0] =='!') {
  while(good && line[0]=='!') {
    good = getline(oldcoeffsfile,line).good();
    nskip++;
    //    cout << line << endl;
  }
  }
  cout << " skipped " << nskip << " comments lines in file : " << oldcoeffsfilename.c_str() << endl;
  cout << " at line = " << line.c_str() << endl;
  nskip=0;
 if ( line.compare(0,4," ---")==0) {
  while(good && line.compare(0,4," ---")==0) {
    good = getline(oldcoeffsfile,line).good();
    nskip++;
  }
 }
  cout << " skipped " << nskip << " separation lines in file : " << oldcoeffsfilename.c_str() << endl;
  cout << line.c_str() << endl;
 
  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_fit;
  vector<double> yptarcoeffs_fit;
  vector<double> ytarcoeffs_fit;
  vector<double> deltacoeffs_fit;
  vector<int> xfpexpon_fit;
  vector<int> xpfpexpon_fit;
  vector<int> yfpexpon_fit;
  vector<int> ypfpexpon_fit;
  vector<int> xtarexpon_fit;

  vector<double> xptarcoeffs_xtar;
  vector<double> yptarcoeffs_xtar;
  vector<double> ytarcoeffs_xtar;
  vector<double> deltacoeffs_xtar;
  vector<int> xfpexpon_xtar;
  vector<int> xpfpexpon_xtar;
  vector<int> yfpexpon_xtar;
  vector<int> ypfpexpon_xtar;
  vector<int> xtarexpon_xtar;

  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue,angle,ztartrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  TString currentline;
  int num_recon_terms_old;
  int num_recon_terms_fit;
  int num_recon_terms_xtar;

  num_recon_terms_old = 0;
  num_recon_terms_fit = 0;
  num_recon_terms_xtar = 0;
  // add zero order term to fit
  xptarcoeffs_fit.push_back(0.0);
  ytarcoeffs_fit.push_back(0.0);
  yptarcoeffs_fit.push_back(0.0);
  deltacoeffs_fit.push_back(0.0);
  xfpexpon_fit.push_back(0);
  xpfpexpon_fit.push_back(0);
  yfpexpon_fit.push_back(0);
  ypfpexpon_fit.push_back(0);
  xtarexpon_fit.push_back(0);
  num_recon_terms_fit = 1;
    //
    Double_t Coeff[4];
    Int_t Exp[5];
  while( good && line.compare(0,4," ---")!=0 ){
    sscanf(line.c_str()," %le %le %le %le %1d%1d%1d%1d%1d"
	   ,&Coeff[0],&Coeff[1]
	   ,&Coeff[2],&Coeff[3]
	   ,&Exp[0]
	   ,&Exp[1]
	   ,&Exp[2]
	   ,&Exp[3]
	   ,&Exp[4]);
    
    xptarcoeffs_old.push_back(Coeff[0]);
    ytarcoeffs_old.push_back(Coeff[1]);
    yptarcoeffs_old.push_back(Coeff[2]);
    deltacoeffs_old.push_back(Coeff[3]);
    

    xfpexpon_old.push_back(Exp[0]);
    xpfpexpon_old.push_back(Exp[1]);
    yfpexpon_old.push_back(Exp[2]);
    ypfpexpon_old.push_back(Exp[3]);
    xtarexpon_old.push_back(Exp[4]);

    num_recon_terms_old++;
    norder= Exp[0]+Exp[1]+Exp[2]+Exp[3]+Exp[4];
    if (Exp[4]==0) {
    xptarcoeffs_fit.push_back(Coeff[0]);
    ytarcoeffs_fit.push_back(Coeff[1]);
    yptarcoeffs_fit.push_back(Coeff[2]);
    deltacoeffs_fit.push_back(Coeff[3]);
    

    xfpexpon_fit.push_back(Exp[0]);
    xpfpexpon_fit.push_back(Exp[1]);
    yfpexpon_fit.push_back(Exp[2]);
    ypfpexpon_fit.push_back(Exp[3]);
    xtarexpon_fit.push_back(Exp[4]);
    num_recon_terms_fit++;
    } else {
    xptarcoeffs_xtar.push_back(Coeff[0]);
    ytarcoeffs_xtar.push_back(Coeff[1]);
    yptarcoeffs_xtar.push_back(Coeff[2]);
    deltacoeffs_xtar.push_back(Coeff[3]);
    

    xfpexpon_xtar.push_back(Exp[0]);
    xpfpexpon_xtar.push_back(Exp[1]);
    yfpexpon_xtar.push_back(Exp[2]);
    ypfpexpon_xtar.push_back(Exp[3]);
    xtarexpon_xtar.push_back(Exp[4]);
    num_recon_terms_xtar++;
    }
    good = getline(oldcoeffsfile,line).good();
  }

  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  cout << "num recon terms in fit matrix = " << num_recon_terms_fit << endl;
  cout << "num recon terms in xtar matrix = " << num_recon_terms_xtar << endl;
  npar= num_recon_terms_fit ;
  //
  TVectorD b_ytar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
  //
  
  
  const Int_t nysieve=9;  
  vector <Double_t> ys_cent;
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=(nys-4)*0.6*2.54;
    ys_cent.push_back(pos);
  }
  
  for (int iSetting=0; iSetting<nSettings; iSetting++){
    //  Get info for that optics run
    Int_t nrun = runTot[iSetting];
    TString OpticsFile = "DATfiles/list_of_optics_run.dat";
    ifstream file_optics(OpticsFile.Data());
    TString opticsline;
    TString OpticsID="";
    Int_t RunNum=0.;
    Double_t CentAngle=0.;
    Int_t SieveFlag=1;
    Double_t ymis =0.0;
    Int_t nfoils=0;
    TString temp;
    //
    vector <Double_t> ztar_foil;
    Int_t ndelcut=-1;
    vector<Double_t > delcut;
   // vector<Double_t > delwidth;
    if (file_optics.is_open()) {
      //
      cout << " Open file = " << OpticsFile << endl;
      while (RunNum!=nrun  ) {
	temp.ReadToDelim(file_optics,',');
//	cout << temp << endl;
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
	nfoils = temp.Atoi();
	temp.ReadToDelim(file_optics,',');
	SieveFlag = temp.Atoi();
	temp.ReadToDelim(file_optics);
	ndelcut = temp.Atoi();
//	temp.ReadToDelim(file_optics);
  //      ymis = temp.Atof();
	for (Int_t nf=0;nf<nfoils-1;nf++) {
	  temp.ReadToDelim(file_optics,',');
	  ztar_foil.push_back(temp.Atof());
	  cout << "ztar foil " << ztar_foil[nf] << endl;
	}
	temp.ReadToDelim(file_optics);
	ztar_foil.push_back(temp.Atof());
	for (Int_t nd=0;nd<ndelcut-1;nd++) {
	  temp.ReadToDelim(file_optics,',');
	  delcut.push_back(temp.Atof());
	  cout << " nd = " << nd << " " << delcut[nd] << endl;
	}
	temp.ReadToDelim(file_optics);
	delcut.push_back(temp.Atof());
}
//what is this below
/*	for (Int_t nw=0;nw<ndelcut-1;nw++) {
	  temp.ReadToDelim(file_optics,',');
	  delwidth.push_back(temp.Atof());
	}
	temp.ReadToDelim(file_optics);
	delwidth.push_back(temp.Atof());
      }
*/

    } else {
      cout << " No file = " << OpticsFile << endl;    
    }
    cout << RunNum << " " << OpticsID << " " << CentAngle << " " << nfoils << " " << SieveFlag << endl;
    
    TString inputroot;
    inputroot = Form("hist/Optics_%s_%d_fit_tree.root",OpticsID.Data(),FileID);
    cout << " INfile = " << inputroot << endl;
    TFile *fsimc = new TFile(inputroot);
    TTree *FitTree = (TTree*)fsimc->Get("TFit");
    //Declaration of leaves types
    Double_t  ys,xtar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT,xtarT;
    Double_t xfp,xpfp,yfp,ypfp,ysieveT,ysieve;
    FitTree->SetBranchAddress("ys",&ysieve);
    FitTree->SetBranchAddress("ysT",&ysieveT);
    FitTree->SetBranchAddress("xtar",&xtar);
    FitTree->SetBranchAddress("xtarT",&xtarT);
    FitTree->SetBranchAddress("xptar",&xptar);
    FitTree->SetBranchAddress("yptar",&yptar);
    FitTree->SetBranchAddress("ytar",&ytar);
    FitTree->SetBranchAddress("xptarT",&xptarT);
    FitTree->SetBranchAddress("yptarT",&yptarT);
    FitTree->SetBranchAddress("ytarT",&ytarT);
    FitTree->SetBranchAddress("ztarT",&ztarT);
    FitTree->SetBranchAddress("delta",&delta);
    FitTree->SetBranchAddress("xpfp",&xpfp);
    FitTree->SetBranchAddress("ypfp",&ypfp);
    FitTree->SetBranchAddress("xfp",&xfp);
    FitTree->SetBranchAddress("yfp",&yfp);
    //

    vector<Int_t > Ztar_Cnts;
    vector<vector<vector<Int_t> > > Ztar_Ys_Delta_Cnts;
    vector<Int_t> Max_Per_Run_Per_Foil;
    Ztar_Cnts.resize(maxFoils);
    Ztar_Ys_Delta_Cnts.resize(maxFoils);
    Max_Per_Run_Per_Foil.resize(maxFoils);
   
    for (Int_t nf=0;nf<maxFoils;nf++) {//max foils
      Ztar_Ys_Delta_Cnts[nf].resize(maxDel);//max del cut
      for (Int_t nd=0;nd<maxDel;nd++) {
	Ztar_Ys_Delta_Cnts[nf][nd].resize(nysieve);
      }
    }
    //
    for (Int_t nf=0;nf<maxFoils;nf++) {//max foils
      Max_Per_Run_Per_Foil[nf] = 30000;//check this number?
      Ztar_Cnts[nf]=0;
      for (Int_t nd=0;nd<maxDel;nd++) {//max del cut
	for (Int_t ny=0;ny<nysieve;ny++) {	
	  Ztar_Ys_Delta_Cnts[nf][nd][ny]=0;
	}}}
    
  Long64_t nentries = FitTree->GetEntries();
  for (int i = 0; i < nentries; i++) {
    FitTree->GetEntry(i);
    //
    if (TMath::Abs(delta) < 100. ) {
      Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
      Double_t etemp;
      for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar/100., xtarexpon_old[icoeffold] );
	deltatemp += deltacoeffs_old[icoeffold] * etemp;
	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
      } // for icoeffold loop
      hytarrecon->Fill(ytartemp*100.);
      hyptarrecon->Fill(yptartemp);
      hxptarrecon->Fill(xptartemp);
      if ( delta>-15. &&  delta<30. ) {
	//
	Int_t found_nf=-1;
	Int_t found_nd=-1;
	Int_t found_ny=-1;
	Bool_t good_bin=kFALSE;
	for (Int_t nf=0;nf<nfoils;nf++) {
	  if (abs(ztarT-ztar_foil[nf])<2.5) found_nf=nf;
	}
	for (Int_t nd=0;nd<ndelcut;nd++) {
	  if (delta >=delcut[nd] && delta <delcut[nd+1]) found_nd=nd;
	}
	for (Int_t ny=0;ny<nysieve;ny++) {	
	  if (abs(ysieveT-ys_cent[ny])<.5) found_ny=ny;
	}
	if (found_nf!=-1 &&found_nd!=-1 && found_ny!=-1)  {
	  good_bin=kTRUE;
	}

	if (good_bin && nfit < nfit_max && Ztar_Ys_Delta_Cnts[found_nf][found_nd][found_ny]< MaxPerBin && Ztar_Cnts[found_nf]< MaxZtarPerBin && Ztar_Cnts[found_nf]<Max_Per_Run_Per_Foil[found_nf]) {

	  //reconstruct it
	  Double_t ytar_xtar = 0.0,yptar_xtar=0.0,xptar_xtar=0.0;
          for( int icoeff_xtar=0; icoeff_xtar<num_recon_terms_xtar; icoeff_xtar++ ){
	    etemp= 
	      pow( xfp / 100.0, xfpexpon_xtar[icoeff_xtar] ) * 
	      pow( yfp / 100.0, yfpexpon_xtar[icoeff_xtar] ) * 
	      pow( xpfp, xpfpexpon_xtar[icoeff_xtar] ) * 
	      pow( ypfp, ypfpexpon_xtar[icoeff_xtar] ) * 
	      pow( xtar/100., xtarexpon_xtar[icoeff_xtar] );
	    ytar_xtar += ytarcoeffs_xtar[icoeff_xtar] * etemp;
	    yptar_xtar += yptarcoeffs_xtar[icoeff_xtar] * etemp;
	    xptar_xtar += xptarcoeffs_xtar[icoeff_xtar] *etemp; 
	  }
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
	    etemp= 
	      pow( xfp / 100.0, xfpexpon_fit[icoeff_fit] ) * 
	      pow( yfp / 100.0, yfpexpon_fit[icoeff_fit] ) * 
	      pow( xpfp, xpfpexpon_fit[icoeff_fit] ) * 
	      pow( ypfp, ypfpexpon_fit[icoeff_fit] ) * 
	      pow( xtar/100., xtarexpon_fit[icoeff_fit] );
	    if (nfit < nfit_max ) {
              lambda[icoeff_fit][nfit] = etemp;
	      b_xptar[icoeff_fit] += (xptarT-xptar_xtar) * etemp;
	      b_yptar[icoeff_fit] += (yptarT-yptar_xtar) * etemp;
	      b_ytar[icoeff_fit] += (ytarT-ytar_xtar*100) /100.0 * etemp;
	    }
	  } // for icoeff_fit loop
	   
	  ////////////////////////////////////////////////////////////

	  hytar->Fill(ytar);
	  hyptar->Fill(yptar);
	  hxptar->Fill(xptar);
	  hytardiff->Fill(ytar-ytarT);
	  hyptardiff->Fill(1000.*(yptar-yptarT));
	  hxptardiff->Fill(1000.*(xptar-xptarT));
	  Ztar_Cnts[found_nf]++;
	  Ztar_Ys_Delta_Cnts[found_nf][found_nd][found_ny]++;
	  nfit++;
	  xfptrue.push_back( xfp );
	  yfptrue.push_back( yfp );
	  xpfptrue.push_back( xpfp );
	  ypfptrue.push_back( ypfp );
	  xtartrue.push_back( xtarT );//used to be xtar
	  xptartrue.push_back( xptarT );
	  ytartrue.push_back( ytarT );
	  yptartrue.push_back( yptarT );
	  ztartrue.push_back(ztarT);
	  angle.push_back(CentAngle*(3.14159)/180.0);
	}
      }
    }
  }
  //
  for (Int_t nf=0;nf<maxFoils;nf++) cout << " counts foil " << nf << " : " << Ztar_Cnts[nf] << endl;
   
 for (Int_t nf=0;nf<nfoils;nf++) {
    cout << " ztar = " << ztar_foil[nf] << endl;
    for (Int_t nd=0;nd<ndelcut-1;nd++) {
      cout << " Ndelta = " << (delcut[nd]+delcut[nd+1])/2 << endl;       
      for (Int_t ny=0;ny<nysieve;ny++) {
	cout <<  Ztar_Ys_Delta_Cnts[nf][nd][ny] << " " ;
      }
      cout << endl;
    }}
  
  //
   //
  ////////////////////
  //end each run loop
  ////////////////////
  }
   //
  if (nfit < nfit_max) {
    cout << " nfit < nfit_max, set nfit_max = " << nfit << endl;
    return;
  }
  //
  cout << " number to fit = " << nfit << " max = " << nfit_max << endl;
  for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
  }
  for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 5000 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
      }
    }
  }
 
  TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_ytar );
  cout << "ytar solution ok = " << ok << endl;
  //b_ytar.Print();
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  //b_yptar.Print();
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  //b_xptar.Print();


  ///////////////////////////////////////////////////////////////////////////////////////////////
    for( int ifit=0; ifit<nfit; ifit++){
          Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
	  Double_t etemp;
     for( int ipar=0; ipar<npar; ipar++){
       etemp=lambda[ipar][ifit];
        	ytarnew += b_ytar[ipar] * etemp;
	        yptarnew += b_yptar[ipar] * etemp;
	         xptarnew += b_xptar[ipar] *etemp;        
    }
          Double_t ytar_xtar = 0.0,yptar_xtar=0.0,xptar_xtar=0.0;
          for( int icoeff_xtar=0; icoeff_xtar<num_recon_terms_xtar; icoeff_xtar++ ){
        	etemp= 
	  pow( xfptrue.at(ifit) / 100.0, xfpexpon_xtar[icoeff_xtar] ) * 
	  pow( yfptrue.at(ifit) / 100.0, yfpexpon_xtar[icoeff_xtar] ) * 
	  pow( xpfptrue.at(ifit), xpfpexpon_xtar[icoeff_xtar] ) * 
	  pow( ypfptrue.at(ifit), ypfpexpon_xtar[icoeff_xtar] ) * 
	  pow( xtartrue.at(ifit)/100., xtarexpon_xtar[icoeff_xtar] );
        	ytar_xtar += ytarcoeffs_xtar[icoeff_xtar] * etemp;
	        yptar_xtar += yptarcoeffs_xtar[icoeff_xtar] * etemp;
		       xptar_xtar += xptarcoeffs_xtar[icoeff_xtar] *etemp; 
	  }
	  hytarnew->Fill((ytarnew+ytar_xtar)*100.);
	  hyptarnew->Fill(yptarnew+yptar_xtar);
	  hxptarnew->Fill(xptarnew);
	  hytarnewdiff->Fill((ytarnew+ytar_xtar)*100.-ytartrue.at(ifit));
	  hyptarnewdiff->Fill(1000*(yptarnew+yptar_xtar-yptartrue.at(ifit)));
	  hxptarnewdiff->Fill(1000*(xptarnew+xptar_xtar-xptartrue.at(ifit)));
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // write out coeff
  char coeffstring[100];
  Double_t tt;
  cout << "writing new coeffs file" << endl;
  //newcoeffsfile << "! new fit to "<<endl;//+fname << endl;
  //newcoeffsfile << " ---------------------------------------------" << endl;
  for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
    newcoeffsfile << " ";
    //      tt=xptarcoeffs_fit[icoeff_fit];
    tt=b_xptar[icoeff_fit] ;
    sprintf( coeffstring, "%16.9g", tt );
    newcoeffsfile << coeffstring; 
    //      newcoeffsfile << " ";
    sprintf( coeffstring, "%16.9g", b_ytar[icoeff_fit] );
    newcoeffsfile << coeffstring;
    sprintf( coeffstring, "%16.9g", b_yptar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    sprintf( coeffstring, "%16.9g", deltacoeffs_fit[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    newcoeffsfile << " ";
    newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_fit[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_fit[icoeff_fit]; 
    newcoeffsfile << endl;
    
  }
  //
  for( int icoeff_fit=0; icoeff_fit<num_recon_terms_xtar; icoeff_fit++ ){
    newcoeffsfile << " ";
    //      tt=xptarcoeffs_fit[icoeff_fit];
    tt=xptarcoeffs_xtar[icoeff_fit] ;
    sprintf( coeffstring, "%16.9g", tt );
    newcoeffsfile << coeffstring; 
    //      newcoeffsfile << " ";
    sprintf( coeffstring, "%16.9g", ytarcoeffs_xtar[icoeff_fit] );
    newcoeffsfile << coeffstring;
    sprintf( coeffstring, "%16.9g", yptarcoeffs_xtar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    sprintf( coeffstring, "%16.9g", deltacoeffs_xtar[icoeff_fit] );
    //newcoeffsfile << " ";
    newcoeffsfile << coeffstring; 
    newcoeffsfile << " ";
    newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_xtar[icoeff_fit]; 
    newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_xtar[icoeff_fit]; 
    newcoeffsfile << endl;
    
  }
  //
  newcoeffsfile << " ---------------------------------------------" << endl;
  
  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
  //
  TCanvas *cdiff = new TCanvas("cdiff","Old matrix Diff target",800,800);
  cdiff->Divide(2,2);
  cdiff->cd(1);
  hytardiff->Draw();
  hytardiff->Fit("gaus");
  TF1 *fitcydiff=hytardiff->GetFunction("gaus");
  cdiff->cd(2);
  hyptardiff->Draw();
  hyptardiff->Fit("gaus");
  TF1 *fitcypdiff=hyptardiff->GetFunction("gaus");
  cdiff->cd(3);
  hxptardiff->Draw();
  hxptardiff->Fit("gaus");
  TF1 *fitcxpdiff=hxptardiff->GetFunction("gaus");
  cdiff->cd(4);
  //  hDeltadiff->Draw();
  //hDeltadiff->Fit("gaus");
  //TF1 *fitcdeldiff=hDeltadiff->GetFunction("gaus");
  //
  TCanvas *cnewdiff = new TCanvas("cnewdiff","Newfit diff target",800,800);
  cnewdiff->Divide(2,2);
  cnewdiff->cd(1);
  hytarnewdiff->Draw();
  hytarnewdiff->Fit("gaus");
  TF1 *fitcynewdiff=hytarnewdiff->GetFunction("gaus");
  cnewdiff->cd(2);
  hyptarnewdiff->Draw();
  hyptarnewdiff->Fit("gaus");
  TF1 *fitcypnewdiff=hyptarnewdiff->GetFunction("gaus");
  cnewdiff->cd(3);
  hxptarnewdiff->Draw();
  hxptarnewdiff->Fit("gaus");
  TF1 *fitcxpnewdiff=hxptarnewdiff->GetFunction("gaus");
  cnewdiff->cd(4);
  //   hDeltanewdiff->Draw();
  //hDeltanewdiff->Fit("gaus");
  //TF1 *fitcdelnewdiff=hDeltanewdiff->GetFunction("gaus");
  //
 
}
