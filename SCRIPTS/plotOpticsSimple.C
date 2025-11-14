void plotOpticsSimple(int nrun, int nEvent, int foilPOS){
  //use nEvent = -1 for full replay
  //round foil positions
  //read the input file

//to do: make a script to loop over delta cuts and fit H.react.z to 1 sigma
//
  //TFile *f = new TFile(Form("../ROOTfiles/OPTICS/6_667GeV/nps_hms_optics_hadd_ang6_667fit_%i_1_%i.root",nrun, nEvent));
 TFile *f = new TFile(Form("../ROOTfiles/OPTICS/nps_hms_optics_pass1Test_%i_1_%i.root",nrun, nEvent));
 //TFile *f = new TFile(Form("../ROOTfiles/OPTICS/5_878GeV/nps_hms_optics_%i_1_%i.root",nrun, nEvent));

 TTree *tt = (TTree*)f->Get("T");


  //foil cut
  int foilmin = foilPOS-2.5;
  int foilmax = foilPOS+2.5;
  cout << foilmin << " is lower foilcut" << endl;
  cout << foilmax << " is upper foilcut" << endl;

  //here's the cut
  TCut cut = Form("H.react.z>%i&&H.react.z<%i&&abs(H.gtr.dp)<8&&H.cal.etracknorm>0.65&&H.cer.npeSum>6",foilmin,foilmax);
  TCut cutCentral = "H.extcor.xsieve<1&&H.extcor.xsieve>-2&&abs(H.extcor.ysieve)<0.6";
  TCut edgeCut = "H.extcor.xsieve<1&&H.extcor.xsieve>-2&&H.extcor.ysieve<-4.8&&H.extcor.ysieve>-5.8";

  //make the output file
  TCanvas *canvas = new TCanvas("canvas","plots",800,800);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  std::string pdf_file_name= Form("Output/carbonSieve/plots/analysis_pass1Test_%i_%ifoil_%i.pdf",nrun,foilPOS,nEvent);
//std::string pdf_file_name= Form("Output/carbonSieve/plots/analysis_SieveSimple_%i_%ifoil_%i.pdf",nrun,foilPOS,nEvent);

//  std::string pdf_file_name= Form("Output/carbonSieve/plots/analysis_SieveSimple_6p59_%i_%ifoil_%i.pdf",nrun,foilPOS,nEvent);

  gROOT->SetBatch(true);
  gROOT->ForceStyle();
  gStyle->SetOptStat(1);
  canvas->SetGridx();
  canvas->SetGridy();

//TFile *fout = new TFile(Form("Output/carbonSieve/hist/analysis_SieveSimple_%i_%ifoil_%i.root",nrun,foilPOS,nEvent),"RECREATE");

  TFile *fout = new TFile(Form("Output/carbonSieve/hist/analysis_pass1Test_%i_%ifoil_%i.root",nrun,foilPOS,nEvent),"RECREATE");
  //TFile *fout = new TFile(Form("Output/carbonSieve/hist/analysis_SieveSimple_6p59_%i_%ifoil_%i.root",nrun,foilPOS,nEvent),"RECREATE");


  //make the plots
  TH1F *h_z = new TH1F("h_z",";H.react.z [cm]",100,-14,14);
  TH1F *h_xsieve = new TH1F("h_xsieve",";xSieve[cm]",200,-12.0,12.0); 
  TH1F *h_ysieve = new TH1F("h_ysieve",";ySieve[cm]",200,-7.0,7.0);
  TH2F *h2_ypVz = new TH2F("h2_ypVz",";zVertex [cm];ypTar",100,-14,14,100,-0.05,0.05);
  TH2F *h2_ypVy = new TH2F("h2_ypVy",";yTar [cm];ypTar",100,-5,5,100,-0.05,0.05);
  TH2F *h2_yfpVxfp = new TH2F("h2_yfpVxfp",";xfp [cm];yfp [cm]",100,0,8,100,-10,10);
  TH2F *h2_dpVz = new TH2F("h2_dpVz",";zVertex [cm];delta",100,-14,14,100,-12,12);
  TH2F *h2_sieve = new TH2F("h2_sieve",";ySieve [cm];xSieve[cm]",200,-7.0,7.0,200,-12.0,12.0);
  TH2F *h2_xpVd = new TH2F("h2_xpVd",";delta;xpfp",100,-12,12,100,-0.15,0.15);
 // TH2F *h2_kinW = new TH2F("h2_kinW",";H.kin.W;H.gtr.dp",100,0,5,100,-12,12);  

  //plot this stuff format Yaxis:Xaxis
  tt->Draw("H.react.z>>h_z",cut);
  tt->Draw("H.extcor.xsieve>>h_xsieve",cut);
  tt->Draw("H.extcor.ysieve>>h_ysieve",cut);
  //looks like ypTar is filled by H.gtr.ph
  tt->Draw("H.gtr.ph:H.react.z>>h2_ypVz",cut);
  tt->Draw("H.gtr.dp:H.react.z>>h2_dpVz",cut);
  tt->Draw("H.gtr.ph:H.gtr.y>>h2_ypVy",cut);
 //looks like yfp filled by H.dc.y_fp
  tt->Draw("H.dc.y_fp:H.dc.x_fp>>h2_yfpVxfp",cut);
  tt->Draw("H.extcor.xsieve:H.extcor.ysieve>>h2_sieve",cut);
  tt->Draw("H.dc.xp_fp:H.gtr.dp>>h2_xpVd",cut);
  
  //save plots
  canvas->Update();

  h_z->SetLineWidth(2);
  h_z->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h_xsieve->SetLineWidth(2);
  h_xsieve->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h_ysieve->SetLineWidth(2);
  h_ysieve->Draw();
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVz->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_dpVz->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_ypVy->Draw("colz");
  canvas->Print((pdf_file_name +"(").c_str());

  h2_yfpVxfp->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());
  
  h2_sieve->Draw("colz");  
  canvas->Print((pdf_file_name +"(").c_str());

  h2_xpVd->Draw("colz");
  canvas->Print((pdf_file_name +")").c_str());  

  fout->Write();
  fout->Close();

}
