
void plotGlamour(int nrun, int FileID, int foilPOS){
// use foilPOS = 8 for 7.5 cm tgt

std::string MatrixEl=/*"angular6_67";*/"coeff2018";
std::string file_format=".pdf";
std::string pdf_fname=Form("Output/carbonSieve/plots/GlamourPlots_");

TFile * f = TFile::Open(Form("./Output/carbonSieve/hist/analysis_SieveData_coeff2018_%i_%ifoil_%i.root",nrun,foilPOS,FileID));

TFile *fout = new TFile(Form("Output/carbonSieve/hist/GlamourPlots_coeff2018_%i_%ifoil_%i.root",nrun, foilPOS, FileID),"RECREATE");

//make the output file

auto *canvas = new TCanvas("canvas", "HMS Sieve",50,50,700,1200);
gROOT->ForceStyle();
gROOT->SetBatch(true);
gStyle->SetOptStat(0);
canvas->SetFillColor(0);
canvas->SetBorderMode(0);
canvas->SetBorderSize(0);
canvas->SetFrameFillColor(0);
canvas->SetFrameBorderMode(0);
canvas->SetGridx();
canvas->SetGridy();

//save plots
canvas->Update();

//2D Histos
TH2F *h1 = (TH2F*)f->Get("h2_sieve");
  h1->SetName("h1");
  h1->SetTitle("HMS Sieve");
  h1->Draw("histo"); 

canvas->Print((pdf_fname + MatrixEl + file_format+")").c_str());

fout->Write();
fout->Close();
}
