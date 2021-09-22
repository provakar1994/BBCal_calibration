#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

static const Int_t kNrows=26;
static const Int_t kNcols=2;
static const Double_t target_FADC_amp=10.; //mV

void targetFADCamp(){

  vector<double> trigtoFADC;
  vector<int> elemID;
  
  string InFile = "trigtoFADCcoef_PS.txt";
  ifstream infile_data;
  infile_data.open(InFile);
  string readline;
  if (infile_data.is_open() ) {
    while (getline(infile_data,readline)) {
      double ratio;
      int elem;
      infile_data >> elem >> ratio;
      elemID.push_back( elem );
      trigtoFADC.push_back( ratio );
    }
  }
  
  TH2F* target_amp = new TH2F("target_amp"," Target FADC Amplitude ; Ncol ; Nrow",kNcols,1,kNcols+1,kNrows,1,kNrows+1);

  string OutFile = "Output/target_FADC_amp_PS.txt";
  cout << " Write target FADC amplitudes to : " << OutFile << endl;
  ofstream outfile_data;
  outfile_data.open(OutFile);
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      Double_t target_FADC_amp_pCh = target_FADC_amp/trigtoFADC.at(r*kNcols+c); 
      target_amp->Fill(float(c+1),float(r+1),target_FADC_amp_pCh);
      outfile_data << elemID.at(r*kNcols+c) << " " << target_FADC_amp_pCh << endl;
    }
  }
  outfile_data.close();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetPaintTextFormat("4.2f");  
  const Int_t Number=3;
  Double_t Red[Number] = { 1.0,0.0,0.0};
  Double_t Blue[Number] = { 1.0,0.0,1.0};
  Double_t Green[Number] = { 0.0,1.0,0.0};
  Double_t Len[Number] = { 0.0,.5,1.0};
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Len,Red,Green,Blue,nb);
  TCanvas* can2d;
  can2d= new TCanvas("can_2d","2d ",700,1000);
  can2d->cd();
  target_amp->SetMaximum(16);
  target_amp->SetMinimum(8);
  target_amp->Draw("text colz");

  can2d->SaveAs("Output/targetFADCamp_PS.pdf");
  //can2d->Close();
    
}
