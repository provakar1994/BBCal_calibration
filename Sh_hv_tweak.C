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
#include <TObjString.h>
#include <TLatex.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>

using namespace std;

static const Int_t shNCol=7;
static const Int_t shNRow=27;

Double_t HV_Value[2][16][12];
Int_t HV_Block[2][16][12];

void SetHVMap();
void ReadHV(Int_t nrun);
void TweakHV();
void ReadAlpha();
void ReadGain(Int_t nrun);

vector<double> gAlpha;
vector<double> adcGain; // Gain = New/Old

Double_t HVUpdate[shNCol*shNRow];
Double_t HV_Crate[shNCol*shNRow];
Double_t HV_Slot[shNCol*shNRow];
Double_t HV_Chan[shNCol*shNRow];
//

void SetHVMap() {
  cout << " Setting HV map " << endl;
  for (Int_t nc=0;nc<2;nc++) {
    for (Int_t ns=0;ns<16;ns++) {
      for (Int_t nch=0;nch<12;nch++) {
	HV_Block[nc][ns][nch]=-1;
      }
    }
  } 
  Int_t nslot=2;// HV slots = 0 to 15
  Int_t nchan=3;// HV channels = 0 to 11
  for (Int_t nr=0;nr<shNRow;nr++) {
    if (nr==12) nslot=5;
    if (nr==12) nchan=0;
    for (Int_t nc=0;nc<shNCol;nc++) {
      if (nr > 11) {
	HV_Block[1][nslot][nchan++] = nr*shNCol+nc;
 	HV_Slot[nr*shNCol+nc] =nslot ;
 	HV_Crate[nr*shNCol+nc] =1 ;
 	HV_Chan[nr*shNCol+nc] =nchan-1 ;
      } else {
	HV_Block[0][nslot][nchan++] = nr*shNCol+nc;
  	HV_Slot[nr*shNCol+nc] =nslot ;
 	HV_Crate[nr*shNCol+nc] =0 ;
 	HV_Chan[nr*shNCol+nc] =nchan-1 ;
      }
      if (nchan==12) {
	nchan=0;
	nslot++;
      }
    }
  }
}//

void ReadHV(Int_t nrun) {
  TString HVfileName = Form("hv_set/run_%d_hv.set",nrun);
  cout << " read file = " << HVfileName << endl;
  ifstream file_hv(HVfileName.Data());
  string hvline;
  string crate;
  string slot;
  Int_t CrateNum;
  Int_t SlotNum;
  if (file_hv.is_open()) {
    while (getline(file_hv,hvline)) {
      if (hvline.compare(0,1,"#")!=0) {	     
	istringstream tokenStream(hvline);
	string token;
	char delimiter = ' ';
	vector<TString> tokens;
	while (getline(tokenStream,token,delimiter))  {
	  tokens.push_back(token);
	}
	CrateNum=-1;
	if (tokens[0] == "rpi17:2001") CrateNum=0;
	if (tokens[0] == "rpi18:2001") CrateNum=1;
	TString SlotNumStr(tokens[1](1,tokens[1].Sizeof()+1));
	SlotNum=SlotNumStr.Atoi();
	for (UInt_t  it=3;it<tokens.size();it++) {
	  HV_Value[CrateNum][SlotNum][it-3] = tokens[it].Atof();
	}
      }
    }
    file_hv.close();
  } else {
    cout << " could not open : " << HVfileName << endl;
  }
}//


void TweakHV(Int_t nrun, Double_t HVShift) {
  SetHVMap();
  string OutFile = Form("hv_set/run_%d_hv_newset.set",nrun);
  ReadHV(nrun);
  ReadAlpha();
  ReadGain();
  ofstream outfile_hv;
  outfile_hv.open(OutFile);
  TString CrateName[2] = {"rpi17:2001","rpi18:2001"};
  for (Int_t nc=0;nc<2;nc++) {
    for (Int_t ns=0;ns<16;ns++) {
      outfile_hv << CrateName[nc] << " S" << ns << " " << "DV" ;
      for (Int_t nch=0;nch<12;nch++) {
	if (HV_Block[nc][ns][nch] == -1 ) {
	  outfile_hv << " " <<  HV_Value[nc][ns][nch] ;
	} else {
	  int blk = HV_Block[nc][ns][nch];
	  double alphaINV = 1./gAlpha.at(blk);
	  double gain = adcGain.at(blk); // Gain = (New/target( or, old))
	  double HV_new = HV_Value[nc][ns][nch]*pow( gain,alphaINV );
	  outfile_hv << " " << HV_new;
	}
      }
      outfile_hv << endl;
    }
  }
  outfile_hv.close();     
}//

void ReadAlpha(){
  string InFile = "../golden/golden_alpha_sh.txt";
  ifstream infile_data;
  infile_data.open(InFile);
  TString currentline;
  if (infile_data.is_open() ) {
    cout << " Reading alphas for SH : " << InFile << endl;
    TString temp;
    while( currentline.ReadLine( infile_data ) ){
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();
      if( ntokens > 1 ){
	// temp = ( (TObjString*) (*tokens)[0] )->GetString();
	// elemID.push_back( temp.Atof() );
	temp = ( (TObjString*) (*tokens)[1] )->GetString();
	Double_t alpha = temp.Atof();
	gAlpha.push_back( alpha );
      }
      delete tokens;
    }
    infile_data.close();
  } else {
    cout << " No file : " << InFile << endl;
  }
}//

void ReadGain(Int_t nrun){
  string InFile Form("Gain/run_%d_adcGain.txt",nrun);
  ifstream infile_data;
  infile_data.open(InFile);
  TString currentline;
  if (infile_data.is_open() ) {
    cout << " Reading adc gains for SH : " << InFile << endl;
    TString temp;
    while( currentline.ReadLine( infile_data ) ){
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();
      if( ntokens > 1 ){
	temp = ( (TObjString*) (*tokens)[0] )->GetString();
	// elemID.push_back( temp.Atof() );
	// temp = ( (TObjString*) (*tokens)[1] )->GetString();
	Double_t ADCgain = temp.Atof();
	adcGain.push_back( ADCgain );
      }
      delete tokens;
    }
    infile_data.close();
  } else {
    cout << " No file : " << InFile << endl;
  }
}//
