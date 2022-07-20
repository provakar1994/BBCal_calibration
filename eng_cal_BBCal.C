#include <iostream>
#include <sstream>
#include <fstream>
#include "TChain.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "gmn_tree.C"

const Int_t ncell = 241; // SH(189) + PS(52), Convention: 0-188: SH; 189-240: PS.
const Int_t kNcolsSH = 7;  // SH columns
const Int_t kNrowsSH = 27;  // SH rows
const Int_t kNcolsPS = 2; // PS columns
const Int_t kNrowsPS = 26; // PS rows

const Double_t Mp = 0.938; // GeV

void ReadGain(TString,bool);
bool SHorPS = 1; // SH = 1, PS = 0
Double_t oldADCgainSH[kNcolsSH*kNrowsSH] = {0.};
Double_t oldADCgainPS[kNcolsPS*kNrowsPS] = {0.};  

void eng_cal_BBCal(const char *configfilename, int iter)
{
  TString outFile = Form("hist/eng_cal_BBCal_%d.root",iter);
  TFile *fout = new TFile(outFile,"RECREATE");
  TChain *C = new TChain("T");
  gmn_tree *T = new gmn_tree(C);

  Int_t Nmin = 10;
  Double_t minMBratio = 0.1;
  Double_t E_beam = 0.;

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }   
  } 
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_beam" ){
  	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
  	E_beam = sval.Atof();
      }
      if( skey == "Min_Event_Per_Channel" ){
  	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
  	Nmin = sval.Atof();
      }
      if( skey == "Min_MB_Ratio" ){
  	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
  	minMBratio = sval.Atoi();
      }
    } 
    delete tokens;
  }

  fout->cd();
  // Physics histos
  TH1D* h_deltaE = new TH1D("h_deltaE","",100,-1.5,1.5);
  TH1D* h_clusE = new TH1D("h_clusE","",100,0.,2.);
  TH2D* h_corPandAng = new TH2D("h_corPandAng","",100,30,60,100,0.4,1.2);
  TH1D* h_W = new TH1D("h_W","W distribution",200,0.7,1.6);
  TH1D* h_Q2 = new TH1D("h_Q2","Q2 distribution",40,0.,4.);
  
  TMatrixD M(ncell,ncell);
  TVectorD B(ncell);
  
  Int_t events_per_cell[ncell];
  memset(events_per_cell, 0, ncell*sizeof(int));
  
  Double_t E_e = 0;
  Double_t p_rec = 0.;
  Double_t A[ncell] = {0.};
  Double_t oldConstants_sum[ncell] = {1.};
  Int_t oldConstants_nevent[ncell] = {0};
  bool badcells[ncell]; // Cells that have events less than Nmin

  Long64_t Nevents = C->GetEntries();  
  for(Long64_t nevent = 0; nevent<Nevents; nevent++){
    if( nevent%1000 == 0){
      cout << nevent << "/" << Nevents << endl;
    }
    T->GetEntry(nevent);
    
    E_e = 0;
    memset(A, 0, ncell*sizeof(double));

    // Choosing track with least chi2 
    Double_t chi2min = 1000.;
    Int_t tr_min = -1;
    for(Int_t tr = 0; tr<T->bb_tr_n; tr++){
      if(T->bb_tr_chi2[tr]<chi2min){
	chi2min = T->bb_tr_chi2[tr];
	tr_min = tr;
      }
    }
    p_rec = (T->bb_tr_p[tr_min]); 
    E_e = p_rec; // Neglecting e- mass. 

    if(T->bb_tr_p[tr_min]==0 || E_e==0 || tr_min<0) continue;

    Double_t P_ang = 57.3*TMath::ACos(T->bb_tr_pz[tr_min]/T->bb_tr_p[tr_min]);
    Double_t Q2 = 4.*E_beam*p_rec*pow( TMath::Sin(P_ang/57.3/2.),2. );
    Double_t W2 = Mp*Mp + 2.*Mp*(E_beam-p_rec) - Q2;

    h_Q2->Fill(Q2);
    if(W2>0.) h_W->Fill( TMath::Sqrt(W2) );

    // Choosing only events which had clusters in both PS and SH. Also avoiding PS cluster with zero size (i.e. idblk = -1)
    if(T->bb_sh_nclus==0 || T->bb_ps_nclus==0 ) continue;
   
    if( T->bb_tr_tg_th[tr_min]>-0.15&&T->bb_tr_tg_th[tr_min]<0.15 &&T->bb_tr_tg_ph[tr_min]>-0.3&&T->bb_tr_tg_ph[tr_min]<0.3 ){ //&& fabs(W-0.945446)<.022498 ){ //cut on tracks and W
      Int_t cl_max = -1;
      Double_t nblk = -1.;
      Double_t Emax = -10.;

      // ****** Shower ******
      // Loop over all the clusters first: select highest energy
      // for(int ihit=0; ihit<T->Ndata_bb_sh_a_p; ihit++ ){
      // 	Int_t blkID = T->bb_sh_adcelemID[ihit];
      // 	Double_t a_c = T->bb_sh_a_c[ihit];
      // 	Double_t a_p = T->bb_sh_a_p[ihit];
      // 	if(a_c>0&&a_p>0){
      // 	  oldConstants_sum[blkID] += a_c/a_p; // a_c = (adc.gain) * a_p
      // 	  oldConstants_nevent[blkID]++;
      // 	}
      // }

      for(Int_t cl = 0; cl<T->bb_sh_nclus; cl++){
	if(T->bb_sh_clus_e[cl]>Emax){
	  Emax = T->bb_sh_clus_e[cl];
	  cl_max = cl;
	}
      }

      // Reject events with max on the edge
      if(T->bb_sh_clus_row[cl_max]==0 || T->bb_sh_clus_row[cl_max]==26 ||
	 T->bb_sh_clus_col[cl_max]==0 || T->bb_sh_clus_col[cl_max]==6) continue;
    
      // Loop over all the blocks in main cluster and fill in A's
      nblk = T->bb_sh_clus_nblk[cl_max];
      for(Int_t blk = 0; blk<nblk; blk++){
	//Int_t blkID = int(T->bb_sh_clus_blk_row[blk]*kNcolsSH + T->bb_sh_clus_blk_col[blk]);
	Int_t blkID = int(T->bb_sh_clus_blk_id[blk]);
	A[blkID] += T->bb_sh_clus_blk_e[blk];
	events_per_cell[ blkID ]++; 
      }
    
      // ****** PreShower ******
      // for(int ihit=0; ihit<T->Ndata_bb_ps_a_p; ihit++ ){
      // 	Int_t blkID = T->bb_ps_adcelemID[ihit];
      // 	Double_t a_c = T->bb_ps_a_c[ihit];
      // 	Double_t a_p = T->bb_ps_a_p[ihit];
      // 	if(a_c>0&&a_p>0){ // Need this condition for digitized data
      // 	  oldConstants_sum[189+blkID] += a_c/a_p;  // a_c = (adc.gain) * a_p
      // 	  oldConstants_nevent[189+blkID]++;     
      // 	}
      // }

      nblk = T->bb_ps_clus_nblk[0];
      for(Int_t blk=0; blk<nblk; blk++){
	//Int_t blkID = int(T->bb_ps_clus_blk_row[blk]*kNcolsPS + T->bb_ps_clus_blk_col[blk]);
	Int_t blkID = int(T->bb_ps_clus_blk_id[blk]);
	A[189+blkID]+=T->bb_ps_clus_blk_e[blk];
	events_per_cell[189+blkID]++;
      }

      // Let's fill some interesting histograms
      Double_t clusEngBBCal = T->bb_sh_e + T->bb_ps_e;
      h_deltaE->Fill( 1.-(clusEngBBCal/p_rec) );
      h_clusE->Fill( clusEngBBCal );
      h_corPandAng->Fill( P_ang, p_rec );

      // Let's costruct the matrix
      for(Int_t icol = 0; icol<ncell; icol++){
	B(icol)+= A[icol];
	for(Int_t irow = 0; irow<ncell; irow++){
	  M(icol,irow)+= A[icol]*A[irow]/E_e;
	} 
      }   
    }
  }
  
  // B.Print();  
  // M.Print();
  
  //Diagnostic histograms
  TH1D *h_constRatio_SH = new TH1D("h_constRatio_SH","",189,0,189);
  TH1D *h_constChan_SH = new TH1D("h_constChan_SH","",189,0,189);
  TH1D *h_neventChan_SH = new TH1D("h_neventChan_SH","",189,0,189);
  TH2D *h_constDV_SH = new TH2D("h_constDV_SH","",kNcolsSH,1,kNcolsSH+1,kNrowsSH,1,kNrowsSH+1);
  TH1D *h_oldConstChan_SH = new TH1D("h_onlConstChan_SH","",189,0,189);

  TH1D *h_constRatio_PS = new TH1D("h_constRatio_PS","",52,0,52);
  TH1D *h_constChan_PS = new TH1D("h_constChan_PS","",52,0,52);
  TH1D *h_neventChan_PS = new TH1D("h_neventChan_PS","",52,0,52);
  TH2D *h_constDV_PS = new TH2D("h_constDV_PS","",kNcolsPS,1,kNcolsPS+1,kNrowsPS,1,kNrowsPS+1);
  TH1D *h_oldConstChan_PS = new TH1D("h_onlConstChan_PS","",52,0,52);

  // Leave the bad channel out of the calculation
  for(Int_t j = 0; j<ncell; j++){
    badcells[j]=false;
    if( events_per_cell[j]<Nmin || M(j,j)< minMBratio*B(j) ){
      B(j) = 1.;
      M(j, j) = 1.;
      for(Int_t k = 0; k<ncell; k++){
	if(k!=j){
	  M(j, k) = 0.;
	  M(k, j) = 0.;
	}
      }
      badcells[j]=true;
    }
  }  
  
  // Getting coefficients
  TMatrixD M_inv = M.Invert();
  TVectorD Coeff = M_inv*B;
  
  // Let's read in old coefficients for both SH and PS
  TString adcGain_SH = Form("Gain/eng_cal_gainCoeff_sh_%d.txt",iter-1);
  TString adcGain_PS = Form("Gain/eng_cal_gainCoeff_ps_%d.txt",iter-1);
  ReadGain(adcGain_SH,SHorPS);
  SHorPS = 0; // Setting the flag for PS
  ReadGain(adcGain_PS,SHorPS);

  // SH : Filling diagnostic histograms
  int cell = 0;
  adcGain_SH = Form("Gain/eng_cal_gainCoeff_sh_%d.txt",iter);
  ofstream adcGainSH_outData;
  adcGainSH_outData.open(adcGain_SH);
  for(Int_t shrow = 0; shrow<kNrowsSH; shrow++){
    for(Int_t shcol = 0; shcol<kNcolsSH; shcol++){
      //Double_t oldConst = oldConstants_sum[cell]/(double)oldConstants_nevent[cell];
      Double_t oldConst = oldADCgainSH[shrow*kNcolsSH+shcol];
      if(!badcells[cell]){
	h_constRatio_SH->Fill( cell, Coeff(cell) );
	h_constChan_SH->Fill( cell, Coeff(cell)*oldConst );
	h_neventChan_SH->Fill( cell, events_per_cell[cell] );
	h_constDV_SH->Fill( shcol+1, shrow+1, Coeff(cell)*oldConst );
	h_oldConstChan_SH->Fill( cell, oldConst );

	cout << Coeff(cell) << "  ";
	adcGainSH_outData << Coeff(cell)*oldConst << " ";
      }else{
	h_constRatio_SH->Fill( cell, 1. );
	h_constChan_SH->Fill( cell, oldConst );
	h_neventChan_SH->Fill( cell, events_per_cell[cell] );
	h_constDV_SH->Fill( shcol+1, shrow+1, oldConst );
	h_oldConstChan_SH->Fill( cell, oldConst );

	cout << 1. << "  ";
	adcGainSH_outData << 1. << " ";
     }
      cell++;
    }
    cout << endl;
    adcGainSH_outData << endl;
  }
  cout << endl;

  // PS : Filling diagnostic histograms
  adcGain_PS = Form("Gain/eng_cal_gainCoeff_ps_%d.txt",iter);
  ofstream adcGainPS_outData;
  adcGainPS_outData.open(adcGain_PS);
  for(Int_t psrow = 0; psrow<kNrowsPS; psrow++){
    for(Int_t pscol = 0; pscol<kNcolsPS; pscol++){
      //Double_t oldConst =  oldConstants_sum[cell]/(double)oldConstants_nevent[cell];
      Double_t oldConst = oldADCgainPS[psrow*kNcolsPS+pscol];
      if(!badcells[cell]){
	h_constRatio_PS->Fill( ncell-cell, Coeff(cell) );
	h_constChan_PS->Fill( ncell-cell, Coeff(cell)*oldConst );
	h_neventChan_PS->Fill( ncell-cell, events_per_cell[cell] );
	h_constDV_PS->Fill( pscol+1, psrow+1, Coeff(cell)*oldConst );
	h_oldConstChan_PS->Fill( ncell-cell, oldConst );

	cout << Coeff(cell) << "  ";
	adcGainPS_outData << Coeff(cell)*oldConst << " ";
      }else{
	h_constRatio_PS->Fill( ncell-cell, 1. );
	h_constChan_PS->Fill( ncell-cell, oldConst );
	h_neventChan_PS->Fill( ncell-cell, events_per_cell[cell] );
	h_constDV_PS->Fill( pscol+1, psrow+1, oldConst );
	h_oldConstChan_PS->Fill( ncell-cell, oldConst );

	cout << 1. << "  ";
	adcGainPS_outData << 1. << " ";
     }
      cell++;
    }
    cout << endl;
    adcGainPS_outData << endl;
  }

  fout->Write();
  adcGainSH_outData.close();
  adcGainPS_outData.close();
}


void ReadGain( TString adcGain, bool SHorPS ){
  ifstream adcGain_data;
  adcGain_data.open(adcGain);
  string readline;
  Int_t elemID=0;
  if( adcGain_data.is_open() ){
    cout << " Reading ADC gain file : "<< adcGain << endl;
    while(getline(adcGain_data,readline)){
      istringstream tokenStream(readline);
      string token;
      char delimiter = ' ';
      while(getline(tokenStream,token,delimiter)){
	TString temptoken=token;
	if(SHorPS){
	  oldADCgainSH[elemID] = temptoken.Atof();
	}else{
	  oldADCgainPS[elemID] = temptoken.Atof();
	}
	elemID++;
      }
    }
  }else{
    cout << " No file : " << adcGain << endl;
  }
  adcGain_data.close();
}

