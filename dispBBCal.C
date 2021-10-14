#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "bbcal.h"

const Int_t kPSNrows = 26;
const Int_t kSHNrows = 27;

const Int_t kPSNcols = 2;
const Int_t kSHNcols = 7;

const Int_t kNumModules_PS = kPSNrows*kPSNcols;
const Int_t kNumModules_SH = kSHNrows*kSHNcols;

const Int_t DISP_MIN_SAMPLE = 0;
const Int_t DISP_MAX_SAMPLE = 50;

const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
const Int_t numSamples = 50;

const Int_t minADC = 0;
const Int_t maxADC = 4000;
const Int_t kCanvSize = 100;

std::string user_input;

Int_t gCurrentEntry = -1;

TChain *T = 0;
Int_t foundModules = 0;
TCanvas *canvas = 0;

TCanvas *subCanv[4];

void clicked_displayEntryButton();
void clicked_displayNextButton();

namespace hcalgui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  TGTextButton *displayEntryButton;
  TGTextButton *displayNextButton;
  TGNumberEntry *entryInput;
  //TGLabel *ledLabel;

  TRootEmbeddedCanvas *canv[4];

  TGCompositeFrame* AddTabSub(Int_t sub) {
    

    if ( sub < 4 ) 
      tf = fTab->AddTab(Form("BBCal SH fADC Sub%d",sub+1));

    if ( sub > 3 )
      tf = fTab->AddTab(Form("BBCal PS fADC Sub%d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("BBcalSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot(), 1000, 900);
      frame1 = new TGHorizontalFrame(main, 150, 20, kFixedWidth);
      //ledLabel = new TGLabel(frame1,"LED Bit:    , Count:      ");
      displayEntryButton = new TGTextButton(frame1,"&Display Entry","clicked_displayEntryButton()");
      entryInput = new TGNumberEntry(frame1,0,5,-1,TGNumberFormat::kNESInteger);
      displayNextButton = new TGTextButton(frame1,"&Next Entry","clicked_displayNextButton()");
      exitButton = new TGTextButton(frame1, "&Exit", 
          "gApplication->Terminate(0)");
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
          kLHintsExpandX,2,2,2,2);
      //frame1->AddFrame(ledLabel,frame1LH);
      frame1->AddFrame(displayEntryButton,frame1LH);
      frame1->AddFrame(entryInput,frame1LH);
      frame1->AddFrame(displayNextButton,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      frame1->Resize(800, displayNextButton->GetDefaultHeight());
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));

      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

      // Create Tab1 (BBCal Sub1)

      for(Int_t i = 0; i < 8; i++) {
        tf = AddTabSub(i);
      }      

      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      for(Int_t i = 0; i < 4; i++) {
        subCanv[i] = canv[i]->GetCanvas();

	subCanv[i]->Divide(7,7,0.001,0.001);
        
      }
      for(Int_t i = 4; i < 8; i++) {
        subCanv[i] = canv[i]->GetCanvas();
        
	subCanv[i]->Divide(2,7,0.001,0.001);
        
      }
      
    }
  }
};


Double_t nhit = 0;
TH1F *gPShistos[kPSNrows][kPSNcols];
TH1F *gSHhistos[kSHNrows][kSHNcols];
Bool_t gPSsaturated[kPSNrows][kPSNcols];
Bool_t gSHsaturated[kSHNrows][kSHNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

bool is_number(const std::string& mystring)
{
  std::string::const_iterator it = mystring.begin();
  while (it != mystring.end() && std::isdigit(*it)) ++it;
  return !mystring.empty() && it == mystring.end();
}

void displayEvent(Int_t entry = -1)
{
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }
  
  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }
  
  T->GetEntry(gCurrentEntry);
  std::cout << "Displaying event " << gCurrentEntry << std::endl;
  //hcalgui::ledLabel->SetText(TString::Format("LED Bit: %02d, Count: %5d",Int_t(hcalt::ledbit),Int_t(hcalt::ledcount)));
  
  Int_t r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kPSNrows; r++) {
    for(c = 0; c < kPSNcols; c++) {
      gPShistos[r][c]->Reset("ICES M");
      gPSsaturated[r][c] = false;
    }
  }
  for(r = 0; r < kSHNrows; r++) {
    for(c = 0; c < kSHNcols; c++) {
      gSHhistos[r][c]->Reset("ICES M");
      gSHsaturated[r][c] = false;
    } 
  }  

  Double_t PSpeak[kPSNrows][kPSNcols];
  Double_t PSadc[kPSNrows][kPSNcols];
  //Double_t PStdc[kPSNrows][kSHNcols];

  Double_t SHpeak[kSHNrows][kSHNcols];
  Double_t SHadc[kSHNrows][kSHNcols];
  //Double_t SHtdc[kSHNrows][kSHNcols];
  
  for(r = 0; r < kSHNrows; r++) {
    for(c  = 0; c < kSHNcols; c++) {
      SHpeak[r][c] = 0.0;
      SHadc[r][c] = 0.0;
      //SHtdc[r][c] = 0.0;

    }
  }
  for(r = 0; r < kPSNrows; r++) {
    for(c  = 0; c < kPSNcols; c++) {
      PSpeak[r][c] = 0.0;
      PSadc[r][c] = 0.0;
      //PStdc[r][c] = 0.0;
    
    }
  }


  // Display shower
  for(Int_t m = 0; m < bbcalt::ndata_sh; m++) {
    r = bbcalt::row_sh[m];
    c = bbcalt::col_sh[m];
    SHpeak[r][c] = bbcalt::a_amp_sh[m];
    if(c < 0) {
      std::cerr << "Error: column negative." << std::endl;
      continue;
    }
    
    if(c >= kSHNcols)
      continue;
    
    // Fill adc, tdc, and cluster arrays
    idx = bbcalt::samps_idx_sh[m];
    n = bbcalt::nsamps_sh[m];
    SHadc[r][c] = bbcalt::a_p_sh[m];
    //SHtdc[r][c] = bbcalt::tdc_sh[m];
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      gSHhistos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,bbcalt::samps_sh[idx+s]); //Prepared for trigger samples to be added to the root tree
      if(SHpeak[r][c]<bbcalt::samps_sh[idx+s])
	SHpeak[r][c]=bbcalt::samps_sh[idx+s];
      if(SHpeak[r][c]>4095) {
	gSHsaturated[r][c] = true;
      }
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
	gSHhistos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  for(r = 0; r < kSHNrows; r++) {
    for(c = 0; c < kSHNcols; c++) {
      sub = r/7;
      subCanv[sub]->cd( (r%7)*kSHNcols + c + 1 );
      gSHhistos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g)",r+1,c+1,SHadc[r][c]));
      gSHhistos[r][c]->SetMaximum(60);
      if(gSHsaturated[r][c])
	gSHhistos[r][c]->SetLineColor(kRed+1);
      else
	gSHhistos[r][c]->SetLineColor(kBlue+1);
      gSHhistos[r][c]->Draw();
      gPad->Update();
      std::cout << "SH Peak at [" << r << "][" << c << "]=" << SHpeak[r][c];
    }
  }

  // Display preshower
  for(Int_t m = 0; m < bbcalt::ndata_ps; m++) {
    r = bbcalt::row_ps[m];
    c = bbcalt::col_ps[m];
    PSpeak[r][c] = bbcalt::a_amp_ps[m];
    if(c < 0) {
      std::cerr << "Error: column negative." << std::endl;
      continue;
    }
    
    if(c >= kPSNcols)
      continue;
    
    // Fill adc, tdc, and cluster arrays
    idx = bbcalt::samps_idx_ps[m];
    n = bbcalt::nsamps_ps[m];
    PSadc[r][c] = bbcalt::a_p_ps[m];
    //PStdc[r][c] = bbcalt::tdc_ps[m];
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      gPShistos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,bbcalt::samps_ps[idx+s]); //Prepared for trigger samples to be added to the root tree
      if(PSpeak[r][c]<bbcalt::samps_ps[idx+s])
	PSpeak[r][c]=bbcalt::samps_ps[idx+s];
      if(PSpeak[r][c]>4095) {
	gPSsaturated[r][c] = true;
      }
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
	gPShistos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  for(r = 0; r < kPSNrows; r++) {
    for(c = 0; c < kPSNcols; c++) {
      sub = r/7 + 4;
      subCanv[sub]->cd( (r%7)*kPSNcols + c + 1 );
      gPShistos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g)",r+1,c+1,PSadc[r][c]));
      gPShistos[r][c]->SetMaximum(60);
      if(gPSsaturated[r][c])
	gPShistos[r][c]->SetLineColor(kRed+1);
      else
	gPShistos[r][c]->SetLineColor(kBlue+1);
      gPShistos[r][c]->Draw();
      gPad->Update();
      std::cout << "PS Peak at [" << r << "][" << c << "]=" << PSpeak[r][c];
    }
  }

  std::cout << std::endl;
}
void clicked_displayNextButton()
{
  hcalgui::entryInput->SetIntNumber(++gCurrentEntry);
  displayEvent(gCurrentEntry);
}

void clicked_displayEntryButton()
{
  gCurrentEntry = hcalgui::entryInput->GetIntNumber();
  displayEvent(gCurrentEntry);
}

Int_t dispBBCal(Int_t run = 1198, Int_t event = -1)
{

  // Initialize function with user commands
  //bool layer;

  cout << "Enter run number for analysis." << endl;
  cin >> run;

  //cout << "Layer? Enter 0 for preshower or 1 for shower." << endl;
  //cin >> layer;

  //if( layer==0 ){
  //  string type = "preshower";
  //}else{
  //  string type = "shower";
  //}

  hcalgui::SetupGUI();
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleFontSize(0.08);  

  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("/volatile/halla/sbs/seeds/rootfiles/bbtotshower_samp_%d*.root",run));
    //T->Add(TString::Format("../Rootfiles/bbshower_%d*.root",run));
    T->SetBranchStatus("*",0);
    
    T->SetBranchStatus("bb.ps.*",1);
    T->SetBranchStatus("bb.sh.*",2);

    // Shower branches
    T->SetBranchAddress("bb.sh.nsamps",bbcalt::nsamps_sh);
    T->SetBranchAddress("bb.sh.a",bbcalt::a_sh);
    T->SetBranchAddress("bb.sh.a_p",bbcalt::a_p_sh);
    T->SetBranchAddress("bb.sh.a_amp",bbcalt::a_amp_sh);
    T->SetBranchAddress("bb.sh.a_amp_p",bbcalt::a_amp_p_sh);
    T->SetBranchAddress("bb.sh.ped",bbcalt::ped_sh);
    //T->SetBranchAddress("bb.sh.tdc",bbcalt::tdc_sh);
    T->SetBranchAddress("bb.sh.samps",bbcalt::samps_sh);
    T->SetBranchAddress("bb.sh.samps_idx",bbcalt::samps_idx_sh);
    T->SetBranchAddress("bb.sh.adcrow",bbcalt::row_sh);
    T->SetBranchAddress("bb.sh.adccol",bbcalt::col_sh);
    T->SetBranchStatus("Ndata.bb.sh.adcrow",1);
    T->SetBranchAddress("Ndata.bb.sh.adcrow",&bbcalt::ndata_sh);

    // Preshower branches
    T->SetBranchAddress("bb.ps.nsamps",bbcalt::nsamps_ps);
    T->SetBranchAddress("bb.ps.a",bbcalt::a_ps);
    T->SetBranchAddress("bb.ps.a_p",bbcalt::a_p_ps);
    T->SetBranchAddress("bb.ps.a_amp",bbcalt::a_amp_ps);
    T->SetBranchAddress("bb.ps.a_amp_p",bbcalt::a_amp_p_ps);
    T->SetBranchAddress("bb.ps.ped",bbcalt::ped_ps);
    //T->SetBranchAddress("bb.ps.tdc",bbcalt::tdc_ps);
    T->SetBranchAddress("bb.ps.samps",bbcalt::samps_ps);
    T->SetBranchAddress("bb.ps.samps_idx",bbcalt::samps_idx_ps);
    T->SetBranchAddress("bb.ps.adcrow",bbcalt::row_ps);
    T->SetBranchAddress("bb.ps.adccol",bbcalt::col_ps);
    T->SetBranchStatus("Ndata.bb.ps.adcrow",1);
    T->SetBranchAddress("Ndata.bb.ps.adcrow",&bbcalt::ndata_ps);

    /*
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchStatus("sbs.trig.*",2);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.a_p",hcalt::a_p);
    T->SetBranchAddress("sbs.hcal.a_amp",hcalt::a_amp);
    T->SetBranchAddress("sbs.hcal.a_amp_p",hcalt::a_amp_p);
    T->SetBranchAddress("sbs.hcal.ped",hcalt::ped);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.adcrow",hcalt::row);
    T->SetBranchAddress("sbs.hcal.adccol",hcalt::col);
    
    // Add trigger generic detector branches
    T->SetBranchAddress("sbs.trig.a_p",hcalt::Ta_p);
    T->SetBranchAddress("sbs.trig.adcelemID",hcalt::TelemID);
    T->SetBranchAddress("sbs.trig.nsamps",hcalt::Tnsamps);
    T->SetBranchAddress("sbs.trig.samps",hcalt::Tsamps);
    T->SetBranchAddress("sbs.trig.samps_idx",hcalt::Tsamps_idx);
    T->SetBranchAddress("sbs.trig.adccol",hcalt::Tcol);
    T->SetBranchAddress("sbs.trig.a_amp_p",hcalt::Ta_amp_p);
    T->SetBranchAddress("sbs.trig.a_amp",hcalt::Ta_amp);
    
    // Add clustering branches
    T->SetBranchAddress("sbs.hcal.clus.id",hcalt::cid);
    T->SetBranchAddress("sbs.hcal.clus.row",hcalt::crow);
    T->SetBranchAddress("sbs.hcal.clus.col",hcalt::ccol);
    T->SetBranchAddress("sbs.hcal.clus.e",hcalt::ce);
    T->SetBranchAddress("sbs.hcal.clus.eblk",hcalt::ceblk);

    T->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
    T->SetBranchAddress("Ndata.sbs.hcal.adcrow",&hcalt::ndata);

    T->SetBranchStatus("Ndata.sbs.trig.adccol",2);
    T->SetBranchAddress("Ndata.sbs.trig.adccol",&hcalt::Tndata);
    */
    
    
    std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    for(Int_t r = 0; r < kPSNrows; r++) {
      for(Int_t c = 0; c < kPSNcols; c++) {
	gPShistos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	gPSsaturated[r][c] = false;
      }
    }
    for(Int_t r = 0; r < kSHNrows; r++) {
      for(Int_t c = 0; c < kSHNcols; c++) {
	gSHhistos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	gSHsaturated[r][c] = false;
      }
    }
  }

  if( T->GetEntries()<=0 ){
    cout << "Error: no data." << endl;
    return 0;
  }

  gCurrentEntry = event;
  while( user_input != "q" ) {
    if(is_number(user_input)) {
      gCurrentEntry = std::stoi(user_input);
    } else {
      gCurrentEntry++;
    }
    displayEvent(gCurrentEntry);
    std::cout << "Display options: <enter> == next event, or q to stop." << std::endl;
    getline(std::cin,user_input);
  }

  return 0;
}

