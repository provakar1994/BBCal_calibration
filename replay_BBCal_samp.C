#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"

#include "THaShower.h"
#include "THaEvent.h"
#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"

#include "SBSBigBite.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"
// Simple example replay script
//
void replay(int run_number = 124, uint nev = -1, uint nseg = 0)
{
  //load SBS-offline
  //gSystem->Load("libsbs.so");
  //--- Define the experimental configuration, i.e. spectrometers, detectors ---

  //THaHRS* bb = new THaHRS("R", "Right HRS" );
  //hrs->AddDetector( new THaVDC("vdc", "Vertical Drift Chambers") );
  // gHaApps->Add(hrs);
  
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new THaShower("ps", "BigBite preshower") );
  /*
  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  //ts->SetDataOutputLevel(0);
  ts->SetStoreRawHits(kTRUE);
  bigbite->AddDetector( ts );
  */
  
  SBSBBShower* ps= new SBSBBShower("ps", "BigBite preshower");
  //ts->SetDataOutputLevel(0);
  ps->SetStoreRawHits(kTRUE);
  bigbite->AddDetector( ps );
  
  /*
  SBSBBShower* sh= new SBSBBShower("sh", "BigBite shower");
  //ts->SetDataOutputLevel(0);
  sh->SetStoreRawHits(kTRUE);
  bigbite->AddDetector( sh );
  */
  SBSGenericDetector* trig= new SBSGenericDetector("trig","BigBite shower trig");
  trig->SetModeADC(SBSModeADC::kWaveform);
  bigbite->AddDetector( trig );
  //SBSGenericDetector* tdctrig= new SBSGenericDetector("tdctrig","BigBite shower TDC trig");
  // tdctrig->SetModeADC(SBSModeADC::kNone);emacs 
  //tdctrig->SetModeTDC(SBSModeTDC::kTDCSimple);
  //bigbite->AddDetector( tdctrig );
  //  bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  // bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  gHaApps->Add(bigbite);
  
  // Ideal beam (perfect normal incidence and centering)
  // THaIdealBeam* ib = new THaIdealBeam("IB", "Ideal beam");
  // gHaApps->Add(ib);

  //--- Set up the run we want to replay ---

  // This often requires a bit of coding to search directories, test
  // for non-existent files, etc.
  //bool seg_ok = true;

  TString exp = "bbpreshower_samp";

  string firstname = "e1209019_trigtest_%d";

  string endname = Form(".evio.%d",nseg);

  string combined(string(firstname)+endname);
   const char* RunFileNamePattern = combined.c_str();
  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/sbs"));

  THaAnalyzer* analyzer = new THaAnalyzer;
  THaEvent* event = new THaEvent;
  THaRun* run = 0;
  int seg = 0;
  bool seg_ok = true;
  while(seg_ok) {
    TString data_fname;
    //data_fname = TString::Format("%s/ts_bbshower_%d.evio.%d",getenv("DATA_DIR"),run_number,seg);
    data_fname = TString::Format("%s/e1209019_trigtest_%d.evio.%d","../data",run_number,seg);
 //new THaRun( pathList, Form(RunFileNamePattern, run_number) );
    std::cout << "Looking for segment " << seg << " file " << data_fname.Data() << std::endl;
    if( gSystem->AccessPathName(data_fname)) {
      seg_ok = false;
      std::cout << "Segment " << seg << " not found. Exiting" << std::endl;
      continue;
    }
   run = new THaRun(data_fname);
    run->SetLastEvent(nev);

  run->SetDataRequired(0);//for the time being
   run->SetDate(TDatime());
  
  

  analyzer->SetEvent( event );
  TString out_dir = "/adaqfs/home/a-onl/sbs/Rootfiles";
  if( out_dir.IsNull() )  out_dir = ".";
  TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);

  analyzer->SetOutFile( out_file );
  
  analyzer->SetCutFile( "replay.cdef" );
  analyzer->SetOdefFile( "replay.odef" );

  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();

  //--- Analyze the run ---
  // Here, one could replay more than one run with
  // a succession of Process calls. The results would all go into the
  // same ROOT output file

  run->Print();
  

  analyzer->Process(run);
    // Cleanup this run segment
    delete run;
    
    seg++; // Increment for next search
 }

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();

  // Open the ROOT file so that a user doing interactive analysis can 
  // immediately look at the results. Of course, don't do this in batch jobs!
  // To close the file later, simply type "delete rootfile" or just exit

  //TFile* rootfile = new TFile( out_file, "READ" );
}

