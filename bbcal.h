#ifndef BBCAL_H
#define BBCAL_H
#include <TPaveStats.h>

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_BBCAL_MODULES = 196;
const Int_t MAX_TRIGGER_MODULES = 16;
const Int_t MAX_BBCAL_TDC_MODULES = 196;

namespace bbcalt {
  // BBCAL shower core vars
  Double_t samps_sh[MAX_BBCAL_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t nsamps_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t row_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t col_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t samps_idx_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_p_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t ped_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_amp_sh[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_amp_p_sh[MAX_BBCAL_MODULES+1000] = {0};
  Int_t ndata_sh = 0;
  Double_t tdc_sh[MAX_BBCAL_TDC_MODULES+100];
  // BBCAL preshower core vars
  Double_t samps_ps[MAX_BBCAL_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t nsamps_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t row_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t col_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t samps_idx_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_p_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t ped_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_amp_ps[MAX_BBCAL_MODULES+1000] = {0};
  Double_t a_amp_p_ps[MAX_BBCAL_MODULES+1000] = {0};
  Int_t ndata_ps = 0;
  Double_t tdc_ps[MAX_BBCAL_TDC_MODULES+100];
  /*
  //Cluster vars
  Double_t nclus[MAX_BBCAL_MODULES+1000] = {0};
  Double_t nblk[MAX_BBCAL_MODULES+1000] = {0};
  Double_t cid[MAX_BBCAL_MODULES+1000] = {0};
  Double_t cblkid[MAX_BBCAL_MODULES+1000] = {0};
  Double_t crow[MAX_BBCAL_MODULES+1000] = {0};
  Double_t ccol[MAX_BBCAL_MODULES+1000] = {0};
  Double_t ce[MAX_BBCAL_MODULES+1000] = {0};
  Double_t ceblk[MAX_BBCAL_MODULES+1000] = {0};
  Double_t cnblk[MAX_BBCAL_MODULES+1000] = {0};
  Double_t cblke[MAX_BBCAL_MODULES+1000] = {0};
  //Trigger generic detector vars
  Double_t Tsamps[MAX_TRIGGER_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t Tnsamps[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Tsamps_idx[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_p[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_amp_p[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Ta_amp[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t TelemID[MAX_TRIGGER_MODULES+1000] = {0};
  Double_t Tcol[MAX_TRIGGER_MODULES+1000] = {0};
  Int_t Tndata = 0;
  */
}

void fixStats()
{
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  if(ps) {
    ps->SetX1NDC(0.6);
    ps->SetY1NDC(0.55);
  }
}



#endif // BBCAL_H
