#ifndef FADC_DATA_H
#define FADC_DATA_H
#include <TPaveStats.h>

const Int_t MAX_FADC_SAMPLES = 250;
const Int_t MAX_FADC_DATA_MODULES = 288;
const Int_t MAX_FADC_DATA_TDC_MODULES = 288;
namespace fadc_datat {
  Double_t samps[MAX_FADC_DATA_MODULES*MAX_FADC_SAMPLES+1000];
  Double_t nsamps[MAX_FADC_DATA_MODULES+1000] = {0};
  Double_t row[MAX_FADC_DATA_MODULES+1000] = {0};
  Double_t col[MAX_FADC_DATA_MODULES+1000] = {0};
  Double_t samps_idx[MAX_FADC_DATA_MODULES+1000] = {0};
  Double_t a[MAX_FADC_DATA_MODULES+1000] = {0};
  Int_t ndata = 0;
  Double_t ledbit = -1;
  Double_t ledcount = 0;
  Double_t tdc[MAX_FADC_DATA_TDC_MODULES+100];
  Double_t amp[MAX_FADC_DATA_MODULES+100];
};

void fixStats()
{
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  if(ps) {
    ps->SetX1NDC(0.6);
    ps->SetY1NDC(0.55);
  }
}



#endif // FADC_DATA_H
