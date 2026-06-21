
#ifndef SPEFITTER_H
#define SPEFITTER_H

#include "TObject.h"
#include "TF1.h"
#include "TH1.h"

#include "NumIntegration.h"
#include "DFTmethod.h"
#include "PMTModel.h"

#include "Minuit2/Minuit2Minimizer.h"

class SPEFitter : public TObject {
private:
   NumIntegration num;
   DFTmethod dft;
   PMTModel mod;
   ROOT::Minuit2::Minuit2Minimizer *mMOD;
   TF1 *ped_func;

public:
   SPEFitter() = default;
   virtual ~SPEFitter() = default;

   ROOT::Minuit2::Minuit2Minimizer *mNum;
   ROOT::Minuit2::Minuit2Minimizer *mFFT;
   Int_t fit_status;
   Double_t vals[20];
   Double_t errs[20];
   Double_t ndof;
   Double_t chi2r;

   Double_t FindMu(TH1D *hspec, Double_t _Q0, Double_t _s0);
   Double_t FindG(TH1 *hspec, Double_t _Q0, Double_t _mu);

   void SetNummethod(NumIntegration _num);
   void FitwNummethod(TH1 *hspec);

   void SetDFTmethod(DFTmethod _dft);
   void FitwDFTmethod(TH1 *hspec);

   void SetPMTModel(PMTModel _mod);
   void FitwPMTModel(TH1 *hspec);

   ClassDef(SPEFitter, 1)
};

#endif
