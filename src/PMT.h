
#ifndef PMT_H
#define PMT_H

#include "PMType.h"
#include "TObject.h"
#include "TH1.h"

#include "Pedestal.h"

class PMT : public TObject {
private:
   TH1D *spectrum;
   Int_t nbins;
   Double_t min;
   Double_t max;
   Pedestal ped;
   TF1 res;

public:
   PMT();
   PMT(Int_t _nbins, Double_t _min, Double_t _max, Pedestal _ped, PMType::Response _res, Double_t *_params);
   virtual ~PMT();

   void GenSpectrum(Int_t ntot, Double_t mu);
   TH1D *GetSpectrum() { return spectrum; };
   void DrawSpectrum() { spectrum->Draw("PEZ"); };

   ClassDef(PMT, 1)
};

#endif
