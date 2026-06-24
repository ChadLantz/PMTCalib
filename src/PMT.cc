#include "PMT.h"


#include "SPEResponseFactory.h"
#include "TRandom.h"

ClassImp(PMT)

   PMT::PMT()
{
}

PMT::~PMT() {}

PMT::PMT(Int_t _nbins, Double_t _min, Double_t _max, Pedestal _ped, PMType::Response _res, Double_t *_params)
{
   nbins = _nbins;
   min = _min;
   max = _max;

   ped = _ped;
   SPEResponseFactory sperf;
   res = sperf.Build(_res, _params);

   spectrum = new TH1D("hspectrum", "PMT spectrum; Charge [AU]; Entries", nbins, min, max);

   spectrum->SetMarkerStyle(20);
   spectrum->SetMarkerSize(0.75);
   spectrum->SetLineColor(kBlack);
   spectrum->SetMarkerColor(kBlack);
}

void PMT::GenSpectrum(Int_t ntot, Double_t mu)
{
   spectrum->Reset();

   gRandom->SetSeed(0);

   for (Int_t i = 0; i < ntot; i++) {
      Double_t q = ped.GenQ();

      Int_t npe = gRandom->Poisson(mu);

      for (Int_t j = 0; j < npe; j++) {
         q += res.GetRandom();
      }

      spectrum->Fill(q);
   }

   return;
}
