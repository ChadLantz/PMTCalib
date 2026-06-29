#include "PMT.h"


#include "SPEResponse.h"
#include "TRandom.h"

ClassImp(PMT)

   PMT::PMT()
{
}

PMT::~PMT() {}

PMT::PMT(Int_t _nbins, Double_t _min, Double_t _max, Double_t Q0, Double_t s0, PMType::Response _res, Double_t *_params)
: m_nbins(_nbins),
   m_min(_min),
   m_max(_max),
   m_ped(Q0, s0),
   m_res(_res)
{
   m_res.SetParameters(_params);
   m_spefunc = new TF1("pmt", m_res, m_min, m_max, m_res.NPar());
   m_spefunc->SetParameters(_params);
}

TH1D *PMT::GenSpectrum(Int_t ntot, Double_t mu)
{
   TH1D *spectrum = new TH1D("hspectrum", "PMT spectrum; Charge [AU]; Entries", m_nbins, m_min, m_max);
   spectrum->SetMarkerStyle(20);
   spectrum->SetMarkerSize(0.75);
   spectrum->SetLineColor(kBlack);
   spectrum->SetMarkerColor(kBlack);
   spectrum->Reset();
   gRandom->SetSeed(0);

   for (Int_t i = 0; i < ntot; i++) {
      Double_t q = m_ped.GenQ();
      Int_t npe = gRandom->Poisson(mu);
      for (Int_t j = 0; j < npe; j++) {
         q += m_spefunc->GetRandom();
      }

      spectrum->Fill(q);
   }

   return spectrum;
}
