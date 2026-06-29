#ifndef PMT_H
#define PMT_H

#include "PMType.h"
#include "SPEResponse.h"
#include "TH1.h"

#include "Pedestal.h"
class PMT {
public:
   PMT();
   PMT(const PMT &other);
   PMT(Int_t _nbins, Double_t _min, Double_t _max, Double_t Q0, Double_t s0, PMType::Response _res, Double_t *_params);
   virtual ~PMT();

   TH1D *GenSpectrum(Int_t ntot, Double_t mu);

private:
   Int_t m_nbins;
   Double_t m_min;
   Double_t m_max;
   Pedestal m_ped;
   SPEResponse m_res;
   TF1 *m_spefunc{nullptr};
};

#endif
