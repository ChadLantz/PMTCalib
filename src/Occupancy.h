
#ifndef OCCUPANCY_H
#define OCCUPANCY_H

#include "TObject.h"
#include "TH1.h"

class Occupancy : public TObject {
private:
   Double_t Q0;
   Double_t s0;

public:
   Occupancy();
   Occupancy(Double_t _Q0, Double_t _s0);
   virtual ~Occupancy();

   Double_t Gauss1(Float_t x);
   Float_t FindG(TH1D *_h, Float_t f);

   Float_t mu_bf;

   ClassDef(Occupancy, 1)
};

#endif
