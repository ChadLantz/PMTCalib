
#ifndef PEDESTAL_H
#define PEDESTAL_H

#include "TObject.h"
#include "TF1.h"
#include "TH1.h"

class Pedestal : public TObject
{
 private:
  TF1 *pedfunc;
  TF1 *pedfit;
      
 public:
  Pedestal();
  Pedestal( Double_t _Q0, Double_t _s0 );
  virtual ~Pedestal();

  Double_t GenQ();
  void LocatePedestal( TH1 *hspec, Double_t _Q0, Double_t _s0 );
  
  Int_t status;
  Double_t chi2;
  Double_t wbin;
  Double_t Q0;
  Double_t s0;
  Double_t dQ0;
  Double_t ds0;
  
  ClassDef( Pedestal, 1 )
        
};

#endif
