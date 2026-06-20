
#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include <math.h>

#include "TObject.h"
#include "TH1.h"

class Deconvolution : public TObject
{
 private:
  Double_t Q0;
  Double_t dQ0;
  Double_t s0;
  Double_t ds0;
  Double_t x1;
  Double_t mu;
  Float_t mu_bf;
  
 public:
  Deconvolution();
  Deconvolution( Double_t _Q0, Double_t _dQ0, Double_t _s0, Double_t _ds0, Double_t _mu );  
  virtual ~Deconvolution();
  
  TH1D* CleanUps( TH1D *h );
  TH1D* Deconvolute( TH1D* h, Double_t _Q0, Double_t _s0, Double_t _mu );
  TH1D* RunSingle( TH1D* h, Double_t _Q0, Double_t _s0 );
  TH1D* Run( TH1D* h, Int_t ntoys );    
  Double_t fftPhase( Double_t vy, Double_t vz );
  Double_t GridMu( TH1D *h, Double_t _Q0, Double_t _s0 );
    
  ClassDef( Deconvolution, 1 )
        
};

#endif
