#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include <math.h>

#include "TObject.h"
#include "TGraph.h"

#include "SPEResponse.h"

class DFTmethod : public TObject
{
 private:
  UInt_t N;
  UInt_t M;
  Int_t nbins;
  Double_t xmin;
  Double_t xmax;
  Double_t step;
  Double_t edge;
  std::vector<Double_t> xvalues;
  std::vector<Double_t> yvalues;

  TGraph *gr;
  
 public:
  
  DFTmethod();
  DFTmethod( Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef );
  virtual ~DFTmethod();
  
  Double_t fftPhase( Double_t vy, Double_t vz );
  void CalculateValues();
  Double_t GetValue( Double_t xx );

  TGraph* GetGraph();
  TGraph* GetGraphN( Int_t n );

  SPEResponse spef;
  Double_t wbin;
  Double_t Norm;
  Double_t Q0;
  Double_t s0;
  Double_t mu;
  
  ClassDef( DFTmethod, 1 )
        
};

#endif
