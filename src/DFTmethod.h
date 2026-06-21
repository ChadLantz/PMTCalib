#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include <RtypesCore.h>
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
  const UInt_t nPars{8};
  Double_t parCache[8];
  TGraph *gr;
  
 public:
  DFTmethod();
  DFTmethod( Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef );
  DFTmethod(DFTmethod &);
  DFTmethod operator=(const DFTmethod &);
  virtual ~DFTmethod();
  
  Double_t fftPhase( Double_t vy, Double_t vz );
  void CalculateValues();
  Double_t GetValue( Double_t xx );
  Double_t Eval(Double_t *xx, Double_t *pars);
  double operator() (double *xx, double *pars){ return Eval(xx, pars); }

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
