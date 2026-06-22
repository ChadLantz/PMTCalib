#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include <RtypesCore.h>
#include <fftw3.h>
#include <math.h>

#include "TObject.h"
#include "TGraph.h"

#include "SPEResponse.h"

struct FFTWState;
class DFTmethod : public TObject {

public:
   DFTmethod();
   DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef);
   virtual ~DFTmethod();

   Double_t fftPhase(Double_t vy, Double_t vz);
   Double_t Eval(Double_t *xx, Double_t *pars);
   double operator()(double *xx, double *pars) { return Eval(xx, pars); }

   TGraph *GetGraph();
   TGraph *GetGraphN(Int_t n);

private:
   void CalculateValues();

   const UInt_t nPars{8};
   UInt_t N;
   UInt_t M;
   UInt_t nCalls;
   Int_t nbins;
   Double_t xmin;
   Double_t xmax;
   Double_t step;
   Double_t edge;
   Double_t wbin;
   Double_t Norm;
   Double_t Q0;
   Double_t s0;
   Double_t mu;
   Double_t parCache[8];
   TGraph *gr;
   SPEResponse spef;

   ClassDef(DFTmethod, 1)
};

#endif
