#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include "SPEResponse.h"

#include "RtypesCore.h"
#include "TObject.h"
#include "TGraph.h"

struct FFTWState;
class DFTmethod : public TObject {

public:
   DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef);
   virtual ~DFTmethod() = default;

   Double_t Eval(Double_t *xx, Double_t *pars);
   Double_t operator()(Double_t *xx, Double_t *pars) { return Eval(xx, pars); }
   void SetParameters(Double_t *pars);
   TGraph *GetGraph();
   TGraph *GetGraphN(Int_t n);
   UInt_t GetNCalls(){ return nCalls; };

private:
   void CalculateValues();

   UInt_t nPars;
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
   SPEResponse spef;

   ClassDef(DFTmethod, 1)
};

#endif
