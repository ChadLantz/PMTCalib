#ifndef NUMINTEGRATION_H
#define NUMINTEGRATION_H

#include <math.h>

#include "PMType.h"
#include "TObject.h"
#include "TGraph.h"

#include "SPEResponseFactory.h"



class NumIntegration : public TObject
{
 private:

  Int_t nbins;
  
  Double_t xmin;
  Double_t xmax;

  Double_t step;
      
  unsigned int N;
    
  std::vector<Double_t> xvalues;
  std::vector<Double_t> yvalues;
    
  Double_t edge;

  TGraph *gr;
  
 public:
  
  NumIntegration();
  
  virtual ~NumIntegration();
  
  NumIntegration( Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Response _sper, Double_t *_params);

  TF1 spef;
  
  Double_t wbin;
  
  Double_t Norm;
  
  Double_t Q0;
  Double_t s0;

  Double_t mu;
  
  void CalculateValues();
  Double_t GetValue( Double_t xx );

  TGraph* GetGraph();

  
  ClassDef( NumIntegration, 1 )
        
};

#endif
