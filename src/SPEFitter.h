#ifndef SPEFITTER_H
#define SPEFITTER_H

#include "PMType.h"

#include "RtypesCore.h"
#include "Fit/FitResult.h"

#include <map>

class IModel;
class DFTmethod;
class PMTModel;
class NumIntegration;
class TH1;

class SPEFitter {
public:
   SPEFitter() = default;
   ~SPEFitter() = default;

   std::map<std::string, Double_t> GenerateSeeds(TH1 *hspec, const Double_t Q0, const Double_t s0);
   NumIntegration *CreateNumethod(TH1 *hspec, PMType::Response sper, const Double_t Q0, const Double_t s0);
   DFTmethod *CreateDFTmethod(TH1 *hspec, PMType::Response sper, const Double_t Q0, const Double_t s0);
   PMTModel *CreatePMTModel(TH1 *hspec, PMType::Model model, const Double_t Q0, const Double_t s0);
   ROOT::Fit::FitResult HybridMinimize(IModel *model, TH1 *hspec, Int_t maxItersGA = 150, Int_t maxItersSimplex = 1e5,
                                       Int_t maxItersMigrad = 1e8, Double_t tolSimplex = 1.0,
                                       Double_t tolMigrad = 1e-2);
   void SetVerbose(Int_t level) { m_verbose = level; }

private:
   Int_t m_verbose{0};
};

#endif