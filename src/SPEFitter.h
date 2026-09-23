#ifndef SPEFITTER_H
#define SPEFITTER_H

#include "PMType.h"

#include "RtypesCore.h"
#include "TFitResultPtr.h"

#include <map>
#include <string>

class IModel;
class DFTmethod;
class PMTModel;
class NumIntegration;
class TH1;
class TF1;

class SPEFitter {
public:
   SPEFitter() = default;
   ~SPEFitter() = default;

   std::map<std::string, Double_t>
   GenerateSeeds(TH1 *hspec, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   NumIntegration *
   MakeNumethod(TH1 *hspec, PMType::Response sper, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   DFTmethod *
   MakeDFTmethod(TH1 *hspec, PMType::Response sper, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   PMTModel *
   MakePMTModel(TH1 *hspec, PMType::Model model, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   IModel *MakeModel(const Int_t nBins, const Double_t wBin, const PMType::Method method,
                     const PMType::Model model = PMType::SIMPLEGAUSS, const PMType::Response response = PMType::GAUSS);
   IModel *MakeModel(const TH1 *hspec, const PMType::Method method, const PMType::Model model = PMType::SIMPLEGAUSS,
                     const PMType::Response response = PMType::GAUSS);
   TF1 *MakeTF1(const TH1 *hspec, const PMType::Method method, const PMType::Model model = PMType::SIMPLEGAUSS,
                const PMType::Response response = PMType::GAUSS);
   TF1 *MakeTF1(std::shared_ptr<IModel> model);
   TF1 *MakeTF1(IModel *model);
   std::map<std::string, Double_t>
   SeedModel(IModel *model, TH1 *hspec, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   std::map<std::string, Double_t>
   SeedModel(TF1 *model, TH1 *hspec, const Double_t Q0, const Double_t s0, const Double_t errTol = 0.1);
   Bool_t SeedModel(IModel *model, std::map<std::string, Double_t> seeds);
   Bool_t SeedModel(TF1 *model, const std::map<std::string, Double_t> seeds);
   TFitResultPtr HybridMinimize(IModel *model, TH1 *hspec, Int_t maxItersGA = 150, Int_t maxItersSimplex = 1e5,
                                Int_t maxItersMigrad = 1e8, Double_t tolSimplex = 1.0, Double_t tolMigrad = 1e-2);
   void SetVerbose(Int_t level) { m_verbose = level; }

private:
   Int_t m_verbose{0};
};

#endif