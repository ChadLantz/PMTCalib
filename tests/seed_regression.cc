#include "DFTmethod.h" // IWYU pragma: keep
#include "NumIntegration.h" // IWYU pragma: keep
#include "PMT.h"
#include "PMTModel.h" // IWYU pragma: keep
#include "SPEFitter.h"

#include "TH1D.h" // IWYU pragma: keep

#include <cassert>
#include <cmath>
#include <map>
#include <string>

namespace {

void CheckSeeds(const std::map<std::string, Double_t> &seeds)
{
   for (const auto &[name, value] : seeds)
      assert(std::isfinite(value) && "seed must be finite");

   assert(seeds.at("Norm") > 0.0);
   assert(seeds.at("#sigma_{0}") > 0.0);
   assert(seeds.at("#mu") > 0.0);
   assert(seeds.at("Q") > seeds.at("Q0"));
   assert(seeds.at("#sigma") > 0.0);
   assert(seeds.at("#alpha") > 0.0);
   assert(seeds.at("w") > 0.0 && seeds.at("w") < 1.0);
   assert(seeds.at("xMin") < seeds.at("xMax"));
}

template <typename Model> void CheckModel(Model *model)
{
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar) {
      const auto &settings = model->ParSettings(ipar);
      assert(std::isfinite(settings.Value()));
      assert(std::isfinite(settings.LowerLimit()));
      assert(std::isfinite(settings.UpperLimit()));
      assert(settings.LowerLimit() < settings.UpperLimit());
      assert(settings.Value() >= settings.LowerLimit());
      assert(settings.Value() <= settings.UpperLimit());
   }
}

TH1D *MakeSparseSpectrum()
{
   auto *spectrum = new TH1D("sparse", "sparse", 100, -20.0, 180.0);
   spectrum->Fill(0.0, 3.0);
   spectrum->Fill(40.0, 2.0);
   spectrum->Fill(80.0, 1.0);
   return spectrum;
}

TH1D *MakeFlatTailSpectrum()
{
   auto *spectrum = new TH1D("flatTail", "flatTail", 100, -20.0, 180.0);
   for (Int_t bin = 1; bin <= spectrum->GetNbinsX(); ++bin)
      spectrum->SetBinContent(bin, bin < 15 ? 10.0 : 1.0);
   return spectrum;
}

TH1D *MakeShiftedSparseSpectrum()
{
   auto *spectrum = new TH1D("shiftedSparse", "shiftedSparse", 100, 0.0, 200.0);
   spectrum->Fill(35.0, 12.0);
   spectrum->Fill(75.0, 4.0);
   spectrum->Fill(120.0, 2.0);
   return spectrum;
}

TH1D *MakeSingleSpikeSpectrum()
{
   auto *spectrum = new TH1D("singleSpike", "singleSpike", 100, -20.0, 180.0);
   spectrum->Fill(80.0, 25.0);
   return spectrum;
}

void CheckSpectrum(TH1 *spectrum, Double_t Q0, Double_t s0)
{
   SPEFitter fitter;
   const auto seeds = fitter.GenerateSeeds(spectrum, Q0, s0);
   CheckSeeds(seeds);

   CheckModel(fitter.CreatePMTModel(spectrum, PMType::Model::SIMPLEGAUSS, Q0, s0));
   CheckModel(fitter.CreateDFTmethod(spectrum, PMType::Response::GAUSS, Q0, s0));
   CheckModel(fitter.CreateNumethod(spectrum, PMType::Response::GAUSS, Q0, s0));
}

} // namespace

int main()
{
   constexpr Int_t nbins = 200;
   constexpr Double_t xmin = -50.0;
   constexpr Double_t xmax = 450.0;
   constexpr Double_t s0 = 2.0;
   constexpr Double_t Q = 40.0;
   constexpr Double_t sigma = 13.0;
   constexpr Double_t alpha = 1.0 / 8.0;
   constexpr Double_t w = 0.2;
   Double_t responseParameters[4] = {Q, sigma, alpha, w};

   PMT lowLight(nbins, xmin, xmax, 20.0, s0, PMType::Response::GAUSS, responseParameters);
   CheckSpectrum(lowLight.GenSpectrum(20000, 0.25), 20.0, s0);

   PMT nominal(nbins, xmin, xmax, 0.0, s0, PMType::Response::GAUSS, responseParameters);
   CheckSpectrum(nominal.GenSpectrum(20000, 1.0), 0.0, s0);

   CheckSpectrum(MakeSparseSpectrum(), 0.0, s0);
   CheckSpectrum(MakeFlatTailSpectrum(), 0.0, s0);
   CheckSpectrum(MakeShiftedSparseSpectrum(), 35.0, s0);
   CheckSpectrum(MakeSingleSpikeSpectrum(), 0.0, s0);
   return 0;
}