#include "SPEFitter.h"
#include "IModel.h"
#include "NumIntegration.h"
#include "DFTmethod.h"
#include "PMTModel.h"

#include "Fit/BinData.h"
#include "Fit/DataOptions.h"
#include "Fit/DataRange.h"
#include "Fit/FitConfig.h"
#include "Fit/FitResult.h"
#include "Fit/Fitter.h"
#include "Fit/ParameterSettings.h"
#include "HFitInterface.h"
#include "Math/Minimizer.h"
#include "TH1.h"
#include "TF1.h"
#include <TMath.h>
#include <RtypesCore.h>

////////////////////////////////////////////////////////////////////////////////
/// Generate seeds used by all spectra fit by the PMTCalib library given the
/// histogram to be fit amd the same PMT's dark current fit. Initially seeds
/// pedestal location and width from the dark current fit, then re-fits the
/// pedestal, re-seeding with the result. The mean photon count is estimated
/// by assuming P(0) = e^{-lambda} (Poisson distributed).
///
/// @param histo Histogram containing the low-light PMT spectrum
/// @param dcFit Gaussian fit to the dark current spectrum of the given PMT
std::map<std::string, Double_t> SPEFitter::GenerateSeeds(TH1 *hspec, const Double_t Q0, const Double_t s0)
{
   std::map<std::string, Double_t> seeds;
   Double_t wbin = hspec->GetBinWidth(1);
   seeds["Q0"] = Q0;
   seeds["#sigma_{0}"] = s0;
   seeds["Q"] = hspec->GetMean();
   // Norm means different things depending on the method. DFT is integral normalized, PMTModel is Q0 peak normalized
   seeds["pedHeight"] = hspec->GetBinContent(hspec->FindBin(Q0));
   seeds["Norm"] = hspec->Integral();
   seeds["xMin"] = hspec->GetBinCenter(hspec->FindFirstBinAbove(2)); //< Assumes raw counts for now
   seeds["xMax"] = hspec->GetBinCenter(hspec->FindLastBinAbove(2));  //< ^^ same

   // If we have the resolution, fit the pedestal
   if (hspec->GetBinWidth(1) < 2.0 * s0) {
      // Fit the pedestal informed by the dark current fit
      TF1 pedFit("ped", "gausn", Q0 - 2.0 * s0, Q0 + 2.0 * s0);
      pedFit.SetParameters(0.5 * wbin * (seeds.at("Norm") + seeds.at("pedHeight") ), Q0, s0);
      pedFit.SetParLimits(0, wbin * seeds.at("pedHeight"), wbin * seeds.at("Norm"));
      pedFit.SetParLimits(1, Q0 - s0, Q0 + s0);
      pedFit.SetParLimits(2, 0.75 * s0 , 1.5 * s0);
      hspec->Fit(&pedFit, "RQ");

      // Adjust Q0 and sigma0 for this spectrum
      seeds["pedPop"] = pedFit.GetParameter(0) / wbin;
      seeds["Q0"] = pedFit.GetParameter(1);
      seeds["#sigma_{0}"] = pedFit.GetParameter(2);

      // Provide an estimate for mu
      seeds["#mu"] = TMath::Log(wbin * seeds.at("Norm") / pedFit.GetParameter(0));
   } else {
      Warning("GenerateSeeds", "binWidth is > 2 sigma0, cannot pre-fit pedestal");
      // Just use the pedestal bin count
      seeds["pedPop"] = seeds.at("pedHeight");
      seeds["#mu"] = TMath::Log(seeds.at("Norm") / seeds.at("pedHeight"));
   }

   // Calculate gain and lambda
   seeds["gain"] = (seeds.at("Q") - seeds.at("Q0")) / seeds.at("#mu");
   seeds["#lambda"] = 1.0 / seeds["Q"];

   // Fit the tail to seed alpha
   Double_t tailStart = std::max(3.0 * seeds.at("gain"), seeds.at("Q"));
   Double_t lnRise = TMath::Log(0.5 * hspec->GetBinContent(hspec->FindBin(tailStart)));
   Double_t run = seeds.at("xMax") - tailStart;
   seeds["#alpha"] = lnRise / run;

   seeds["w"] = 0.2;      // Fraction of PEs that miss the first dynode
   seeds["#theta"] = 7.0; // Only used in test function
   seeds["#kappa"] = 7.0; // Peak prominance. Not sure how to seed yet

   if (m_verbose > 0) {
      Info("SPEFitter::GenerateSeeds", "Initial seeds:");
      for (auto [key, val] : seeds) {
         Info("SPEFitter::GenerateSeeds", "\t%s\t: %.2e", key.c_str(), val);
      }
   }

   return seeds;
}

////////////////////////////////////////////////////////////////////////////////
/// Three step minimization algorithm. Start with Genetic up to maxCallsGA, then
/// pass the result to Minuit2 Simplex up to maxCallsLocal, before finally
/// passing the result to Minuit2 Migrad for the final fit
///
/// @param model Templated fit model. Can be DFTmethod, NumIntegration, or PMTModel
/// @param hspec Histogram containing PMT spectrum to be fit
/// @param maxCallsGA Maximum number of calls for the Genetic algorithm
/// @param maxCallsLocal Maximum number of calls for Minuit2 (Simplex and Migrad)
/// @param tolGA Genetic algorithm minimization tolerance
/// @param tolSimplex Simplex algorithm minimization tolerance
/// @param tolMigrad Migrad algorithm minimization tolerance
/// @return Fit result pointer containing the results of the Migrad minimization
ROOT::Fit::FitResult SPEFitter::HybridMinimize(IModel *model, TH1 *hspec, Int_t maxItersGA, Int_t maxItersSimplex,
                                               Int_t maxItersMigrad, Double_t tolSimplex, Double_t tolMigrad)
{
   // Construct the BinData object with options and range
   Double_t xmin = 0.0, xmax = 0.0;
   model->GetRange(xmin, xmax);
   ROOT::Fit::DataRange range(xmin, xmax);
   ROOT::Fit::DataOptions opt;
   opt.fUseRange = kTRUE;
   ROOT::Fit::BinData data(opt, range);
   ROOT::Fit::FillData(data, hspec);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(m_verbose);
   fitter.SetFunction(*model);
   fitter.Config().SetParamsSettings(model->ParamsSettings());
   fitter.Config().SetMinimizer("Genetic");
   fitter.Config().MinimizerOptions().SetMaxIterations(maxItersGA);
   fitter.Fit(data);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Simplex");
   fitter.Config().MinimizerOptions().SetMaxIterations(maxItersSimplex);
   fitter.Config().MinimizerOptions().SetTolerance(tolSimplex);
   fitter.Fit(data);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Migrad");
   fitter.Config().MinimizerOptions().SetMaxIterations(maxItersMigrad);
   fitter.Config().MinimizerOptions().SetTolerance(tolMigrad);
   fitter.Fit(data);

   // Set the model parameters with the fitted result
   model->SetFitResult(fitter.Result());

   return fitter.Result();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a NumIntegration object with initial parameters and limits set by
/// examining the spectrum to be fit
///
/// @param hspec Histogram containing PMT spectrum to be fit
/// @param sper Single Photo-Electron Response function (enum) to be used
/// @param Q0 Pedestal mean (usually from dark current data)
/// @param s0 Pedestal width (usually from dark current data)
/// @return ROOT::Fit::FitResult containing fit results
NumIntegration *SPEFitter::CreateNumethod(TH1 *hspec, PMType::Response sper, Double_t Q0, Double_t s0)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0);
   std::map<std::string, std::pair<Double_t, Double_t>> limits;

   limits["Norm"] = {0.75 * seeds.at("Norm"), 1.25 * seeds.at("Norm") };
   limits["Q0"] = {seeds.at("Q0") - seeds.at("#sigma_{0}"), seeds.at("Q0") + seeds.at("#sigma_{0}")};
   limits["#sigma_{0}"] = {0.5 * seeds.at("#sigma_{0}"), 3.0 * seeds.at("#sigma_{0}")};
   limits["#mu"] = {0.02, 2.0 * seeds.at("#mu")};
   limits["#lambda"] = {0, 1.0};
   limits["#theta"] = {0.2 * seeds.at("#theta"), 5.0 * seeds.at("#theta")};
   limits["#kappa"] = {0.2 * seeds.at("#kappa"), 5.0 * seeds.at("#kappa")};
   limits["Q"] = {0.2 * seeds.at("Q"), 5.0 * seeds.at("Q")};
   limits["#sigma"] = {seeds.at("#sigma_{0}"), 15.0 * seeds.at("#sigma_{0}")};
   limits["#alpha"] = limits["#alpha_{1}"] =
      limits["#alpha_{2}"] = {0.1 * seeds.at("#alpha"), 5.0 * seeds.at("#alpha")};
   limits["w"] = limits["w_{1}"] = limits["w_{2}"] = {0.01, 0.75};

   NumIntegration *num =
      new NumIntegration(hspec->GetNbinsX(), hspec->GetBinWidth(1), seeds.at("xMin"), seeds.at("xMax"), sper);
   for (UInt_t ipar = 0; ipar < num->NPar(); ++ipar) {
      std::string parName = num->ParSettings(ipar).Name();
      num->ParSettings(ipar).SetValue(seeds.at(parName));
      num->ParSettings(ipar).SetLimits(limits.at(parName).first, limits.at(parName).second);
   }

   return num;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a DFTmethod object with initial parameters and limits set by examining
/// the spectrum to be fit
///
/// @param hspec Histogram containing PMT spectrum to be fit
/// @param sper Single Photo-Electron Response function (enum) to be used
/// @param Q0 Pedestal mean (usually from dark current data)
/// @param s0 Pedestal width (usually from dark current data)
/// @return ROOT::Fit::FitResult containing fit results
DFTmethod *SPEFitter::CreateDFTmethod(TH1 *hspec, PMType::Response sper, Double_t Q0, Double_t s0)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0);
   std::map<std::string, std::pair<Double_t, Double_t>> limits;

   limits["Norm"] = {0.9 * seeds.at("Norm"), 1.1 * seeds.at("Norm")};
   limits["Q0"] = {seeds.at("Q0") - seeds.at("#sigma_{0}"), seeds.at("Q0") + seeds.at("#sigma_{0}")};
   limits["#sigma_{0}"] = {0.0, 3.0 * seeds.at("#sigma_{0}")};
   limits["#mu"] = {0.02, 2.0 * seeds.at("#mu")};
   limits["#lambda"] = {0.0, 1.0};
   limits["#kappa"] = {0.0, 5.0 * seeds.at("#kappa")};
   limits["#theta"] = {0.0, 5.0 * seeds.at("#theta")};
   limits["Q"] = {0.8 * seeds.at("gain"), 1.2 * seeds.at("gain")};
   limits["#lambda"] = {0.5 / seeds.at("Q"), 3.0 / seeds.at("Q")};
   limits["#sigma"] = {seeds.at("#sigma_{0}"), 10.0 * seeds.at("#sigma_{0}")};
   limits["#alpha"] = limits["#alpha_{1}"] =
      limits["#alpha_{2}"] = {0.1 * seeds.at("#alpha"), 5.0 * seeds.at("#alpha")};
   limits["w"] = limits["w_{1}"] = limits["w_{2}"] = {0.01, 0.75};
   seeds["#sigma"] = 3.0 * seeds.at("#sigma_{0}");

   DFTmethod *dft = new DFTmethod(hspec->GetNbinsX(), hspec->GetBinWidth(1), seeds.at("xMin"), seeds.at("xMax"), sper);
   for (UInt_t ipar = 0; ipar < dft->NPar(); ++ipar) {
      std::string parName = dft->ParSettings(ipar).Name();
      dft->ParSettings(ipar).SetValue(seeds.at(parName));
      dft->ParSettings(ipar).SetLimits(limits.at(parName).first, limits.at(parName).second);
   }
   return dft;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a PMTModel object with initial parameters and limits set by examining
/// the spectrum to be fit
///
/// @param hspec Histogram containing PMT spectrum to be fit
/// @param model PMT model function (enum) to be used
/// @param Q0 Pedestal mean (usually from dark current data)
/// @param s0 Pedestal width (usually from dark current data)
/// @return ROOT::Fit::FitResult containing fit results
PMTModel *SPEFitter::CreatePMTModel(TH1 *hspec, PMType::Model model, Double_t Q0, Double_t s0)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0);
   PMTModel *pmt = new PMTModel(hspec->GetNbinsX(), hspec->GetBinWidth(1), seeds.at("xMin"), seeds.at("xMax"), model);

   pmt->ParSettings(0).SetValue(0.1 * seeds.at("Norm")); 
   pmt->ParSettings(0).SetLimits(0.05 * seeds.at("Norm"), 0.15 * seeds.at("Norm"));

   pmt->ParSettings(1).SetValue(seeds.at("Q0"));
   pmt->ParSettings(1).SetLimits(seeds.at("Q0") - seeds.at("#sigma_{0}"), seeds.at("Q0") + seeds.at("#sigma_{0}"));

   pmt->ParSettings(2).SetValue(seeds.at("#sigma_{0}"));
   pmt->ParSettings(2).SetLimits(0.5 * seeds.at("#sigma_{0}"), 1.5 * seeds.at("#sigma_{0}"));

   pmt->ParSettings(3).SetValue(seeds.at("#mu"));
   pmt->ParSettings(3).SetLimits(0.02, std::max(1.5 * seeds.at("#mu"), 1.0));

   pmt->ParSettings(4).SetValue(seeds.at("gain"));
   pmt->ParSettings(4).SetLimits(0.5 * seeds.at("gain"), 1.5 * seeds.at("gain"));

   pmt->ParSettings(5).SetValue(2.0 * seeds.at("#sigma_{0}"));
   pmt->ParSettings(5).SetLimits(seeds.at("#sigma_{0}"), 5.0 * seeds.at("#sigma_{0}"));

   pmt->ParSettings(6).SetValue(seeds.at("#alpha"));
   pmt->ParSettings(6).SetLimits(0.01 * seeds.at("#alpha"), 3.0 * seeds.at("#alpha"));

   pmt->ParSettings(7).SetValue(seeds.at("w"));
   pmt->ParSettings(7).SetLimits(0.00, 0.65);

   return pmt;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TF1 containing the given model and transfer the parameter settings
TF1 *SPEFitter::MakeTF1(IModel *model){
   Double_t xMin = 0.0, xMax = 0.0;
   model->GetRange(xMin, xMax);
   TF1 *fit = new TF1("model", model, xMin, xMax, model->NPar());
   for(UInt_t ipar = 0; ipar < model->NPar(); ++ipar){
      ROOT::Fit::ParameterSettings par = model->ParSettings(ipar);
      fit->SetParameter(ipar, par.Value());
      fit->SetParLimits(ipar, par.LowerLimit(), par.UpperLimit());
      fit->SetParName(ipar, par.Name().c_str());
   }
   return fit;
}