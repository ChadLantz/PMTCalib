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
#include "TMath.h"
#include "ROOT/EExecutionPolicy.hxx"
#include <algorithm>
#include <limits>
#include <memory>
#include <cmath>
#include "RtypesCore.h"
#include "TFitResult.h"
#include <TFitResultPtr.h>

////////////////////////////////////////////////////////////////////////////////
/// Generate seeds used by all methods in this library given the
/// histogram to be fit amd the same PMT's dark current fit. Initially seeds
/// pedestal location and width from the dark current fit, then re-fits the
/// pedestal, re-seeding with the result. The mean photon count is estimated
/// by assuming P(0) = e^{-lambda} (Poisson distributed).
///
/// @param histo Histogram containing the low-light PMT spectrum
/// @param Q0 Estimated pedistal position
/// @param s0 Esitmated pedistal width
/// @param errTol Maximum bin error for selecting bin range
std::map<std::string, Double_t>
SPEFitter::GenerateSeeds(TH1 *hspec, const Double_t Q0, const Double_t s0, const Double_t errTol)
{
   std::map<std::string, Double_t> seeds;
   seeds["Q0"] = Q0;
   seeds["#sigma_{0}"] = s0;
   seeds["Norm"] = hspec->Integral();
   Double_t wbin = hspec->GetBinWidth(1);
   Double_t mean = hspec->GetMean();
   Int_t pedBin = hspec->FindBin(Q0);
   Int_t firstPedBin = std::max(1, hspec->FindBin(Q0 - s0));
   Int_t lastPedBin = std::min(hspec->GetNbinsX(), hspec->FindBin(Q0 + 0.5 * s0));
   Double_t pedAmp = 0.0;
   for (Int_t bin = firstPedBin; bin <= lastPedBin; ++bin)
      pedAmp = std::max(pedAmp, hspec->GetBinContent(bin));
   if (!std::isfinite(seeds.at("Norm")) || seeds.at("Norm") <= 0.0)
      seeds["Norm"] = 1.0;

   auto validPositive = [](const Double_t value) { return std::isfinite(value) && value > 0.0; };
   auto estimateMu = [&](const Double_t population) {
      if (validPositive(population) && seeds.at("Norm") > population) {
         const Double_t estimate = TMath::Log(seeds.at("Norm") / population);
         if (std::isfinite(estimate) && estimate > 0.0)
            return std::clamp(estimate, 0.02, 10.0);
      }
      return 0.25;
   };
   auto estimateGain = [&](const Double_t pedestal, const Double_t occupancy) {
      const Double_t observedGain = (mean - pedestal) / occupancy;
      if (validPositive(observedGain))
         return observedGain;
      return std::max(2.0 * std::abs(s0), wbin);
   };
   auto fitIsUsable = [](const TF1 *fit) {
      if (!fit->IsValid())
         return false;
      for (Int_t ipar = 0; ipar < fit->GetNpar(); ++ipar) {
         if (!std::isfinite(fit->GetParameter(ipar)))
            return false;
      }
      return true;
   };

   // Provide first estimates
   Double_t pedPop = pedAmp * s0 * TMath::Sqrt(TMath::TwoPi()) / wbin;
   Double_t mu = estimateMu(pedPop);

   // If we have the resolution, fit the pedestal
   if (hspec->GetBinWidth(1) < 2.0 * s0) {
      // Fit with a Gaussian first to estimate mu
      TF1 *gaus = new TF1("gaus", "gaus", Q0 - 2.5 * s0, Q0 + 1.5 * s0);
      gaus->SetParameters(pedAmp, Q0, s0);
      gaus->SetParLimits(0, 0.5 * pedAmp, seeds.at("Norm"));
      gaus->SetParLimits(1, Q0 - 0.5 * s0, Q0 + 0.5 * s0);
      gaus->SetParLimits(2, 0.5 * s0, 2.0 * s0);
      hspec->Fit(gaus, "LRQ");
      if (!fitIsUsable(gaus)) {
         Warning("GenerateSeeds", "Pedestal fit failed for %s; retaining dark-current seeds", hspec->GetName());
         seeds["pedPop"] = pedPop;
         seeds["#mu"] = mu;
         seeds["Q"] = Q0 + estimateGain(Q0, mu);
         seeds["#sigma"] = 2.0 * s0;
      } else {
      // Adjust Q0 and sigma0 for this spectrum
      pedAmp = gaus->GetParameter(0);
      Double_t Q0Fit = gaus->GetParameter(1);
      Double_t s0Fit = gaus->GetParameter(2);
      pedPop = pedAmp * s0Fit * TMath::Sqrt(TMath::TwoPi()) / wbin; // Overestimates population as mu approaches zero
      // Refine the estimate for mu and gain
      mu = estimateMu(pedPop);
      Double_t gain = estimateGain(Q0Fit, mu);
      Double_t Q1 = Q0Fit + gain; // First estimate at 1PE peak location. Usually an overestimate
      Double_t Q1amp = hspec->GetBinContent(hspec->FindBin(Q0Fit + gain));
      if (!validPositive(Q1amp))
         Q1amp = std::max(pedAmp, 1.0);

      if (m_verbose > 1) {
         Info("GenerateSeeds",
              "Fit pedestal of %s. Provided Q0 = %.2e, s0 = %.2e. Gave:\n\tQ0 = %.2e\n\ts0 = %.2e\n\tmu = %.2f\n\tgain "
              "= %.2e\n\tQ1 = %.2e\n\tQ1amp = %.2e",
              hspec->GetName(), Q0, s0, Q0Fit, s0Fit, mu, gain, Q1, Q1amp);
      }

      // Fit with a
      TF1 *dblGaus =
         new TF1("dblgaus", "gaus(0) + gaus(3)", Q0Fit - s0Fit, std::max(Q0Fit + 2.0 * s0Fit, Q1 + s0Fit));
      dblGaus->SetParameter(0, 0.8 * pedAmp); // Overlap will reduce population
      dblGaus->SetParLimits(0, 0.5 * pedAmp, 1.25 * pedAmp);
      dblGaus->SetParameter(1, Q0Fit);
      dblGaus->SetParLimits(1, Q0Fit - 0.25 * s0Fit, Q0Fit + 0.25 * s0Fit);
      dblGaus->SetParameter(2, s0Fit);
      dblGaus->SetParLimits(2, 0.8 * s0Fit, 1.2 * s0Fit);
      dblGaus->SetParameter(3, Q1amp);
      dblGaus->SetParLimits(3, 0.5 * Q1amp, 2.0 * Q1amp);
      dblGaus->SetParameter(4, 0.8 * Q1);
      dblGaus->SetParLimits(4, std::max(Q0Fit + 0.5 * s0Fit, Q1 - 5.0 * s0Fit), Q1 + 2.5 * s0Fit);
      dblGaus->SetParameter(5, 2.0 * s0Fit);
      dblGaus->SetParLimits(5, 1.5 * s0Fit, 10.0 * s0Fit);
      hspec->Fit(dblGaus, "LRQ");

      if (fitIsUsable(dblGaus)) {
         pedAmp = dblGaus->GetParameter(0);
         seeds["Q0"] = Q0Fit = dblGaus->GetParameter(1);
         seeds["#sigma_{0}"] = s0Fit = dblGaus->GetParameter(2);
         seeds["Q"] = Q1 = dblGaus->GetParameter(4);
         seeds["#sigma"] = dblGaus->GetParameter(5);
         seeds["pedPop"] = pedPop = pedAmp * s0Fit * TMath::Sqrt(TMath::TwoPi()) / wbin;
         seeds["#mu"] = estimateMu(pedPop);
      } else {
         Warning("GenerateSeeds", "Double-pedestal fit failed for %s; retaining single-Gaussian seeds",
                 hspec->GetName());
         seeds["pedPop"] = pedPop;
         seeds["#mu"] = mu;
         seeds["Q"] = Q1;
         seeds["#sigma"] = 2.0 * s0Fit;
      }

      if (m_verbose > 1) {
         Info("GenerateSeeds",
              "Fit double Gauss of %s. Provided Q0 = %.2e, s0 = %.2e. Gave:\n\tQ0 = %.2e\n\ts0 = %.2e\n\tQ = "
              "%.2e\n\tsigma = %.2e\n\tmu = %.2f",
              hspec->GetName(), Q0, s0, Q0Fit, s0Fit, Q1, seeds.at("#sigma"), seeds.at("#mu"));
      }
               }

   } else {
      Warning("GenerateSeeds", "binWidth is > 2 sigma0, cannot pre-fit pedestal");
      // Just use the initial estimates
      seeds["pedPop"] = pedPop;
      seeds["#mu"] = mu;
      seeds["Q"] = Q0 + estimateGain(Q0, mu);
      seeds["#sigma"] = 2.0 * s0;
   }

   // Some spectra won't have a pronounced peak and are just a shoulder on the pedestal, so use that initial value
   Double_t pePeakVal = hspec->GetBinContent(hspec->FindBin(seeds.at("Q")));
   for (UInt_t bin = 1; bin <= hspec->GetNbinsX(); ++bin) {
      Double_t content = hspec->GetBinContent(bin);
      Double_t error = hspec->GetBinError(bin);
      Double_t binLowEdge = hspec->GetBinLowEdge(bin);
      Double_t binUpEdge = binLowEdge + wbin;
      // fill xMin and xMax
      if (content > 0.0 && error / content < errTol) {
         if (!seeds.contains("xMin")) {
            seeds["xMin"] = binLowEdge;
         }
         seeds["xMax"] = binUpEdge;
      }
   }
      if (!seeds.contains("xMin") || !seeds.contains("xMax") || !std::isfinite(seeds["xMin"]) ||
         !std::isfinite(seeds["xMax"]) ||
       ((seeds.contains("xMin") && seeds.contains("xMax") && seeds["xMin"] >= seeds["xMax"]))) {
      Error("GenerateSeeds", "Failed to find valid xMin and xMax for histogram %s", hspec->GetName());
      seeds["xMin"] = hspec->GetXaxis()->GetXmin();
      seeds["xMax"] = hspec->GetXaxis()->GetXmax();
   }

   if (seeds.at("xMax") - seeds.at("xMin") < 2.0 * wbin) {
      seeds["xMin"] = hspec->GetXaxis()->GetXmin();
      seeds["xMax"] = hspec->GetXaxis()->GetXmax();
   }

   if (!std::isfinite(seeds["Q0"]))
      seeds["Q0"] = 0.0;
   if (!validPositive(seeds["#sigma_{0}"]))
      seeds["#sigma_{0}"] = validPositive(wbin) ? wbin : 1.0;
   if (!validPositive(seeds["#mu"]))
      seeds["#mu"] = 0.25;
   const Double_t gainFloor = std::max(2.0 * seeds["#sigma_{0}"], wbin);
   if (!validPositive(seeds["Q"]) || seeds["Q"] <= seeds["Q0"])
      seeds["Q"] = std::max(seeds["Q0"] + gainFloor, gainFloor);
   if (!validPositive(seeds["#sigma"]))
      seeds["#sigma"] = 2.0 * seeds["#sigma_{0}"];
   seeds["#sigma"] = std::max(seeds["#sigma"], 0.5 * seeds["#sigma_{0}"]);
   if (!validPositive(seeds["pedPop"]))
      seeds["pedPop"] = seeds["Norm"] * TMath::Exp(-seeds["#mu"]);

   // Use a small weighted neighborhood to identify a resolved SPE peak without
   // making the seed path depend on an optional ROOT spectrum component.
   const Int_t firstPeakBin = std::max(1, hspec->FindBin(seeds["Q0"] + 2.0 * seeds["#sigma_{0}"]));
   const Int_t lastPeakBin = std::min(hspec->GetNbinsX(), hspec->FindBin(seeds["xMax"]));
   Int_t peakBin = -1;
   Bool_t resolvedPeak = kFALSE;
   const Double_t expectedPeak = seeds["Q0"] + seeds["Q"];
   Double_t peakScore = std::numeric_limits<Double_t>::infinity();
   for (Int_t bin = firstPeakBin + 1; bin < lastPeakBin; ++bin) {
      const Double_t left = hspec->GetBinContent(bin - 1);
      const Double_t center = hspec->GetBinContent(bin);
      const Double_t right = hspec->GetBinContent(bin + 1);
      const Double_t smoothed = (left + 2.0 * center + right) / 4.0;
      const Double_t leftSmoothed = (hspec->GetBinContent(bin - 2) + 2.0 * left + center) / 4.0;
      const Double_t rightSmoothed = (center + 2.0 * right + hspec->GetBinContent(bin + 2)) / 4.0;
      const Double_t prominence = smoothed - std::max(leftSmoothed, rightSmoothed);
      if (smoothed <= leftSmoothed || smoothed < rightSmoothed || prominence < std::max(2.0, 0.05 * smoothed))
         continue;

      const Double_t distance = std::abs(hspec->GetBinCenter(bin) - expectedPeak);
      if (distance < peakScore) {
         peakScore = distance;
         peakBin = bin;
      }
   }
   if (peakBin > 0 && peakScore <= std::max(3.0 * seeds["#sigma_{0}"], 0.25 * seeds["Q"])) {
      resolvedPeak = kTRUE;
      seeds["Q"] = hspec->GetBinCenter(peakBin) - seeds["Q0"];
      if (!validPositive(seeds["Q"]))
         seeds["Q"] = hspec->GetBinCenter(peakBin);

      // Estimate the observed SPE width locally when no double-Gaussian
      // refinement supplied a more informative width.
      if (seeds["#sigma"] <= 2.0 * seeds["#sigma_{0}"]) {
         const Int_t firstWidthBin = std::max(1, peakBin - 3);
         const Int_t lastWidthBin = std::min(hspec->GetNbinsX(), peakBin + 3);
         const Double_t baseline = std::min(hspec->GetBinContent(firstWidthBin),
                                            hspec->GetBinContent(lastWidthBin));
         Double_t weightSum = 0.0;
         Double_t weightedMean = 0.0;
         for (Int_t bin = firstWidthBin; bin <= lastWidthBin; ++bin) {
            const Double_t weight = std::max(0.0, hspec->GetBinContent(bin) - baseline);
            weightSum += weight;
            weightedMean += weight * hspec->GetBinCenter(bin);
         }

         if (validPositive(weightSum)) {
            weightedMean /= weightSum;
            Double_t observedVariance = 0.0;
            for (Int_t bin = firstWidthBin; bin <= lastWidthBin; ++bin) {
               const Double_t weight = std::max(0.0, hspec->GetBinContent(bin) - baseline);
               observedVariance += weight * TMath::Sq(hspec->GetBinCenter(bin) - weightedMean);
            }
            observedVariance /= weightSum;
            const Double_t intrinsicVariance = observedVariance - TMath::Sq(seeds["#sigma_{0}"]);
            const Double_t localSigma = intrinsicVariance > 0.0 ? TMath::Sqrt(intrinsicVariance) : 0.0;
            if (validPositive(localSigma) && localSigma >= 0.5 * seeds["#sigma_{0}"] &&
                localSigma <= 10.0 * seeds["#sigma_{0}"])
               seeds["#sigma"] = localSigma;
         }
      }

      if (m_verbose > 1)
         Info("GenerateSeeds", "Resolved SPE peak for %s near Q = %.2e", hspec->GetName(), seeds["Q"]);
   }

   // Calculate the slope of the tail (#alpha)
   pePeakVal = hspec->GetBinContent(hspec->FindBin(seeds.at("Q0") + seeds.at("Q")));
   const Double_t tailContent = hspec->GetBinContent(hspec->FindBin(seeds.at("xMax")));
   const Double_t tailExcess = pePeakVal - tailContent;
   Double_t run = seeds.at("xMax") - hspec->GetBinCenter(hspec->GetMaximumBin());
   Double_t valleyAlpha = 0.0;
   if (resolvedPeak) {
      Int_t valleyBin = -1;
      Double_t valleyValue = std::numeric_limits<Double_t>::infinity();
      for (Int_t bin = firstPeakBin; bin < peakBin; ++bin) {
         const Double_t smoothed =
            (hspec->GetBinContent(bin - 1) + 2.0 * hspec->GetBinContent(bin) + hspec->GetBinContent(bin + 1)) / 4.0;
         if (smoothed < valleyValue) {
            valleyValue = smoothed;
            valleyBin = bin;
         }
      }
      if (valleyBin > 0) {
         const Double_t valleyDistance = hspec->GetBinCenter(valleyBin) - seeds["Q0"];
         if (validPositive(valleyDistance))
            valleyAlpha = 1.0 / valleyDistance;
      }
   }
   const Double_t fallbackAlpha = validPositive(valleyAlpha)
                                    ? valleyAlpha
                                    : 1.0 / std::max(std::abs(seeds.at("Q")), seeds.at("#sigma_{0}"));
   if (validPositive(tailExcess) && validPositive(run)) {
      const Double_t alpha = TMath::Log(tailExcess) / run;
      seeds["#alpha"] = validPositive(alpha) ? alpha : fallbackAlpha;
   } else {
      seeds["#alpha"] = fallbackAlpha;
   }
   if (!validPositive(seeds["#alpha"]))
      seeds["#alpha"] = 1.0 / std::max(seeds["Q"], seeds["#sigma_{0}"]);

   seeds["w"] = 0.2; // Fraction of PEs that miss the first dynode
   seeds["#lambda"] = 1.0 / seeds.at("Q");
   seeds["#theta"] = 7.0; // Only used in test function
   seeds["#kappa"] = 7.0; // Peak prominance. Not sure how to seed yet
   seeds["#alpha_{1}"] = seeds["#alpha_{2}"] = seeds.at("#alpha");
   seeds["w_{1}"] = seeds["w_{2}"] = seeds.at("w");

   if (m_verbose > 0) {
      Info("SPEFitter::GenerateSeeds", "Initial seeds:");
      for (auto [key, val] : seeds) {
         Info("SPEFitter::GenerateSeeds", "\t%-12s\t: %.2e", key.c_str(), val);
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
TFitResultPtr SPEFitter::HybridMinimize(IModel *model, TH1 *hspec, Int_t maxItersGA, Int_t maxItersSimplex,
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
   fitter.LeastSquareFit(data, ROOT::EExecutionPolicy::kMultiThread);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Simplex");
   fitter.Config().MinimizerOptions().SetMaxIterations(maxItersSimplex);
   fitter.Config().MinimizerOptions().SetTolerance(tolSimplex);
   fitter.LeastSquareFit(data, ROOT::EExecutionPolicy::kMultiThread);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Migrad");
   fitter.Config().MinimizerOptions().SetMaxIterations(maxItersMigrad);
   fitter.Config().MinimizerOptions().SetTolerance(tolMigrad);
   fitter.LeastSquareFit(data, ROOT::EExecutionPolicy::kMultiThread);

   // Create a shared pointer to the result (copy)
   TFitResultPtr result(std::make_shared<TFitResult>(fitter.Result()));
   // Set the model parameters with the fitted result
   model->SetFitResult(result);
   return result;
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
NumIntegration *
SPEFitter::CreateNumethod(TH1 *hspec, PMType::Response sper, Double_t Q0, Double_t s0, const Double_t errTol)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0, errTol);
   std::map<std::string, std::pair<Double_t, Double_t>> limits;

   limits["Norm"] = {0.75 * seeds.at("Norm"), 1.25 * seeds.at("Norm")};
   limits["Q0"] = {seeds.at("Q0") - seeds.at("#sigma_{0}"), seeds.at("Q0") + seeds.at("#sigma_{0}")};
   limits["#sigma_{0}"] = {0.5 * seeds.at("#sigma_{0}"), 3.0 * seeds.at("#sigma_{0}")};
   limits["#mu"] = {std::max(0.5 * seeds.at("#mu"), 0.02), std::max(2.0 * seeds.at("#mu"), 1.0)};
   limits["#lambda"] = {0, 1.0};
   limits["#theta"] = {0.2 * seeds.at("#theta"), 5.0 * seeds.at("#theta")};
   limits["#kappa"] = {0.2 * seeds.at("#kappa"), 5.0 * seeds.at("#kappa")};
   limits["Q"] = {0.2 * seeds.at("Q"), 5.0 * seeds.at("Q")};
   limits["#sigma"] = {seeds.at("#sigma_{0}"), 15.0 * seeds.at("#sigma_{0}")};
   limits["#alpha"] = limits["#alpha_{1}"] =
      limits["#alpha_{2}"] = {0.1 * seeds.at("#alpha"), 5.0 * seeds.at("#alpha")};
   limits["w"] = limits["w_{1}"] = limits["w_{2}"] = {0.01, 1.0 - std::numeric_limits<Double_t>::epsilon()};

   Int_t minBin = hspec->GetBin(seeds.at("xMin"));
   Int_t maxBin = hspec->GetBin(seeds.at("xMax"));
   UInt_t nBins = maxBin - minBin;
   NumIntegration *num = new NumIntegration(nBins, hspec->GetBinWidth(minBin), hspec->GetBinLowEdge(minBin),
                                            hspec->GetBinLowEdge(maxBin) + hspec->GetBinWidth(maxBin), sper);
   for (UInt_t ipar = 0; ipar < num->NPar(); ++ipar) {
      std::string parName = num->ParSettings(ipar).Name();
      num->SetParameter(ipar, seeds.at(parName));
      num->SetParLimits(ipar, limits.at(parName).first, limits.at(parName).second);
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
DFTmethod *
SPEFitter::CreateDFTmethod(TH1 *hspec, PMType::Response sper, Double_t Q0, Double_t s0, const Double_t errTol)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0, errTol);
   std::map<std::string, std::pair<Double_t, Double_t>> limits;

   limits["Norm"] = {0.9 * seeds.at("Norm"), 1.1 * seeds.at("Norm")};
   limits["Q0"] = {seeds.at("Q0") - seeds.at("#sigma_{0}"), seeds.at("Q0") + seeds.at("#sigma_{0}")};
   limits["#sigma_{0}"] = {0.0, 3.0 * seeds.at("#sigma_{0}")};
   limits["#mu"] = {std::max(0.5 * seeds.at("#mu"), 0.02), std::max(2.0 * seeds.at("#mu"), 1.0)};
   limits["#lambda"] = {0.0, 1.0};
   limits["#kappa"] = {0.0, 25.0 * seeds.at("#kappa")};
   limits["Q"] = {0.75 * seeds.at("Q"), 1.5 * seeds.at("Q")};
   limits["#lambda"] = limits["#theta"] = {1.0 / limits.at("Q").second, 1.0 / limits.at("Q").first};
   seeds["#sigma"] = std::max(seeds.at("#sigma"), 3.0 * seeds.at("#sigma_{0}"));
   Bool_t hasS1seed = seeds.find("#sigma") != seeds.end();
   Double_t s1 = hasS1seed ? seeds.at("#sigma") : 2.5 * seeds.at("#sigma_{0}");
   limits["#sigma"] = {s0, 1.5 * s1};
   limits["#alpha"] = limits["#alpha_{1}"] =
      limits["#alpha_{2}"] = {0.1 * seeds.at("#alpha"), 5.0 * seeds.at("#alpha")};
   limits["w"] = limits["w_{1}"] = limits["w_{2}"] = {0.01, 1.0 - std::numeric_limits<Double_t>::epsilon()};
   Int_t minBin = hspec->GetBin(seeds.at("xMin"));
   Int_t maxBin = hspec->GetBin(seeds.at("xMax"));
   UInt_t nBins = maxBin - minBin;
   DFTmethod *dft = new DFTmethod(nBins, hspec->GetBinWidth(minBin), hspec->GetBinLowEdge(minBin),
                                  hspec->GetBinLowEdge(maxBin) + hspec->GetBinWidth(maxBin), sper);
   for (UInt_t ipar = 0; ipar < dft->NPar(); ++ipar) {
      std::string parName = dft->ParSettings(ipar).Name();
      dft->SetParameter(ipar, seeds.at(parName));
      dft->SetParLimits(ipar, limits.at(parName).first, limits.at(parName).second);
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
PMTModel *SPEFitter::CreatePMTModel(TH1 *hspec, PMType::Model model, Double_t Q0, Double_t s0, const Double_t errTol)
{
   std::map<std::string, Double_t> seeds = GenerateSeeds(hspec, Q0, s0, errTol);
   Int_t minBin = hspec->GetBin(seeds.at("xMin"));
   Int_t maxBin = hspec->GetBin(seeds.at("xMax"));
   UInt_t nBins = maxBin - minBin;
   PMTModel *pmt = new PMTModel(nBins, hspec->GetBinWidth(minBin), hspec->GetBinLowEdge(minBin),
                                hspec->GetBinLowEdge(maxBin) + hspec->GetBinWidth(maxBin), model);

   pmt->SetParameter(0, seeds.at("Norm"));
   pmt->ParSettings(0).SetLimits(0.9 * seeds.at("Norm"), 1.1 * seeds.at("Norm"));

   pmt->SetParameter(1, seeds.at("Q0"));
   pmt->SetParLimits(1, seeds.at("Q0") - 0.25 * seeds.at("#sigma_{0}"), seeds.at("Q0") + 0.25 * seeds.at("#sigma_{0}"));

   pmt->SetParameter(2, seeds.at("#sigma_{0}"));
   pmt->SetParLimits(2, 0.75 * seeds.at("#sigma_{0}"), 1.5 * seeds.at("#sigma_{0}"));

   pmt->SetParameter(3, seeds.at("#mu"));
   pmt->SetParLimits(3, std::max(0.5 * seeds.at("#mu"), 0.02), std::max(2.0 * seeds.at("#mu"), 1.0));

   pmt->SetParameter(4, seeds.at("Q"));
   pmt->SetParLimits(4, 0.8 * seeds.at("Q"), 1.2 * seeds.at("Q"));

   Bool_t hasS1seed = seeds.find("#sigma_{1}") != seeds.end();
   Double_t s1 = hasS1seed ? seeds.at("#sigma_{1}") : 2.5 * seeds.at("#sigma_{0}");
   pmt->SetParameter(5, s1);
   pmt->SetParLimits(5, hasS1seed ? 0.75 * s1 : seeds.at("#sigma_{0}"),
                     hasS1seed ? 1.5 * s1 : 5.0 * seeds.at("#sigma_{0}"));

   pmt->SetParameter(6, seeds.at("#alpha"));
   pmt->SetParLimits(6, 0.5 * seeds.at("#alpha"), 2.0 * seeds.at("#alpha"));

   pmt->SetParameter(7, seeds.at("w"));
   pmt->SetParLimits(7, 0.01, 1.0 - std::numeric_limits<Double_t>::epsilon());

   return pmt;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TF1 containing the given model and transfer the parameter settings
TF1 *SPEFitter::MakeTF1(IModel *model)
{
   Double_t xMin = 0.0, xMax = 0.0;
   model->GetRange(xMin, xMax);
   TF1 *fit = new TF1("model", model, &IModel::EvalPar, xMin, xMax, model->NPar());
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar) {
      ROOT::Fit::ParameterSettings par = model->ParSettings(ipar);
      fit->SetParameter(ipar, par.Value());
      fit->SetParLimits(ipar, par.LowerLimit(), par.UpperLimit());
      fit->SetParName(ipar, par.Name().c_str());
   }
   fit->SetChisquare(model->GetChiSquare());
   fit->SetNDF(model->GetNDF());
   return fit;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TF1 containing the given model and transfer the parameter settings
TF1 *SPEFitter::MakeTF1(std::shared_ptr<IModel> model)
{
   Double_t xMin = 0.0, xMax = 0.0;
   model->GetRange(xMin, xMax);
   TF1 *fit = new TF1("model", model, &IModel::EvalPar, xMin, xMax, model->NPar());
   for (UInt_t ipar = 0; ipar < model->NPar(); ++ipar) {
      ROOT::Fit::ParameterSettings par = model->ParSettings(ipar);
      fit->SetParameter(ipar, par.Value());
      fit->SetParLimits(ipar, par.LowerLimit(), par.UpperLimit());
      fit->SetParName(ipar, par.Name().c_str());
   }
   fit->SetChisquare(model->GetChiSquare());
   fit->SetNDF(model->GetNDF());
   return fit;
}