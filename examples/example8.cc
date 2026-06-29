#include "PMT.h"
#include "PMTStyle.h"
#include "PMType.h"
#include "PMTModel.h"
#include "SPEFitter.h"

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
#include "TApplication.h"
#include "TH1D.h"
#include "TError.h"
#include "TCanvas.h"
#include "TStopwatch.h"

Int_t example8()
{
   Info("example8.C", "\n\tThe macro starts ...\n");
   gROOT->Reset();
   PMTStyle::SetDefaultStyle();
   ROOT::EnableImplicitMT();

   // Spectrum generation parameters
   Int_t nbins = 250;
   Double_t xmin = -50.0;
   Double_t xmax = 450.0;
   Double_t Norm = 2.0e+5;
   Double_t Q0 = 0.0;
   Double_t s0 = 2.0;
   Double_t mu = 1.2;
   Double_t Q = 40.0;
   Double_t s = 13.0;
   Double_t alpha = 1.0 / 8.0;
   Double_t w = 0.2;
   Double_t Gtrue = w / alpha + (1.0 - w) * Q;
   Double_t p[4] = {Q, s, alpha, w};

   // Generate the spectrum to be fit
   PMT specimen(nbins, xmin, xmax, Q0, s0, PMType::Response::GAUSS, p);
   TH1D *hSpec = specimen.GenSpectrum(Norm, mu);

   // Start the clock before generating the model
   TStopwatch sw;
   sw.Start();

   // Create the fitter/model. The fitter generates seeds from the histogram and
   // the given Pedestal mean and width which is gathered from dark current data
   SPEFitter speFitter;
   speFitter.SetVerbose(2);
   PMTModel *mod = speFitter.CreatePMTModel(hSpec, PMType::Model::TRUNCGAUSS, Q0, s0);

   // The model has been created and seeded.
   // More parameter tuning can be done here

   //Create the BinData object filled from the histogram
   ROOT::Fit::DataRange range(xmin, xmax);
   ROOT::Fit::DataOptions opt;
   opt.fUseRange = kTRUE;
   ROOT::Fit::BinData data(opt, range);
   ROOT::Fit::FillData(data, hSpec);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(1);
   fitter.SetFunction(*mod);
   fitter.Config().SetParamsSettings(mod->ParamsSettings());
   fitter.Config().SetMinimizer("Genetic");
   fitter.Config().MinimizerOptions().SetMaxIterations(30);
   fitter.Fit(data);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < mod->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Simplex");
   fitter.Config().MinimizerOptions().SetMaxIterations(1e6);
   fitter.Config().MinimizerOptions().SetTolerance(1e-1);
   fitter.Fit(data);
   // Feed the parameters back into the fitter
   for (UInt_t ipar = 0; ipar < mod->NPar(); ++ipar)
      fitter.Config().ParSettings(ipar).SetValue(fitter.Result().Parameter(ipar));

   fitter.Config().SetMinimizer("Minuit2", "Migrad");
   fitter.Config().MinimizerOptions().SetMaxIterations(1e6);
   fitter.Config().MinimizerOptions().SetTolerance(1e-2);
   fitter.Fit(data);

   // Stop the clock after fiting
   sw.Stop();
   Info("example8.C", "Fit took %.0fms Real time and %.0fms CPU time", sw.RealTime() * 1000, sw.CpuTime() * 1000);

   // Calculate the gain and print the result
   ROOT::Fit::FitResult result = fitter.Result();
   Double_t Qfit = result.Parameter(4);
   Double_t alphaFit = result.Parameter(6);
   Double_t wFit = result.Parameter(7);
   Double_t Gfit = wFit / alphaFit + (1.0 - wFit) * Qfit;
   Info("example8.C", "\tTrue Gain : %.2f\n\tBF Gain   : %.2f\n\tDeviation : %.2f%%", Gtrue, Gfit,
        (Gfit / Gtrue - 1.0) * 100.0);

   // Wrap the fit in a TF1 for display
   TF1 *displayFit = new TF1("display", mod, xmin, xmax, mod->NPar());
   displayFit->SetFitResult(fitter.Result());

   // Display the result
   ROOT::DisableImplicitMT();
   TCanvas *c1 = new TCanvas("c1", "");
   c1->cd();
   c1->SetLogy();
   hSpec->Draw("PEZ");
   displayFit->Draw("SAME L");
   c1->Update();
   c1->WaitPrimitive();
   return 0;
}

int main()
{
   TApplication theApp("App", 0, 0);

   Int_t status = example8();

   return status;
}
