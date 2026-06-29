#include "PMTModel.h"
#include "PMT.h"
#include "PMTStyle.h"
#include "PMType.h"
#include "SPEFitter.h"

#include "TApplication.h"
#include "TH1D.h"
#include "TError.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TROOT.h"

Int_t example7()
{
   Info("example7.C", "\n\tThe macro starts ...\n");
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
   SPEFitter fitter;
   PMTModel *mod = fitter.CreatePMTModel(hSpec, PMType::Model::TRUNCGAUSS, Q0, s0);

   // The model has been created and seeded.
   // More parameter tuning can be done here

   // Load the model into a TF1
   TF1 *fit = new TF1("dfTF1", mod, xmin, xmax, mod->NPar());
   for (UInt_t ipar = 0; ipar < mod->NPar(); ++ipar) {
      ROOT::Fit::ParameterSettings par = mod->ParSettings(ipar);
      fit->SetParName(ipar, par.Name().c_str());
      fit->SetParameter(ipar, par.Value());
      fit->SetParLimits(ipar, par.LowerLimit(), par.UpperLimit());
   }
   fit->SetLineColor(kBlue);
   hSpec->Fit(fit, "R");

   // Stop the clock after fiting
   sw.Stop();
   Info("example7.C", "Fit took %.0fms Real time and %.0fms CPU time", sw.RealTime() * 1000, sw.CpuTime() * 1000);

   // Calculate the gain and print the result
   Double_t Qfit = mod->ParSettings(4).Value();
   Double_t alphaFit = mod->ParSettings(6).Value();
   Double_t wFit = mod->ParSettings(7).Value();
   Double_t Gfit = wFit / alphaFit + (1.0 - wFit) * Qfit;
   Info("example7.C", "\tTrue Gain : %.2f\n\tBF Gain   : %.2f\n\tDeviation : %.2f%%", Gtrue, Gfit,
        (Gfit / Gtrue - 1.0) * 100.0);

   // Display the result
   ROOT::DisableImplicitMT();
   TCanvas *c1 = new TCanvas("c1", "");
   c1->cd();
   c1->SetLogy();
   hSpec->Draw("PEZ");
   fit->Draw("SAME L");
   c1->Update();
   c1->WaitPrimitive();
   return 0;
}

int main()
{
   TApplication theApp("App", 0, 0);

   Int_t status = example7();

   return status;
}
