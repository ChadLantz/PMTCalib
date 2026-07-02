#include "DFTmethod.h"
#include "PMT.h"
#include "PMTStyle.h"
#include "PMType.h"
#include "SPEFitter.h"

#include "TApplication.h"
#include "TH1D.h"
#include "TError.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <Fit/FitResult.h>
#include <Rtypes.h>

Int_t example4()
{
   Info("example4.C", "\n\tThe macro starts ...\n");
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
   Double_t lambda = 1.0 / 40.0;
   Double_t theta = 8.4;
   Double_t alpha = 1.0 / 8.0;
   Double_t w = 0.2;
   Double_t Gtrue = w / alpha + (1.0 - w) / lambda;
   Double_t p[4] = {lambda, theta, alpha, w};

   // Generate the spectrum to be fit
   PMT specimen(nbins, xmin, xmax, Q0, s0, PMType::Response::GAMMA, p);
   TH1D *hSpec = specimen.GenSpectrum(Norm, mu);

   // Start the clock before generating the model
   TStopwatch sw;
   sw.Start();

   // Create the fitter/model. The fitter generates seeds from the histogram and
   // the given Pedestal mean and width which is gathered from dark current data
   SPEFitter fitter;
   DFTmethod *method = fitter.CreateDFTmethod(hSpec, PMType::Response::GAMMA, Q0, s0);

   // The model has been created and seeded by the SPEFitter for convenience
   // More parameter tuning can be done here

   // Minimize
   ROOT::Fit::FitResult result = fitter.HybridMinimize(method, hSpec);

   // Stop the clock after fiting
   sw.Stop();
   Info("example4.C", "Fit took %.0fms Real time and %.0fms CPU time", sw.RealTime() * 1000, sw.CpuTime() * 1000);

   // Calculate the gain and print the result
   Double_t Gfit = method->Gain(method->Parameters()); // Gotta figure out how to avoid this before committing
   Info("example4.C", "\tTrue Gain : %.2f\n\tBF Gain   : %.2f\n\tDeviation : %.2f%%", Gtrue, Gfit,
        (Gfit / Gtrue - 1.0) * 100.0);

   // Display the result
   TF1 *displayFit = fitter.MakeTF1(method); // Wrap the model in a TF1 with the helper
   displayFit->SetFitResult(result);
   displayFit->SetLineColor(kAzure + 1);
   result.Print(std::cout);
   
   TCanvas *c1 = new TCanvas("c1", "");
   c1->cd();
   c1->SetLogy();
   hSpec->Draw("PEZ");
   displayFit->Draw("SAME L");
   TGraph *grPE[25];
   for (Int_t i = 0; i < 25; i++) {
      grPE[i] = method->GetGraphN(i);
      grPE[i]->Draw("SAME L");
   }
   c1->Update();
   c1->WaitPrimitive();
   return 0;
}

int main()
{
   TApplication theApp("App", 0, 0);

   Int_t status = example4();

   return status;
}
