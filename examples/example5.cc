#include "NumIntegration.h"
#include "PMT.h"
#include "PMTStyle.h"
#include "PMType.h"
#include "SPEFitter.h"

#include "TApplication.h"
#include "TH1D.h"
#include "TError.h"
#include "TCanvas.h"
#include "TStopwatch.h"

Int_t example5()
{
   Info("example5.C", "\n\tThe macro starts ...\n");
   gROOT->Reset();
   PMTStyle::SetDefaultStyle();
   ROOT::EnableImplicitMT();

   // Spectrum generation parameters
   Int_t nbins = 200;
   Double_t xmin = 0.0;
   Double_t xmax = 800.0;
   Double_t Norm = 2.0e+5;
   Double_t Q0 = 20.0;
   Double_t s0 = 2.0;
   Double_t mu = 0.25;
   Double_t Q = 40.0;
   Double_t s = 13.0;
   Double_t alpha = 1.0 / 8.0;
   Double_t w = 0.2;
   Double_t Gtrue = w * alpha + (1.0 - w) * Q;
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
   NumIntegration *num = fitter.CreateNumethod(hSpec, PMType::Response::GAUSS, Q0, s0);

   // The model has been created and seeded.
   // More parameter tuning can be done here

   // Minimize
   fitter.HybridMinimize(num, hSpec);

   // Stop the clock after fiting
   sw.Stop();
   Info("example5.C", "Fit took %.0fms Real time and %.0fms CPU time", sw.RealTime() * 1000, sw.CpuTime() * 1000);

   // Calculate the gain and print the result
   Double_t Qfit = num->ParSettings(4).Value();
   Double_t alphaFit = num->ParSettings(6).Value();
   Double_t wFit = num->ParSettings(7).Value();
   Double_t Gfit = wFit / alphaFit + (1.0 - wFit) * Qfit;
   Info("example5.C", "\tTrue Gain : %.2f\n\tBF Gain   : %.2f\n\tDeviation : %.2f%%", Gtrue, Gfit,
        (Gfit / Gtrue - 1.0) * 100.0);

   TCanvas *c1 = new TCanvas("c1", "");
   c1->cd();
   c1->SetLogy();
   hSpec->Draw("PEZ");
   TGraph *grBF = num->GetGraph();
   grBF->Draw("SAME,L");

   c1->Update();
   c1->WaitPrimitive();

   return 0;
}

int main()
{
   TApplication theApp("App", 0, 0);

   Int_t status = example5();

   return status;
}
