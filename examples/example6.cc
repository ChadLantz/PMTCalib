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
#include <TRandom.h>
#include <Rtypes.h>

Int_t example6()
{
   Info("example6.C", "\n\tThe macro starts ...\n");
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
   Double_t lambda = 1.0 / Q;
   Double_t theta = 8.4;
   Double_t alpha = 1.0 / 8.0;
   Double_t w = 0.2;
   Double_t Gtrue = w / alpha + (1.0 - w) / lambda;
   Double_t p[8] = {Norm, Q0, s0, mu, lambda, theta, alpha, w};

   // Generate the spectrum to be fit
   PMT specimen(nbins, xmin, xmax, Q0, s0, PMType::Response::GAMMA, &p[4]);
   TH1D *hSpec = specimen.GenSpectrum(Norm, mu);

   // Start the clock before generating the model
   TStopwatch sw;
   sw.Start();

   // Create the fitter/model. The fitter generates seeds from the histogram and
   // the given Pedestal mean and width which is gathered from dark current data
   DFTmethod dft(nbins, hSpec->GetBinWidth(1), xmin, xmax, PMType::Response::GAMMA);
   
   // The model has been created
   // More parameter tuning can be done here
   dft.SetParameters(p); // Here we know the parameters, but this is just an example

   // Load the model into a TF1
   TF1 *fit = new TF1("dfTF1", dft, xmin, xmax, dft.NPar());
   fit->SetParNames("Norm", "Q0", "#sigma_{0}", "#mu", "#lambda", "#kappa", "#alpha", "w");
   // Set the parameters on the TF1, not the DFTmethod since TF1 doesn't import
   // parameter settings from IParametricFunctionOneDim for reasons

   // Scramble them up to test how robust the fit is
   fit->SetParameter(0, Norm); // Norm
   // fit->SetParameter(0, gRandom->Poisson(Norm)); // Norm
   fit->SetParLimits(0, 0.0, 3 * Norm);

   fit->SetParameter(1, Q0); // Q 0
   // fit->SetParameter(1, gRandom->Gaus(Q0, s0)); // Q 0
   fit->SetParLimits(1, Q0 - 3 * s0, Q0 + 3 * s0);
   
   fit->SetParameter(2, s0); // sigma 0
   // fit->SetParameter(2, gRandom->Landau(s0)); // sigma 0
   fit->SetParLimits(2, 0.0, 3 * s0);

   fit->SetParameter(3, mu); // mu
   // fit->SetParameter(3, gRandom->PoissonD(mu)); // mu
   fit->SetParLimits(3, 0.02, 3.0 * mu);

   fit->SetParameter(4, lambda); // lambda
   // fit->SetParameter(4, 0.1 / gRandom->Gaus(Q, 0.5 * Q )); // lambda
   fit->SetParLimits(4, 0.0, 0.9);
   
   fit->SetParameter(5, theta); // kappa
   // fit->SetParameter(5, gRandom->Gaus(kappa, 0.25 * kappa)); // kappa
   fit->SetParLimits(5, 0.0, 3.0 * theta);

   fit->SetParameter(6, alpha); // alpha
   // fit->SetParameter(6, 1.0 / gRandom->Gaus(8.0, 2.0)); // alpha
   fit->SetParLimits(6, 0.0, 0.9);

   fit->SetParameter(7, w); // w
   // fit->SetParameter(7, gRandom->Uniform(0.0, 0.7)); // w
   fit->SetParLimits(7, 0.0, 0.7);

   fit->SetLineColor(kAzure + 1);
   hSpec->Fit(fit, "R");

   // Stop the clock after fiting
   sw.Stop();
   Info("example6.C", "Fit took %.0fms Real time and %.0fms CPU time", sw.RealTime() * 1000, sw.CpuTime() * 1000);

   // Calculate the gain and print the result
   Double_t lambdaFit = fit->GetParameter(4);
   Double_t alphaFit = fit->GetParameter(6);
   Double_t wFit = fit->GetParameter(7);
   Double_t Gfit = wFit / alphaFit + (1.0 - wFit) / lambdaFit;
   Info("example6.C", "\tTrue Gain : %.2f\n\tBF Gain   : %.2f\n\tDeviation : %.2f%%", Gtrue, Gfit,
        (Gfit / Gtrue - 1.0) * 100.0);

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

   Int_t status = example6();

   return status;
}
