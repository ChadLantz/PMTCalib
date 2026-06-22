
#include "DFTmethod.h"
#include "TError.h" // IWYU pragma: keep

#include <RtypesCore.h>
#include <fftw3.h>

struct FFTWState {
   std::vector<Double_t> xvalues;
   std::vector<Double_t> yvalues;
   std::vector<Double_t> wfinBG;
   std::vector<Double_t> wfinSG;
   std::vector<Double_t> fftout;
   fftw_complex *wfoutBG;
   fftw_complex *wfoutSG;
   fftw_complex *wfout;
   fftw_plan FWfftBG;
   fftw_plan FWfftSG;
   fftw_plan BWfft;

   FFTWState(){};

   void init(size_t N, Double_t xmin = 0.0, Double_t step = 0.0)
   {
      xvalues.resize(N);
      for (UInt_t i = 0; i < N; ++i) {
         xvalues[i] = xmin + 1.0 * i * step;
      }
      yvalues.resize(N);
      wfinBG.resize(N);
      wfinSG.resize(N);
      fftout.resize(N);
      wfoutBG = fftw_alloc_complex(N / 2 + 1);
      wfoutSG = fftw_alloc_complex(N / 2 + 1);
      wfout = fftw_alloc_complex(N / 2 + 1);
      FWfftBG = fftw_plan_dft_r2c_1d(N, wfinBG.data(), wfoutBG, FFTW_ESTIMATE);
      FWfftSG = fftw_plan_dft_r2c_1d(N, wfinSG.data(), wfoutSG, FFTW_ESTIMATE);
      BWfft = fftw_plan_dft_c2r_1d(N, wfout, fftout.data(), FFTW_ESTIMATE);
   }

   ~FFTWState()
   {
      if (FWfftBG)
         fftw_destroy_plan(FWfftBG);
      if (FWfftSG)
         fftw_destroy_plan(FWfftSG);
      if (BWfft)
         fftw_destroy_plan(BWfft);
   }
};

thread_local FFTWState g_fftwState;

ClassImp(DFTmethod)

   DFTmethod::DFTmethod()
   : gr(nullptr)
{
}

DFTmethod::DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef)
   : N(2 * _nbins + 60),
     M(N / 2 + 1),
     nbins(_nbins),
     xmin(_xmin),
     xmax(_xmax),
     step((xmax - xmin) / (1.0 * nbins * 1.0)),
     edge(xmin),
     nCalls(0),
     gr(nullptr),
     spef(_spef)
{

}

DFTmethod::~DFTmethod()
{
}

void DFTmethod::CalculateValues()
{
   FFTWState& state = g_fftwState;
   if(!state.FWfftBG){
      state.init(N, xmin, step);
   }

   for (UInt_t i = 0; i < N; i++) {
      Double_t xx = state.xvalues[i] - edge;
      Double_t arg = 0.0;

      if (s0 != 0.0)
         arg = (xx - Q0 + edge) / s0;
      else
         Error("CalculateValues", "Division by zero: sigma0");

      state.wfinBG[i] = 1.0 / (sqrt(2.0 * TMath::Pi()) * s0) * TMath::Exp(-0.5 * arg * arg);
      state.wfinSG[i] = spef.GetValue(xx);
   }

   fftw_execute(state.FWfftBG);
   fftw_execute(state.FWfftSG);

   for (UInt_t i = 0; i < M; i++) {
      Double_t amp_BG = sqrt(pow(state.wfoutBG[i][0], 2.0) + pow(state.wfoutBG[i][1], 2.0));
      Double_t ph_BG = fftPhase(state.wfoutBG[i][1], state.wfoutBG[i][0]);

      double ph = (ph_BG + mu * state.wfoutSG[i][1] * step);
      state.wfout[i][0] = amp_BG * TMath::Exp(mu * state.wfoutSG[i][0] * step) * TMath::Cos(ph);
      state.wfout[i][1] = amp_BG * TMath::Exp(mu * state.wfoutSG[i][0] * step) * TMath::Sin(ph);
   }

   fftw_execute(state.BWfft);

   for (UInt_t i = 0; i < N; i++) {
      state.yvalues[i] = Norm * wbin * TMath::Exp(-mu) * state.fftout[i] / Double_t(N);
   }

   if (gr) {
      for (Int_t i = 0; i < nbins; i++)
         gr->SetPoint(i, state.xvalues[i], state.yvalues[i]);
   } else {
      gr = new TGraph(nbins, state.xvalues.data(), state.yvalues.data());
   }

   return;
}

Double_t DFTmethod::Eval(Double_t *xx, Double_t *pars)
{
   Bool_t parsChanged = kFALSE;
   for (UInt_t par = 0; par < nPars; ++par) {
      if (pars[par] != parCache[par]) {
         parCache[par] = pars[par];
         parsChanged = kTRUE;
      }
   }

   if (parsChanged) {
      nCalls++;
      Norm = pars[0];
      Q0 = pars[1];
      s0 = pars[2];
      mu = pars[3];
      spef.SetParams(&pars[4]);
      CalculateValues();
   }

   return gr->Eval(xx[0]);
}

TGraph *DFTmethod::GetGraph()
{
   CalculateValues();

   TGraph *_gr = static_cast<TGraph *>(gr->Clone());
   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);
   return _gr;
}

TGraph *DFTmethod::GetGraphN(Int_t n)
{
   FFTWState& state = g_fftwState;
   if(!state.FWfftBG){
      state.init(N, xmin, step);
   }

   CalculateValues();
   Double_t dblN = n;

   for (UInt_t i = 0; i < N; i++) {
      Double_t xx = state.xvalues.at(i) - edge;

      Double_t arg = 0.0;
      if (s0 != 0.0)
         arg = (xx - Q0 + edge) / s0;
      else
         Error("GetGraphN", "Division by zero: sigma0");

      state.wfinBG[i] = 1.0 / (sqrt(2.0 * TMath::Pi()) * s0) * TMath::Exp(-0.5 * arg * arg);
      state.wfinSG[i] = spef.GetValue(xx);
   }

   fftw_execute(state.FWfftBG);
   fftw_execute(state.FWfftSG);

   for (UInt_t i = 0; i < M; i++) {
      Double_t amp_BG = sqrt(pow(state.wfoutBG[i][0], 2.0) + pow(state.wfoutBG[i][1], 2.0));
      Double_t ph_BG = fftPhase(state.wfoutBG[i][1], state.wfoutBG[i][0]);

      Double_t amp_SG = sqrt(pow(state.wfoutSG[i][0], 2.0) + pow(state.wfoutSG[i][1], 2.0));
      Double_t ph_SG = fftPhase(state.wfoutSG[i][1], state.wfoutSG[i][0]);

      double ph = (ph_BG + dblN * ph_SG);

      state.wfout[i][0] = amp_BG * pow(step * amp_SG, dblN) * TMath::Cos(ph);
      state.wfout[i][1] = amp_BG * pow(step * amp_SG, dblN) * TMath::Sin(ph);
   }

   fftw_execute(state.BWfft);
   fftw_destroy_plan(state.BWfft);

   Double_t x[nbins];
   Double_t y[nbins];

   for (Int_t i = 0; i < nbins; i++) {
      x[i] = state.xvalues.at(i);
      Double_t y_ = Norm * wbin * TMath::Exp(-mu) / TMath::Factorial(n) * pow(mu, dblN) * state.fftout[i] / Double_t(N);

      if (y_ < 1.0e-10)
         y[i] = 1.e-4; // I don't like this. I'd rather use max to set a floor. This might hide issues
      else
         y[i] = y_;
   }

   TGraph *_gr = new TGraph(nbins, x, y);

   _gr->SetLineWidth(2);
   _gr->SetLineStyle(3);
   _gr->SetLineColor(kBlack);
   _gr->SetMarkerColor(kBlack);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

// Atan2? Try replacing once things are stable
Double_t DFTmethod::fftPhase(Double_t vy, Double_t vz)
{
   Double_t thetayz = -999.0;

   Double_t pi = TMath::Pi();

   if (vz > 0 && vy > 0) {
      Double_t ratio = TMath::Abs(vy / vz);
      thetayz = TMath::ATan(ratio);
   } else if (vz < 0 && vy > 0) {
      Double_t ratio = TMath::Abs(vy / vz);
      thetayz = TMath::ATan(ratio);
      thetayz = pi - thetayz;
   } else if (vz < 0 && vy < 0) {
      Double_t ratio = TMath::Abs(vy / vz);
      thetayz = TMath::ATan(ratio);
      thetayz = thetayz + pi;
   } else if (vz > 0 && vy < 0) {
      Double_t ratio = TMath::Abs(vy / vz);
      thetayz = TMath::ATan(ratio);
      thetayz = 2.0 * pi - thetayz;
   } else if (vz == 0 && vy > 0) {
      thetayz = pi / 2.0;
   } else if (vz == 0 && vy < 0) {
      thetayz = 3.0 * pi / 2.0;
   } else if (vz > 0 && vy == 0) {
      thetayz = 0.0;
   } else if (vz < 0 && vy == 0) {
      thetayz = pi;
   }

   thetayz = fmod(thetayz, 2.0 * pi);

   return thetayz;
}
