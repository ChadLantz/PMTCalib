#include "DFTmethod.h"
#include "SPEResponse.h"

#include <RtypesCore.h>
#include "TError.h" // IWYU pragma: keep
#include <TMath.h>

#include <fftw3.h>

struct FFTWState {
   size_t N;
   size_t M;
   std::vector<Double_t> wfinBG; // FFT BG in
   std::vector<Double_t> wfinSG; // FFT SG in
   std::vector<Double_t> fftout; // IFFT real result
   fftw_complex *wfoutBG;
   fftw_complex *wfoutSG;
   fftw_complex *wfout; //
   fftw_plan FWfftBG;   // FFT BG?
   fftw_plan FWfftSG;   // FFT SG?
   fftw_plan BWfft;     //< IFFT
   TGraph *gr;          //< This state owns its graph and we will return clones
   Bool_t INITIALIZED; 

   FFTWState():INITIALIZED(kFALSE) {};

   Bool_t IsInitialized(){ return INITIALIZED; }

   void init(size_t _nBins, Double_t _xmin = 0.0, Double_t _xmax = 0.0)
   {
      N = 2 * _nBins + 60;
      M = N / 2 + 1;
      Double_t step = (_xmax - _xmin) / Double_t(_nBins);
      gr = new TGraph(_nBins);
      for (UInt_t i = 0; i < _nBins; ++i)
         gr->SetPointX(i, _xmin + i * step);

      wfinBG.resize(N);
      wfinSG.resize(N);
      fftout.resize(N);
      wfoutBG = fftw_alloc_complex(M);
      wfoutSG = fftw_alloc_complex(M);
      wfout = fftw_alloc_complex(M);
      FWfftBG = fftw_plan_dft_r2c_1d(N, wfinBG.data(), wfoutBG, FFTW_ESTIMATE);
      FWfftSG = fftw_plan_dft_r2c_1d(N, wfinSG.data(), wfoutSG, FFTW_ESTIMATE);
      BWfft = fftw_plan_dft_c2r_1d(N, wfout, fftout.data(), FFTW_ESTIMATE);
      INITIALIZED = kTRUE;
   }

   ~FFTWState()
   {
      if (gr)
         delete gr;

      if (FWfftBG)
         fftw_destroy_plan(FWfftBG);
      if (FWfftSG)
         fftw_destroy_plan(FWfftSG);
      if (BWfft)
         fftw_destroy_plan(BWfft);

      fftw_free(wfoutBG);
      fftw_free(wfoutSG);
      fftw_free(wfout);
   }
};

thread_local FFTWState g_fftwState;

ClassImp(DFTmethod)

DFTmethod::DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef)
: nPars(4 + _spef.nparams),
   nbins(_nbins),
   xmin(_xmin),
   xmax(_xmax),
   step((_xmax - _xmin) / Double_t(_nbins)),
   edge(_xmin),
   nCalls(0),
   spef(_spef)
{
}

void DFTmethod::CalculateValues()
{
   FFTWState &state = g_fftwState;
   if (!state.IsInitialized())
      state.init(nbins, xmin, step);

   Double_t preFactor = 1.0 / (TMath::Sqrt(TMath::TwoPi()) * s0);
   for (UInt_t i = 0; i < state.N; i++) {
      Double_t xx = xmin + i * step - edge;
      Double_t arg = s0 == 0.0 ? 0.0 : (xx - Q0 + edge) / s0;
      state.wfinBG[i] = preFactor * TMath::Exp(-0.5 * arg * arg);
      state.wfinSG[i] = spef.GetValue(xx);
   }

   fftw_execute(state.FWfftBG);
   fftw_execute(state.FWfftSG);
   for (UInt_t i = 0; i < state.M; i++) {
      Double_t amp_BG = TMath::Hypot(state.wfoutBG[i][0], state.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(state.wfoutBG[i][1], state.wfoutBG[i][0]);
      double ph = (ph_BG + mu * state.wfoutSG[i][1] * step);
      state.wfout[i][0] = amp_BG * TMath::Exp(mu * state.wfoutSG[i][0] * step) * TMath::Cos(ph);
      state.wfout[i][1] = amp_BG * TMath::Exp(mu * state.wfoutSG[i][0] * step) * TMath::Sin(ph);
   }

   fftw_execute(state.BWfft);
   preFactor = Norm * wbin * TMath::Exp(-mu) / Double_t(state.N); // Calculate this once
   for (Int_t i = 0; i < nbins; i++){
      Double_t y = preFactor * state.fftout[i];
      state.gr->SetPointY(i, y);
   }

   return;
}

Double_t DFTmethod::Eval(Double_t *xx, Double_t *pars)
{
   FFTWState &state = g_fftwState;
   if (!state.IsInitialized())
      state.init(nbins, xmin, step);

   SetParameters(pars);
   return state.gr->Eval(xx[0]);
}

void DFTmethod::SetParameters(Double_t *pars)
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
      spef.SetParameters(&pars[4]);
      CalculateValues();
   }
}

TGraph *DFTmethod::GetGraph()
{
   FFTWState &state = g_fftwState;
   if (!state.IsInitialized())
      state.init(nbins, xmin, step);

   CalculateValues();

   TGraph *_gr = static_cast<TGraph *>(state.gr->Clone());
   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);
   return _gr;
}

TGraph *DFTmethod::GetGraphN(Int_t n)
{
   FFTWState &state = g_fftwState;
   if (!state.IsInitialized())
      state.init(nbins, xmin, step);

   CalculateValues();
   Double_t dblN = n;

   Double_t preFactor = 1.0 / (sqrt(2.0 * TMath::Pi()) * s0);
   for (UInt_t i = 0; i < state.N; i++) {
      Double_t xx = xmin + i * step - edge;
      Double_t arg = s0 == 0.0 ? 0.0 : (xx - Q0 + edge) / s0;
      state.wfinBG[i] = preFactor * TMath::Exp(-0.5 * arg * arg);
      state.wfinSG[i] = spef.GetValue(xx);
   }

   fftw_execute(state.FWfftBG);
   fftw_execute(state.FWfftSG);
   for (UInt_t i = 0; i < state.M; i++) {
      Double_t amp_BG = TMath::Hypot(state.wfoutBG[i][0], state.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(state.wfoutBG[i][1], state.wfoutBG[i][0]);
      Double_t amp_SG = TMath::Hypot(state.wfoutSG[i][0], state.wfoutSG[i][1]);
      Double_t ph_SG = TMath::ATan2(state.wfoutSG[i][1], state.wfoutSG[i][0]);
      Double_t ph = (ph_BG + dblN * ph_SG);
      state.wfout[i][0] = amp_BG * pow(step * amp_SG, dblN) * TMath::Cos(ph);
      state.wfout[i][1] = amp_BG * pow(step * amp_SG, dblN) * TMath::Sin(ph);
   }

   fftw_execute(state.BWfft);
   TGraph *_gr = new TGraph(nbins);
   preFactor =  Norm * wbin * TMath::Exp(-mu) * pow(mu, dblN) / (Double_t(state.N) * TMath::Factorial(n)); // Calculate this once
   for (Int_t i = 0; i < nbins; i++) {
      Double_t y = preFactor * state.fftout[i];

      // I don't like this. I'd rather use max to set a floor. This might hide issues
      if (y < 1e-10)
         y = 1e-4;

      _gr->SetPoint(i, xmin + i * step - edge, y);
   }

   _gr->SetLineWidth(2);
   _gr->SetLineStyle(3);
   _gr->SetLineColor(kBlack);
   _gr->SetMarkerColor(kBlack);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

