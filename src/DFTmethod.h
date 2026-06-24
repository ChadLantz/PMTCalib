#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include "PMType.h"

#include "RtypesCore.h"
#include "TF1.h"
#include "TObject.h"
#include "TGraph.h"

#include <fftw3.h>

#include <vector>


struct FFTWState {
   size_t N;                     ///< number of time domain values
   size_t M;                     ///< number of frequency domain values
   std::vector<Double_t> wfinBG; ///< FFT Background input
   std::vector<Double_t> wfinSG; ///< FFT SPEResponse input
   std::vector<Double_t> fftout; ///< IFFT real result
   fftw_complex *wfoutBG;        ///< Holds the FFT of the background
   fftw_complex *wfoutSG;        ///< Holds the FFT of the SPEResponse
   fftw_complex *wfout;          ///< Holds the convolution of the above two
   fftw_plan FWfftBG;            ///< FFT Background
   fftw_plan FWfftSG;            ///< FFT Single photon response
   fftw_plan BWfft;              ///< IFFT of the convolution
   TGraph *gr;                   ///< This state owns its graph and we will return clones

   FFTWState() {};

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

class DFTmethod : public TObject {

public:
   DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Response _sper, Double_t *_params);
   DFTmethod(const DFTmethod &other);
   virtual ~DFTmethod();

   Double_t Eval(Double_t x) { return fftwState.gr ? fftwState.gr->Eval(x) : 0.0; }
   Double_t EvalPar(Double_t *xx, Double_t *pars);
   Double_t operator()(Double_t *xx, Double_t *pars) { return EvalPar(xx, pars); }
   void SetParameters(Double_t *pars);
   TGraph *GetGraph();
   TGraph *GetGraphN(Int_t n);
   UInt_t GetNCalls() { return nCalls; };
   UInt_t GetNpar() { return nPars; };
   const char *GetParName(UInt_t par);

private:
   void CalculateValues();

   UInt_t nPars;
   UInt_t nCalls;
   Int_t nbins;
   Double_t xmin;
   Double_t xmax;
   Double_t step;
   Double_t edge;
   Double_t wbin;
   Double_t Norm;
   Double_t Q0;
   Double_t s0;
   Double_t mu;
   Double_t parCache[11];
   TF1 spef;
   FFTWState fftwState;
   std::vector<const char *> parNames;

   ClassDef(DFTmethod, 1)
};

#endif
