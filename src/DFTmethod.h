#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include "PMType.h"
#include "IModel.h"

#include "RtypesCore.h"
#include "SPEResponse.h"
#include "TObject.h"
#include "TGraph.h"

#include <fftw3.h>

#include <vector>

struct FFTWState {
private:
   Bool_t isInit{kFALSE}; ///< Is the state initialized

public:
   size_t N{0};                  ///< number of time domain values
   size_t M{0};                  ///< number of frequency domain values
   std::vector<Double_t> wfinBG; ///< FFT Background input
   std::vector<Double_t> wfinSG; ///< FFT SPEResponse input
   std::vector<Double_t> fftout; ///< IFFT real result
   fftw_complex *wfoutBG;        ///< Holds the FFT of the background
   fftw_complex *wfoutSG;        ///< Holds the FFT of the SPEResponse
   fftw_complex *wfout;          ///< Holds the convolution of the above two
   fftw_plan FWfftBG;            ///< FFT Background
   fftw_plan FWfftSG;            ///< FFT Single photon response
   fftw_plan BWfft;              ///< IFFT of the convolution
   TGraph gr;                    ///< TGraph as buffer to store x, y and extrapolate with Eval

   FFTWState() {};

   Bool_t IsInitialized() { return isInit; }

   void init(size_t nBins, Double_t xmin = 0.0, Double_t xmax = 0.0)
   {
      N = 2 * nBins + 60;
      M = N / 2 + 1;
      Double_t step = (xmax - xmin) / Double_t(nBins);
      gr.Set(N);
      for (UInt_t i = 0; i < N; ++i)
         gr.SetPointX(i, xmin + i * step);

      wfinBG.resize(N);
      wfinSG.resize(N);
      fftout.resize(N);
      wfoutBG = fftw_alloc_complex(M);
      wfoutSG = fftw_alloc_complex(M);
      wfout = fftw_alloc_complex(M);
      FWfftBG = fftw_plan_dft_r2c_1d(N, wfinBG.data(), wfoutBG, FFTW_ESTIMATE);
      FWfftSG = fftw_plan_dft_r2c_1d(N, wfinSG.data(), wfoutSG, FFTW_ESTIMATE);
      BWfft = fftw_plan_dft_c2r_1d(N, wfout, fftout.data(), FFTW_ESTIMATE);
      isInit = kTRUE;
   }

   ~FFTWState()
   {
      if (isInit) {
         fftw_destroy_plan(FWfftBG);
         fftw_destroy_plan(FWfftSG);
         fftw_destroy_plan(BWfft);
         fftw_free(wfoutBG);
         fftw_free(wfoutSG);
         fftw_free(wfout);
      }
   }
};

class DFTmethod : public IModel {

public:
   DFTmethod(PMType::Response sper);
   DFTmethod(Int_t nbins, Double_t wbin, Double_t xmin, Double_t xmax, PMType::Response sper);
   DFTmethod(const DFTmethod &other);
   virtual DFTmethod *Clone() const override { return new DFTmethod(*this); }
   virtual ~DFTmethod() = default;

   virtual void SetParameters(const Double_t *pars) override;
   UInt_t GetNCalls() { return m_nCalls; };
   TGraph *GetGraph();
   TGraph *GetGraphN(Int_t n);
   virtual Double_t Gain(const Double_t *pars) const override { return m_resp.Gain(&pars[4]); }
   virtual Double_t GainError(const Double_t *pars, const Double_t *errs) const override
   {
      return m_resp.GainError(&pars[4], &errs[4]);
   }

private:
   virtual Double_t DoEvalPar(Double_t x, const Double_t *pars) const override;
   void CalculateValues(const Double_t *pars) const;

   mutable UInt_t m_nCalls{0};
   mutable SPEResponse m_resp;
   mutable FFTWState fftwState;
};

#endif
