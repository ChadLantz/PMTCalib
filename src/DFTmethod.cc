#include "DFTmethod.h"
#include "SPEResponse.h"
#include "PMType.h"

#include "RtypesCore.h"
#include "TMath.h"

#include <fftw3.h>

////////////////////////////////////////////////////////////////////////////////
/** \class DFTmethod
    \brief Convolution based fit method for low light PMT spectra

## DFTmethod: Convolution based fit method for low light PMT spectra

DFTmethod is the implementation of a FT convolution fitting algorithm. Convolves
a background function with a single photo-electron response function, finally
inverting the transform and storing real values in a TGraph for later lookup.

 */
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Construct with a response type. This is the minimum to use Gain and
/// GainError
DFTmethod::DFTmethod() : IModel(), m_resp()
{
   SetResponse(m_resp.GetResponse());
}

////////////////////////////////////////////////////////////////////////////////
/// Construct with a response type. This is the minimum to use Gain and
/// GainError
DFTmethod::DFTmethod(PMType::Response sper) : IModel(), m_resp(sper)
{
   SetResponse(sper);
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor. Initialize
DFTmethod::DFTmethod(Int_t nbins, Double_t wbin, Double_t xmin, Double_t xmax, PMType::Response sper)
   : IModel(4, nbins, wbin, xmin, xmax), m_nCalls(0), m_resp(sper)
{
   SetNpar(4 + m_resp.NPar());
   SetParNames({"Norm", "Q0", "#sigma_{0}", "#mu"});
   for (UInt_t par = 4; par < NPar(); ++par) {
      ParSettings(par).SetName(m_resp.ParSettings(par - 4).Name());
   }
   fftwState.init(nbins, m_xMin, m_xMax);
}

////////////////////////////////////////////////////////////////////////////////
///
DFTmethod::DFTmethod(const DFTmethod &other) : IModel(other), m_nCalls(0), m_resp(other.m_resp)
{
   fftwState.init(m_nBins, m_xMin, m_xMax);
}

////////////////////////////////////////////////////////////////////////////////
///
void DFTmethod::CalculateValues(const Double_t *pars) const
{
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];

   for (UInt_t i = 0; i < fftwState.N; i++) {
      Double_t x = fftwState.gr.GetPointX(i);
      fftwState.wfinBG[i] = TMath::Gaus(x, Q0, s0, kTRUE);
      fftwState.wfinSG[i] = m_resp(x - m_xMin, &pars[4]);
   }

   fftw_execute(fftwState.FWfftBG);
   fftw_execute(fftwState.FWfftSG);
   for (UInt_t i = 0; i < fftwState.M; i++) {
      Double_t amp_BG = TMath::Hypot(fftwState.wfoutBG[i][0], fftwState.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(fftwState.wfoutBG[i][1], fftwState.wfoutBG[i][0]);
      double ph = (ph_BG + mu * fftwState.wfoutSG[i][1] * m_step);
      fftwState.wfout[i][0] = amp_BG * TMath::Exp(mu * fftwState.wfoutSG[i][0] * m_step) * TMath::Cos(ph);
      fftwState.wfout[i][1] = amp_BG * TMath::Exp(mu * fftwState.wfoutSG[i][0] * m_step) * TMath::Sin(ph);
   }

   fftw_execute(fftwState.BWfft);
   Double_t normConst = Norm * m_wBin * TMath::Exp(-mu) / Double_t(fftwState.N);
   for (Int_t i = 0; i < fftwState.N; i++) {
      fftwState.gr.SetPointY(i, normConst * fftwState.fftout[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
///
void DFTmethod::SetParameters(const Double_t *pars)
{
   Bool_t parsChanged = kFALSE;
   for (UInt_t ipar = 0; ipar < NPar(); ++ipar) {
      if (pars[ipar] != m_parameters[ipar]) {
         m_parameters[ipar] = pars[ipar];
         m_parSettings[ipar].SetValue(pars[ipar]);
         parsChanged = kTRUE;
      }
   }

   if (parsChanged) {
      if (!fftwState.IsInitialized())
         fftwState.init(m_nBins, m_xMin, m_xMax);
      m_nCalls++;
      m_resp.SetParameters(&pars[4]);
      CalculateValues(pars);
   }
}

void DFTmethod::SetResponse(PMType::Response resp)
{
   m_resp.SetResponse(resp);
   SetNpar(4 + m_resp.NPar());
   for (UInt_t par = 4; par < NPar(); ++par) {
      ParSettings(par).SetName(m_resp.ParSettings(par - 4).Name());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate the model with the given parameters
Double_t DFTmethod::DoEvalPar(const Double_t x, const Double_t *pars) const
{
   Bool_t parsChanged = kFALSE;
   for (UInt_t ipar = 0; ipar < NPar(); ++ipar) {
      if (pars[ipar] != m_parameters[ipar]) {
         m_parameters[ipar] = pars[ipar];
         m_parSettings[ipar].SetValue(pars[ipar]);
         parsChanged = kTRUE;
      }
   }

   if (parsChanged) {
      if (!fftwState.IsInitialized())
         fftwState.init(m_nBins, m_xMin, m_xMax);
      m_nCalls++;
      m_resp.SetParameters(&pars[4]);
      CalculateValues(pars);
   }
   return fftwState.gr.Eval(x);
};

////////////////////////////////////////////////////////////////////////////////
///
TGraph *DFTmethod::GetGraph()
{
   CalculateValues(Parameters());

   TGraph *_gr = static_cast<TGraph *>(fftwState.gr.Clone());
   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);
   return _gr;
}

////////////////////////////////////////////////////////////////////////////////
///
TGraph *DFTmethod::GetGraphN(Int_t n)
{
   CalculateValues(Parameters());
   Double_t dblN = n;
   Double_t Norm = m_parameters[0];
   Double_t Q0 = m_parameters[1];
   Double_t s0 = m_parameters[2];
   Double_t mu = m_parameters[3];

   for (UInt_t i = 0; i < fftwState.N; i++) {
      Double_t x = fftwState.gr.GetPointX(i);
      fftwState.wfinBG[i] = TMath::Gaus(x, Q0, s0, kTRUE);
      fftwState.wfinSG[i] = m_resp(x, &m_parameters[4]);
   }

   fftw_execute(fftwState.FWfftBG);
   fftw_execute(fftwState.FWfftSG);
   for (UInt_t i = 0; i < fftwState.M; i++) {
      Double_t amp_BG = TMath::Hypot(fftwState.wfoutBG[i][0], fftwState.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(fftwState.wfoutBG[i][1], fftwState.wfoutBG[i][0]);
      Double_t amp_SG = TMath::Hypot(fftwState.wfoutSG[i][0], fftwState.wfoutSG[i][1]);
      Double_t ph_SG = TMath::ATan2(fftwState.wfoutSG[i][1], fftwState.wfoutSG[i][0]);
      Double_t ph = (ph_BG + dblN * ph_SG);
      fftwState.wfout[i][0] = amp_BG * pow(m_step * amp_SG, dblN) * TMath::Cos(ph);
      fftwState.wfout[i][1] = amp_BG * pow(m_step * amp_SG, dblN) * TMath::Sin(ph);
   }

   fftw_execute(fftwState.BWfft);
   TGraph *_gr = new TGraph(fftwState.N);
   Double_t normConst = Norm * m_wBin * TMath::Exp(-mu) * pow(mu, dblN) /
                        (Double_t(fftwState.N) * TMath::Factorial(n)); // Calculate this once
   for (Int_t i = 0; i < fftwState.N; i++) {
      Double_t y = normConst * fftwState.fftout[i];

      // I don't like this. I'd rather use max to set a floor. This might hide issues
      if (y < 1e-10)
         y = 1e-4;

      _gr->SetPoint(i, fftwState.gr.GetPointX(i), y);
   }

   _gr->SetLineWidth(2);
   _gr->SetLineStyle(3);
   _gr->SetLineColor(kBlack);
   _gr->SetMarkerColor(kBlack);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

Double_t DFTmethod::Gain(const Double_t *pars) const
{
   Double_t Q0 = pars[1];
   Double_t mu = pars[3];
   Double_t Q = 0.0;
   switch (m_resp.GetResponse()) {
   case PMType::Response::GAMMA:
   case PMType::Response::GAMMA2EXP: Q = 1.0 / pars[4]; break;
   default: Q = pars[4]; break;
   }
   return Q;
}


Double_t DFTmethod::GainError(const Double_t *pars, const Double_t *errs) const
{
   Double_t Q = 0.0;
   switch (m_resp.GetResponse()) {
   case PMType::Response::GAMMA:
   case PMType::Response::WEIBULL:
   case PMType::Response::GAMMA2EXP: Q = 1.0 / errs[4]; break;
   default: Q = errs[4]; break;
   }
   return Q;
}