#include "DFTmethod.h"
#include "SPEResponseFactory.h"
#include "PMType.h"

#include "RtypesCore.h"
#include "SPEResponseFactory.h"
#include "TError.h" // IWYU pragma: keep
#include "TF1.h"
#include "TMath.h"

#include <fftw3.h>

////////////////////////////////////////////////////////////////////////////////
/** \class DFTmethod
    \brief Convolution based fit method for low light PMT spectra

## DFTmethod: Convolution based fit method for low light PMT spectra

DFTmethod is the implementation of a FT convolution fitting algorithm. Convolves
a background function with a single photo-electron response function before
applying another transfer function(?), finally inverting the transform and
storing real values in a TGraph for later lookup.
 */
////////////////////////////////////////////////////////////////////////////////

ClassImp(DFTmethod)

DFTmethod::DFTmethod(Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Response _sper, Double_t *_params)
:  nPars(5),
   nbins(_nbins),
   xmin(_xmin),
   xmax(_xmax),
   step((_xmax - _xmin) / Double_t(_nbins)),
   edge(_xmin),
   nCalls(0),
   parNames({"wbin", "Norm", "Q0", "#sigma_{0}", "#mu"})
{
   SPEResponseFactory sperf;
   spef = sperf.Build(_sper, _params);
   nPars = 5 + spef.GetNpar();
   fftwState.init(nbins, xmin, xmax);
}

DFTmethod::DFTmethod(const DFTmethod &other)
:  nPars(other.nPars),
   nCalls(0),
   nbins(other.nbins),
   xmin(other.xmin),
   xmax(other.xmax),
   step(other.step),
   edge(other.edge),
   wbin(other.wbin),
   Norm(other.Norm),
   Q0(other.Q0),
   s0(other.s0),
   mu(other.mu),
   spef(other.spef),
   parNames(other.parNames)
{
   // Don't copy the fftwState or parCache
}

DFTmethod::~DFTmethod(){
}

void DFTmethod::CalculateValues()
{
   Double_t normConst = 1.0 / (TMath::Sqrt(TMath::TwoPi()) * s0);
   for (UInt_t i = 0; i < fftwState.N; i++) {
      Double_t xx = xmin + i * step - edge;
      Double_t arg = s0 == 0.0 ? 0.0 : (xx - Q0 + edge) / s0;
      fftwState.wfinBG[i] = normConst * TMath::Exp(-0.5 * arg * arg);
      fftwState.wfinSG[i] = spef.Eval(xx);
   }

   fftw_execute(fftwState.FWfftBG);
   fftw_execute(fftwState.FWfftSG);
   for (UInt_t i = 0; i < fftwState.M; i++) {
      Double_t amp_BG = TMath::Hypot(fftwState.wfoutBG[i][0], fftwState.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(fftwState.wfoutBG[i][1], fftwState.wfoutBG[i][0]);
      double ph = (ph_BG + mu * fftwState.wfoutSG[i][1] * step);
      fftwState.wfout[i][0] = amp_BG * TMath::Exp(mu * fftwState.wfoutSG[i][0] * step) * TMath::Cos(ph);
      fftwState.wfout[i][1] = amp_BG * TMath::Exp(mu * fftwState.wfoutSG[i][0] * step) * TMath::Sin(ph);
   }

   fftw_execute(fftwState.BWfft);
   normConst = Norm * wbin * TMath::Exp(-mu) / Double_t(fftwState.N); // Calculate this once
   for (Int_t i = 0; i < nbins; i++){
      fftwState.gr->SetPointY(i, normConst * fftwState.fftout[i]);
   }
}

Double_t DFTmethod::EvalPar(Double_t *xx, Double_t *pars)
{
   SetParameters(pars);
   return Eval(xx[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Tests if parameters have been changed since the last CalculateValues call.
/// If so, set the parameters and recalculate values. Parameter order is
/// BinWidth - Width of the histogram bins
/// Norm - Histogram population
/// Q0 - Pedestal mean
/// s0 - Pedestal width
/// mu - Mean PE count
/// pars[par > 4] - SPEResponse parameters
/// 
/// @param pars Fit parameters for DFT method and its SPEResponse
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
      wbin = pars[0];
      Norm = pars[1];
      Q0 = pars[2];
      s0 = pars[3];
      mu = pars[4];
      spef.SetParameters(&pars[5]);
      CalculateValues();
   }
}

TGraph *DFTmethod::GetGraph()
{
   CalculateValues();

   TGraph *_gr = static_cast<TGraph *>(fftwState.gr->Clone());
   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);
   return _gr;
}

TGraph *DFTmethod::GetGraphN(Int_t n)
{
   CalculateValues();
   Double_t dblN = n;

   Double_t normConst = 1.0 / (sqrt(2.0 * TMath::Pi()) * s0);
   for (UInt_t i = 0; i < fftwState.N; i++) {
      Double_t xx = xmin + i * step - edge;
      Double_t arg = s0 == 0.0 ? 0.0 : (xx - Q0 + edge) / s0;
      fftwState.wfinBG[i] = normConst * TMath::Exp(-0.5 * arg * arg);
      fftwState.wfinSG[i] = spef.Eval(xx);
   }

   fftw_execute(fftwState.FWfftBG);
   fftw_execute(fftwState.FWfftSG);
   for (UInt_t i = 0; i < fftwState.M; i++) {
      Double_t amp_BG = TMath::Hypot(fftwState.wfoutBG[i][0], fftwState.wfoutBG[i][1]);
      Double_t ph_BG = TMath::ATan2(fftwState.wfoutBG[i][1], fftwState.wfoutBG[i][0]);
      Double_t amp_SG = TMath::Hypot(fftwState.wfoutSG[i][0], fftwState.wfoutSG[i][1]);
      Double_t ph_SG = TMath::ATan2(fftwState.wfoutSG[i][1], fftwState.wfoutSG[i][0]);
      Double_t ph = (ph_BG + dblN * ph_SG);
      fftwState.wfout[i][0] = amp_BG * pow(step * amp_SG, dblN) * TMath::Cos(ph);
      fftwState.wfout[i][1] = amp_BG * pow(step * amp_SG, dblN) * TMath::Sin(ph);
   }

   fftw_execute(fftwState.BWfft);
   TGraph *_gr = new TGraph(nbins);
   normConst =  Norm * wbin * TMath::Exp(-mu) * pow(mu, dblN) / (Double_t(fftwState.N) * TMath::Factorial(n)); // Calculate this once
   for (Int_t i = 0; i < nbins; i++) {
      Double_t y = normConst * fftwState.fftout[i];

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

const char *DFTmethod::GetParName(UInt_t par){
   if(par < parNames.size()){
      return parNames[par];
   }else if (par < nPars){
      return spef.GetParName(par - 5);
   }else{
      Error("GetParName", "Requested parameter %u is greater than the number of parameters %u", par, nPars);
   }
   return "";
}