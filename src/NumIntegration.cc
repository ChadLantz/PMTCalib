#include "NumIntegration.h"
#include "PMType.h"

#include "TObject.h"
#include "TMath.h"
#include "RtypesCore.h"

NumIntegration::NumIntegration(Int_t nbins, Double_t wbin, Double_t xmin, Double_t xmax, PMType::Response sper)
   : IModel(8, wbin, nbins, xmin, xmax),
     m_resp(sper),
     m_gr(new TGraph(nbins))
{
   SetParNames({"Norm", "Q0", "#sigma_{0}", "#mu"});
   SetNpar(4 + m_resp.NPar());
   for (UInt_t par = 4; par < NPar(); ++par) {
      ParSettings(par).SetName(m_resp.ParSettings(par - 4).Name());
   }

   for (UInt_t i = 0; i < m_nBins; i++) {
      m_gr->SetPointX(i, m_xMin + 0.5 * m_wBin + Double_t(i) * m_step);
   }
}

NumIntegration::NumIntegration(const NumIntegration &other)
   : IModel(other),
     m_nCalls(0),
     m_resp(other.m_resp),
     m_gr(static_cast<TGraph *>(other.m_gr->Clone()))
{
}

NumIntegration::~NumIntegration()
{
   if (m_gr)
      delete m_gr;
}

void NumIntegration::CalculateValues(const Double_t *pars) const
{
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];
   Double_t Q = pars[4];
   Double_t s = pars[5];
   Double_t alpha = pars[6];
   Double_t w = pars[7];

   for (UInt_t i = 0; i < m_nBins; i++) {
      Double_t xx = m_gr->GetPointX(i);

      Double_t result = 0.0;

      Double_t arg = 0.0;
      if (s0 != 0.0)
         arg = (xx - Q0) / s0;
      else
         Error("CalculateValues", "Division by zero(i): sigma0");
      Double_t SR0 = TMath::Exp(-mu) / (sqrt(2.0 * TMath::Pi()) * s0) * TMath::Exp(-0.5 * arg * arg);
      result += SR0;

      Double_t x_lo = Q0 - 7.0 * s0;
      Double_t x_hi = Q0 + 7.0 * s0;
      UInt_t nds = 100;
      Double_t ds = (x_hi - x_lo) / nds;

      Double_t SR1 = 0.0;

      for (UInt_t j = 0; j < nds; j++) {
         Double_t tt = x_lo + j * ds;
         Double_t xt = xx - tt;

         Double_t arg = 0.0;
         if (s0 != 0.0)
            arg = (tt - Q0) / s0;
         else
            Error("CalculateValues", "Division by zero(j): sigma0");
         Double_t yy0 = 1.0 / (sqrt(2.0 * TMath::Pi()) * s0) * TMath::Exp(-0.5 * arg * arg);

         if (xt >= 0.0)
            SR1 += yy0 * m_resp(xt, &pars[4]) * ds;
      }

      SR1 *= TMath::Poisson(1, mu);
      result += SR1;

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));
      Double_t k = s / gn / sqrt(2.0 * TMath::Pi()) * TMath::Exp(-pow(Q, 2.0) / (2.0 * pow(s, 2.0)));
      Double_t Qg = Q + k;
      Double_t sg2 = pow(s, 2.0) - (Q + k) * k;

      Double_t Qs = w / alpha + (1.0 - w) * Qg;
      Double_t ss2 = w / pow(alpha, 2.0) + (1 - w) * sg2 + w * (1.0 - w) * pow(Qg - 1.0 / alpha, 2.0);

      for (Int_t n = 2; n < 50; n++) {
         Double_t Qn = Q0 + 1.0 * n * Qs;
         Double_t sn2 = pow(s0, 2.0) + 1.0 * n * ss2;
         Double_t sn = sqrt(sn2);

         Double_t argn = 0.0;
         if (sn != 0.0)
            argn = (xx - Qn) / sn;
         else
            Error("CalculateValues", "Division by zero: sigmaN");
         Double_t SRn = 1.0 / (sqrt(2.0 * TMath::Pi()) * sn) * TMath::Exp(-0.5 * argn * argn);

         SRn *= TMath::Poisson(n, mu);
         result += SRn;
      }

      result *= m_wBin * Norm;
      m_gr->SetPointY(i, result);
   }
}

void NumIntegration::SetParameters(const Double_t *pars)
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
      m_nCalls++;
      m_resp.SetParameters(&pars[4]);
      CalculateValues(pars);
   }
}

Double_t NumIntegration::DoEvalPar(Double_t x, const Double_t *pars) const
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
      m_nCalls++;
      m_resp.SetParameters(&pars[4]);
      CalculateValues(pars);
   }
   return m_gr->Eval(x);
}

TGraph *NumIntegration::GetGraph()
{
   CalculateValues(Parameters());

   for (UInt_t i = 0; i < m_nBins; i++) {
      Double_t y_ = m_gr->GetPointY(i);

      if (y_ < 1.0e-10)
         m_gr->SetPointY(i, 1.e-4);
   }

   TGraph *_gr = static_cast<TGraph *>(m_gr->Clone());
   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);

   return _gr;
}