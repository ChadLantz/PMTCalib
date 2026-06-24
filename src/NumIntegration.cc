
#include "NumIntegration.h"
#include "PMType.h"

ClassImp(NumIntegration)

   NumIntegration::NumIntegration()
{
}

NumIntegration::~NumIntegration() {}

NumIntegration::NumIntegration(Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Response _sper, Double_t *_params)
{
   N = nbins = _nbins;
   xmin = _xmin;
   xmax = _xmax;
   step = (xmax - xmin) / (1.0 * nbins * 1.0);
   SPEResponseFactory sperf;
   spef = sperf.Build(_sper, _params);

   xvalues.clear();
   for (UInt_t i = 0; i < N; i++) {
      Double_t xx = xmin + wbin / 2.0 + 1.0 * i * step;
      xvalues.push_back(xx);
   }
}

void NumIntegration::CalculateValues()
{
   Double_t x[N];
   Double_t y[N];

   yvalues.clear();
   for (UInt_t i = 0; i < N; i++) {
      Double_t xx = xvalues.at(i);
      x[i] = xx;

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
            SR1 += yy0 * spef.Eval(xt) * ds;
      }

      SR1 *= TMath::Poisson(1, mu);
      result += SR1;

      Double_t Q = spef.GetParameter(0);
      Double_t s = spef.GetParameter(1);
      Double_t alpha = spef.GetParameter(2);
      Double_t w = spef.GetParameter(3);

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

      result *= wbin * Norm;

      yvalues.push_back(result);
      y[i] = result;
   }

   gr = new TGraph(N, x, y);

   return;
}

Double_t NumIntegration::GetValue(Double_t xx)
{
   Double_t y_ = gr->Eval(xx);

   return y_;
}

TGraph *NumIntegration::GetGraph()
{
   CalculateValues();

   Double_t x[N];
   Double_t y[N];

   for (UInt_t i = 0; i < N; i++) {
      x[i] = xvalues.at(i);
      Double_t y_ = GetValue(x[i]);

      if (y_ < 1.0e-10)
         y[i] = 1.e-4;
      else
         y[i] = y_;
   }

   TGraph *_gr = new TGraph(N, x, y);

   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);

   return _gr;
}