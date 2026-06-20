
#include "SPEResponse.h"

#include "TMath.h"

ClassImp(SPEResponse)

Double_t _gausexpfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (xx - Q) / s;
      else
         Error("_gausexpfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      result = w * alpha * TMath::Exp(-xx * alpha) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

Double_t _gaus2expfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha1 = par[2];
   Double_t w1 = par[3];
   Double_t alpha2 = par[4];
   Double_t w2 = par[5];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (xx - Q) / s;
      else
         Error("_gaus2expfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      result = w1 * alpha1 * TMath::Exp(-xx * alpha1) + w2 * alpha2 * TMath::Exp(-xx * alpha2) +
               (1.0 - w1 - w2) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

Double_t _gammaexpfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t lambda = par[0];
   Double_t theta = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * xx;

      result = w * 1.0 * alpha * TMath::Exp(-xx * alpha) +
               (1.0 - w) * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
   }

   return result;
}

Double_t _gamma2expfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t lambda = par[0];
   Double_t theta = par[1];
   Double_t alpha1 = par[2];
   Double_t w1 = par[3];
   Double_t alpha2 = par[4];
   Double_t w2 = par[5];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * xx;

      result = w1 * 1.0 * alpha1 * TMath::Exp(-xx * alpha1) + w2 * 1.0 * alpha2 * TMath::Exp(-xx * alpha2) +
               (1.0 - w1 - w2) * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
   }

   return result;
}

Double_t _weibullexpfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t lambda = par[0];
   Double_t kappa = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t f = kappa / lambda;
      Double_t fx = xx / lambda;

      result = w * alpha * TMath::Exp(-xx * alpha) + (1.0 - w) * f * pow(fx, kappa - 1.0) * TMath::Exp(-pow(fx, kappa));
   }

   return result;
}

Double_t _lognormalexpfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (xx > 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (TMath::Log(xx) - Q) / s;
      else
         Error("_lognormalexpfunc", "Division by zero: sigma");

      result = w * alpha * TMath::Exp(-xx * alpha) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * xx) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

Double_t _testfunc(Double_t *x, Double_t *par)
{
   Double_t xx = x[0];
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t lambda = par[2];
   Double_t theta = par[4];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (xx >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (xx - Q) / s;
      else
         Error("_testfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * xx;

      result = w * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

SPEResponse::SPEResponse(): spefunc(nullptr) {}

SPEResponse::~SPEResponse() {}

SPEResponse::SPEResponse(PMType::Response _spetype, Double_t _params[])
{
   spetype = _spetype;
   spefunc = GetFitFunc();
   SetParams(_params);

   spefunc->SetLineColor(kBlue);
   spefunc->SetLineWidth(2.0);
   spefunc->SetNpx(10000);
}

TF1 *SPEResponse::GetFitFunc()
{
  if (spefunc)
    return spefunc;

   Double_t (*funcPtr)(Double_t *, Double_t *);
   Double_t xMin = 0.0, xMax = 0.0;

   switch (spetype) {
   case PMType::GAUSS:
      funcPtr = _gausexpfunc;
      nparams = 4;
      xMin = params[0] - 80.0 * params[1];
      xMax = params[0] + 80.0 * params[1];
      spefunc->SetParNames("Q", "s", "#alpha", "w");
      break;

   case PMType::GAMMA:
      funcPtr = _gammaexpfunc;
      nparams = 4;
      xMin = 1.0 / params[0] - 80.0 * 1.0 / params[0] / sqrt(1.0 + params[1]);
      xMax = 1.0 / params[0] + 80.0 * 1.0 / params[0] / sqrt(1.0 + params[1]);
      spefunc->SetParNames("#lambda", "#kappa", "#alpha", "w");
      break;

   case PMType::WEIBULL:
      funcPtr = _weibullexpfunc;
      nparams = 4;
      xMin = params[0] - 80.0 * params[0];
      xMax = params[0] + 80.0 * params[0];
      spefunc->SetParNames("#lambda", "#kappa", "#alpha", "w");
      break;

   case PMType::LOGNORMAL:
      funcPtr = _lognormalexpfunc;
      nparams = 4;
      xMin = TMath::Log(params[0]) - 80.0 * TMath::Log(params[0]);
      xMax = TMath::Log(params[0]) + 80.0 * TMath::Log(params[0]);
      spefunc->SetParNames("Q", "s", "#alpha", "w");
      break;

   case PMType::GAUSS2EXP:
      funcPtr = _gaus2expfunc;
      nparams = 6;
      xMin = params[0] - 80.0 * params[1];
      xMax = params[0] + 80.0 * params[1];
      spefunc->SetParNames("Q", "s", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}");
      break;

   case PMType::GAMMA2EXP:
      funcPtr = _gamma2expfunc;
      nparams = 6;
      xMin = 1.0 / params[0] - 80.0 * 1.0 / params[0] / sqrt(1.0 + params[1]);
      xMax = 1.0 / params[0] + 80.0 * 1.0 / params[0] / sqrt(1.0 + params[1]);
      spefunc->SetParNames("#lambda", "#kappa", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}");
      break;

   case PMType::TEST:
      funcPtr = _testfunc;
      nparams = 5;
      xMin = params[0] - 80.0 * params[1];
      xMax = params[0] + 80.0 * params[1];
      spefunc->SetParNames("Q", "s", "#lambda", "#theta", "w");
      break;
   }

   return new TF1("spefunc", funcPtr, xMin, xMax, nparams);
}

void SPEResponse::SetParams(Double_t _params[])
{
   for (Int_t i = 0; i < nparams; i++) {
      params[i] = _params[i];
      spefunc->SetParameter(i, params[i]);
   }
}

Double_t SPEResponse::GetValue(Double_t xx)
{
   Double_t result = spefunc->Eval(xx);

   return result;
}

Double_t SPEResponse::GenQ()
{
   Double_t _x = spefunc->GetRandom();

   return _x;
}
