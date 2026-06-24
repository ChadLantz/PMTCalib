
#include "SPEResponseFactory.h"

#include "TMath.h"

ClassImp(SPEResponseFactory)

static Double_t _gausexpfunc(Double_t *x, Double_t *par)
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

static Double_t _gaus2expfunc(Double_t *x, Double_t *par)
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

static Double_t _gammaexpfunc(Double_t *x, Double_t *par)
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

static Double_t _gamma2expfunc(Double_t *x, Double_t *par)
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

static Double_t _weibullexpfunc(Double_t *x, Double_t *par)
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

static Double_t _lognormalexpfunc(Double_t *x, Double_t *par)
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

static Double_t _testfunc(Double_t *x, Double_t *par)
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


////////////////////////////////////////////////////////////////////////////////
/// Main constructor. Initializes member TF1 with the appropriate fit function,
/// parameters, parameter names, and range
///
/// @param _spetype Single photo-electron response function enum
/// @param _params Array of initial parameters
TF1 SPEResponseFactory::Build(PMType::Response _spetype, Double_t _params[])
{
   Double_t (*funcPtr)(Double_t *, Double_t *);
   Double_t xMin = 0.0, xMax = 0.0;
   std::vector<std::string> parNames;

   switch (_spetype) {
   case PMType::GAUSS:
      funcPtr = _gausexpfunc;
      xMin = _params[0] - 80.0 * _params[1];
      xMax = _params[0] + 80.0 * _params[1];
      parNames = {"Q", "#sigma", "#alpha", "w"};
      break;

   case PMType::GAMMA:
      funcPtr = _gammaexpfunc;
      xMin = 1.0 / _params[0] - 80.0 * 1.0 / _params[0] / sqrt(1.0 + _params[1]);
      xMax = 1.0 / _params[0] + 80.0 * 1.0 / _params[0] / sqrt(1.0 + _params[1]);
      parNames = {"#lambda", "#kappa", "#alpha", "w"};
      break;

   case PMType::WEIBULL:
      funcPtr = _weibullexpfunc;
      xMin = _params[0] - 80.0 * _params[0];
      xMax = _params[0] + 80.0 * _params[0];
      parNames = {"#lambda", "#kappa", "#alpha", "w"};
      break;

   case PMType::LOGNORMAL:
      funcPtr = _lognormalexpfunc;
      xMin = TMath::Log(_params[0]) - 80.0 * TMath::Log(_params[0]);
      xMax = TMath::Log(_params[0]) + 80.0 * TMath::Log(_params[0]);
      parNames = {"Q", "#sigma", "#alpha", "w"};
      break;

   case PMType::GAUSS2EXP:
      funcPtr = _gaus2expfunc;
      xMin = _params[0] - 80.0 * _params[1];
      xMax = _params[0] + 80.0 * _params[1];
      parNames = {"Q", "#sigma", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}"};
      break;

   case PMType::GAMMA2EXP:
      funcPtr = _gamma2expfunc;
      xMin = 1.0 / _params[0] - 80.0 * 1.0 / _params[0] / sqrt(1.0 + _params[1]);
      xMax = 1.0 / _params[0] + 80.0 * 1.0 / _params[0] / sqrt(1.0 + _params[1]);
      parNames = {"#lambda", "#kappa", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}"};
      break;

   case PMType::TEST:
      funcPtr = _testfunc;
      xMin = _params[0] - 80.0 * _params[1];
      xMax = _params[0] + 80.0 * _params[1];
      parNames = {"Q", "#sigma", "#lambda", "#theta", "w"};
      break;
   }

   TF1 spefunc("spefunc", funcPtr, xMin, xMax, parNames.size());
   for(UInt_t par = 0; par < parNames.size(); ++par){
      spefunc.SetParName(par, parNames[par].c_str());
   }
   spefunc.SetParameters(_params);
   return spefunc;
}