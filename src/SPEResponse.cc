#include "SPEResponse.h"

#include "PMType.h"
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////
/// Default constructor
SPEResponse::SPEResponse()
: IModel(), m_resp(PMType::Response::GAUSS)
{
   SetResponse(m_resp);
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor
SPEResponse::SPEResponse(PMType::Response resp)
: IModel()
{
   SetResponse(resp);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor
SPEResponse::SPEResponse(const SPEResponse &other)
: IModel(other), m_resp(other.m_resp)
{
   SetResponse(m_resp);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the response function and parameter names accordingly
void SPEResponse::SetResponse(PMType::Response resp)
{
   m_resp = resp;
   switch (m_resp) {
   case PMType::Response::GAUSS:
      SetNpar(4);
      SetParNames({"Q", "#sigma", "#alpha", "w"});
      break;

   case PMType::Response::GAMMA:
      SetNpar(4);
      SetParNames({"#lambda", "#kappa", "#alpha", "w"});
      break;

   case PMType::Response::WEIBULL:
      SetNpar(4);
      SetParNames({"#lambda", "#kappa", "#alpha", "w"});
      break;

   case PMType::Response::LOGNORMAL:
      SetNpar(4);
      SetParNames({"Q", "#sigma", "#alpha", "w"});
      break;

   case PMType::Response::GAUSS2EXP:
      SetNpar(6);
      SetParNames({"Q", "#sigma", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}"});
      break;

   case PMType::Response::GAMMA2EXP:
      SetNpar(6);
      SetParNames({"#lambda", "#kappa", "#alpha_{1}", "w_{1}", "#alpha_{2}", "w_{2}"});
      break;

   case PMType::Response::TEST:
      SetNpar(5);
      SetParNames({"Q", "#sigma", "#lambda", "#theta", "w"});
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::gausexpfunc(Double_t x, const Double_t *par) const
{
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (x - Q) / s;
      else
         Error("_gausexpfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      result = w * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::gaus2expfunc(Double_t x, const Double_t *par) const
{
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha1 = par[2];
   Double_t w1 = par[3];
   Double_t alpha2 = par[4];
   Double_t w2 = par[5];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (x - Q) / s;
      else
         Error("_gaus2expfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      result = w1 * alpha1 * TMath::Exp(-x * alpha1) + w2 * alpha2 * TMath::Exp(-x * alpha2) +
               (1.0 - w1 - w2) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::gammaexpfunc(Double_t x, const Double_t *par) const
{
   Double_t lambda = par[0];
   Double_t theta = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * x;

      result = w * 1.0 * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::gamma2expfunc(Double_t x, const Double_t *par) const
{
   Double_t lambda = par[0];
   Double_t theta = par[1];
   Double_t alpha1 = par[2];
   Double_t w1 = par[3];
   Double_t alpha2 = par[4];
   Double_t w2 = par[5];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * x;

      result = w1 * 1.0 * alpha1 * TMath::Exp(-x * alpha1) + w2 * 1.0 * alpha2 * TMath::Exp(-x * alpha2) +
               (1.0 - w1 - w2) * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::weibullexpfunc(Double_t x, const Double_t *par) const
{
   Double_t lambda = par[0];
   Double_t kappa = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t f = kappa / lambda;
      Double_t fx = x / lambda;

      result = w * alpha * TMath::Exp(-x * alpha) + (1.0 - w) * f * pow(fx, kappa - 1.0) * TMath::Exp(-pow(fx, kappa));
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::lognormalexpfunc(Double_t x, const Double_t *par) const
{
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t alpha = par[2];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (x > 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (TMath::Log(x) - Q) / s;
      else
         Error("_lognormalexpfunc", "Division by zero: sigma");

      result = w * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * x) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::testfunc(Double_t x, const Double_t *par) const
{
   Double_t Q = par[0];
   Double_t s = par[1];
   Double_t lambda = par[2];
   Double_t theta = par[4];
   Double_t w = par[3];

   Double_t result = 0.0;
   if (x >= 0.0) {
      Double_t arg = 0.0;
      if (s != 0.0)
         arg = (x - Q) / s;
      else
         Error("_testfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * x;

      result = w * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx) +
               (1.0 - w) / (sqrt(2.0 * TMath::Pi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// 
Double_t SPEResponse::DoEvalPar(Double_t x, const Double_t *pars) const
{
   Double_t retVal = 0.0;
   switch (m_resp) {
   case PMType::Response::GAUSS: retVal = gausexpfunc(x, pars); break;
   case PMType::Response::GAMMA: retVal = gammaexpfunc(x, pars); break;
   case PMType::Response::WEIBULL: retVal = weibullexpfunc(x, pars); break;
   case PMType::Response::LOGNORMAL: retVal = lognormalexpfunc(x, pars); break;
   case PMType::Response::GAUSS2EXP: retVal = gaus2expfunc(x, pars); break;
   case PMType::Response::GAMMA2EXP: retVal = gamma2expfunc(x, pars); break;
   case PMType::Response::TEST: retVal = testfunc(x, pars); break;
   default:
      Error("SPEResponse::DoEvalPar", "Invalid response type selected m_resp = %d", m_resp);
      break;
   }

   // Info("SPEResponse::gammaexpfunc", "Returned %.2f for x = %.2f", retVal, x);
   // for(UInt_t ipar = 0; ipar < NPar(); ++ipar){
   //    Info("SPEResponse::gammaexpfunc", "%s = %.2f", ParameterName(ipar).c_str(), pars[ipar]);
   // }
   return retVal;
}
