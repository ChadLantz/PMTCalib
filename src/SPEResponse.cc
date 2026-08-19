#include "SPEResponse.h"

#include "PMType.h"
#include "TMath.h"

////////////////////////////////////////////////////////////////////////////////
/// Default constructor
SPEResponse::SPEResponse() : IModel(), m_resp(PMType::Response::GAUSS)
{
   SetResponse(m_resp);
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor
SPEResponse::SPEResponse(PMType::Response resp) : IModel()
{
   SetResponse(resp);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor
SPEResponse::SPEResponse(const SPEResponse &other) : IModel(other), m_resp(other.m_resp)
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
         Error("gausexpfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));

      result = w * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) / (TMath::Sqrt(TMath::TwoPi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
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
         Error("gaus2expfunc", "Division by zero: sigma");

      Double_t gn = 0.5 * TMath::Erfc(-Q / (sqrt(2.0) * s));

      result = w1 * alpha1 * TMath::Exp(-x * alpha1) + w2 * alpha2 * TMath::Exp(-x * alpha2) +
               (1.0 - w1 - w2) / (TMath::Sqrt(TMath::TwoPi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
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

      result = w * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) * f * TMath::Power(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
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

      result = w1 * alpha1 * TMath::Exp(-x * alpha1) + w2 * alpha2 * TMath::Exp(-x * alpha2) +
               (1.0 - w1 - w2) * f * TMath::Power(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx);
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

      result = w * alpha * TMath::Exp(-x * alpha) +
               (1.0 - w) * f * TMath::Power(fx, kappa - 1.0) * TMath::Exp(-TMath::Power(fx, kappa));
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
               (1.0 - w) / (TMath::Sqrt(TMath::TwoPi()) * s * x) * TMath::Exp(-0.5 * arg * arg);
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

      Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));

      Double_t f = lambda * (1.0 + theta);
      Double_t fx = lambda * (1.0 + theta) * x;

      result = w * f * pow(fx, theta) / TMath::Gamma(1.0 + theta) * TMath::Exp(-fx) +
               (1.0 - w) / (TMath::Sqrt(TMath::TwoPi()) * s * gn) * TMath::Exp(-0.5 * arg * arg);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the gain
Double_t SPEResponse::Gain(const Double_t *pars) const
{
   Double_t result = 0.0;
   Double_t Q = pars[0];
   Double_t lambda = pars[0];
   Double_t alpha = pars[2];
   Double_t w = pars[3];
   switch (m_resp) {
   case PMType::Response::GAUSS:
   case PMType::Response::LOGNORMAL:
   case PMType::Response::GAUSS2EXP:      // Until I know what the calculation is
      result = w / alpha + (1.0 - w) * Q; // From PMTCalib example2

   case PMType::Response::GAMMA:
   case PMType::Response::WEIBULL:
   case PMType::Response::GAMMA2EXP: result = w / alpha + (1.0 - w) / lambda;

   case PMType::Response::TEST: result = 0.0;
   }
   Info("SPEResponse::GainError", "\n\tQ/theta\t: %.2e\nalpha\t: %.2e\nw\t: %.2e\ngainError\t: %.2e\n", Q, alpha, w,
        result);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the error on the gain const method
Double_t SPEResponse::GainError(const Double_t *pars, const Double_t *errs) const
{
   Double_t result = 0.0;
   Double_t Q = pars[0];
   Double_t lambda = pars[0];
   Double_t alpha = pars[2];
   Double_t w = pars[3];
   Double_t Qerr = errs[0];
   Double_t lambdaErr = errs[0];
   Double_t alphaErr = errs[2];
   Double_t wErr = errs[3];
   switch (m_resp) {
   case PMType::Response::GAUSS:
   case PMType::Response::LOGNORMAL:
   case PMType::Response::GAUSS2EXP: // Until I know what the calculation is
      result = TMath::Sqrt(TMath::Power((1.0 - w) * Qerr, 2.0) + TMath::Power((Q - alpha) * wErr, 2.0) +
                           TMath::Power(w / (alpha * alpha) * alphaErr, 2));
      break;

   case PMType::Response::GAMMA:
   case PMType::Response::WEIBULL:
   case PMType::Response::GAMMA2EXP:
      result = TMath::Sqrt(TMath::Power(wErr * (1 / alpha - 1 / lambda), 2.0) +
                           TMath::Power(-alphaErr * w / TMath::Power(alpha, 2.0), 2.0) +
                           TMath::Power(-lambdaErr * (1 - w) / TMath::Power(lambda, 2.0), 2.0));

      break;

   case PMType::Response::TEST: result = 0.0;
   }

   return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
///
Double_t SPEResponse::DoEvalPar(const Double_t x, const Double_t *pars) const
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
   default: Error("SPEResponse::DoEvalPar", "Invalid response type selected m_resp = %d", m_resp); break;
   }

   return retVal;
}
