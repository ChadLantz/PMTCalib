#include "PMTModel.h"

#include "IModel.h"
#include "PMType.h"
#include "TMath.h"
#include "Math/SpecFuncMathMore.h"

#include <vector>

PMTModel::PMTModel(Int_t nbins, Double_t wbin, Double_t xmin, Double_t xmax, PMType::Model modtype)
   : IModel(8, nbins, wbin, xmin, xmax), m_modtype(modtype)
{
   SetParNames({"Norm", "Q_{0}", "#sigma_{0}", "#mu", "Q", "#sigma", "#alpha", "w"});
}

PMTModel::PMTModel(const PMTModel &other) : IModel(other), m_modtype(other.m_modtype) {}

Double_t PMTModel::am(const Int_t m, const Double_t o) const
{
   Double_t a = 0.0;

   if (m == 1)
      a = +pow(o, 0.0);
   else if (m == 2)
      a = +2.0 * pow(o, 1.0);
   else if (m == 3)
      a = +2.0 * pow(o, 0.0) + 4.0 * pow(o, 2.0);
   else if (m == 4)
      a = +12.0 * pow(o, 1.0) + 8.0 * pow(o, 3.0);
   else if (m == 5)
      a = +12.0 * pow(o, 0.0) + 48.0 * pow(o, 2.0) + 16.0 * pow(o, 4.0);
   else if (m == 6)
      a = +120.0 * pow(o, 1.0) + 160.0 * pow(o, 3.0) + 32.0 * pow(o, 5.0);
   else if (m == 7)
      a = +120.0 * pow(o, 0.0) + 720.0 * pow(o, 2.0) + 480.0 * pow(o, 4.0) + 64.0 * pow(o, 6.0);
   else if (m == 8)
      a = +1680.0 * pow(o, 1.0) + 3360.0 * pow(o, 3.0) + 1344.0 * pow(o, 5.0) + 128.0 * pow(o, 7.0);
   else if (m == 9)
      a = +1680.0 * pow(o, 0.0) + 13440.0 * pow(o, 2.0) + 13440.0 * pow(o, 4.0) + 3584.0 * pow(o, 6.0) +
          256.0 * pow(o, 8.0);
   else if (m == 10)
      a = +30240.0 * pow(o, 1.0) + 80640.0 * pow(o, 3.0) + 48384.0 * pow(o, 5.0) + 9216.0 * pow(o, 7.0) +
          512.0 * pow(o, 9.0);
   else if (m == 11)
      a = +30240.0 * pow(o, 0.0) + 302400.0 * pow(o, 2.0) + 403200.0 * pow(o, 4.0) + 161280.0 * pow(o, 6.0) +
          23040.0 * pow(o, 8.0) + 1024.0 * pow(o, 10.0);
   else if (m == 12)
      a = +665280.0 * pow(o, 1.0) + 2217600.0 * pow(o, 3.0) + 1774080.0 * pow(o, 5.0) + 506880.0 * pow(o, 7.0) +
          56320.0 * pow(o, 9.0) + 2048.0 * pow(o, 11.0);
   else if (m == 13)
      a = +665280.0 * pow(o, 0.0) + 7983360.0 * pow(o, 2.0) + 13305600.0 * pow(o, 4.0) + 7096320.0 * pow(o, 6.0) +
          1520640.0 * pow(o, 8.0) + 135168.0 * pow(o, 10.0) + 4096.0 * pow(o, 12.0);
   else if (m == 14)
      a = +17297280.0 * pow(o, 1.0) + 69189120.0 * pow(o, 3.0) + 69189120.0 * pow(o, 5.0) + 26357760.0 * pow(o, 7.0) +
          4392960.0 * pow(o, 9.0) + 319488.0 * pow(o, 11.0) + 8192.0 * pow(o, 13.0);
   else if (m == 15)
      a = +17297280.0 * pow(o, 0.0) + 242161920.0 * pow(o, 2.0) + 484323840.0 * pow(o, 4.0) +
          322882560.0 * pow(o, 6.0) + 92252160.0 * pow(o, 8.0) + 12300288.0 * pow(o, 10.0) + 745472.0 * pow(o, 12.0) +
          16384.0 * pow(o, 14.0);
   else if (m == 16)
      a = +518918400.0 * pow(o, 1.0) + 2421619200.0 * pow(o, 3.0) + 2905943040.0 * pow(o, 5.0) +
          1383782400.0 * pow(o, 7.0) + 307507200.0 * pow(o, 9.0) + 33546240.0 * pow(o, 11.0) +
          1720320.0 * pow(o, 13.0) + 32768.0 * pow(o, 15.0);
   else if (m == 17)
      a = +518918400.0 * pow(o, 0.0) + 8302694400.0 * pow(o, 2.0) + 19372953600.0 * pow(o, 4.0) +
          15498362880.0 * pow(o, 6.0) + 5535129600.0 * pow(o, 8.0) + 984023040.0 * pow(o, 10.0) +
          89456640.0 * pow(o, 12.0) + 3932160.0 * pow(o, 14.0) + 65536.0 * pow(o, 16.0);
   else if (m == 18)
      a = +17643225088.0 * pow(o, 1.0) + 94097203200.0 * pow(o, 3.0) + 131736084480.0 * pow(o, 5.0) +
          75277762560.0 * pow(o, 7.0) + 20910489600.0 * pow(o, 9.0) + 3041525760.0 * pow(o, 11.0) +
          233963520.0 * pow(o, 13.0) + 8912896.0 * pow(o, 15.0) + 131072.0 * pow(o, 17.0);
   else if (m == 19)
      a = +17643225088.0 * pow(o, 0.0) + 317578051584.0 * pow(o, 2.0) + 846874869760.0 * pow(o, 4.0) +
          790416523264.0 * pow(o, 6.0) + 338749947904.0 * pow(o, 8.0) + 75277762560.0 * pow(o, 10.0) +
          9124577280.0 * pow(o, 12.0) + 601620480.0 * pow(o, 14.0) + 20054016.0 * pow(o, 16.0) +
          262144.0 * pow(o, 18.0);
   else if (m == 20)
      a = +670442586112.0 * pow(o, 1.0) + 4022655516672.0 * pow(o, 3.0) + 6436248616960.0 * pow(o, 5.0) +
          4290832760832.0 * pow(o, 7.0) + 1430277586944.0 * pow(o, 9.0) + 260050452480.0 * pow(o, 11.0) +
          26671841280.0 * pow(o, 13.0) + 1524105216.0 * pow(o, 15.0) + 44826624.0 * pow(o, 17.0) +
          524288.0 * pow(o, 19.0);
   return a;
}

Double_t PMTModel::bm(const Int_t m, const Double_t o) const
{
   Double_t b = 0.0;

   if (m == 1)
      b = 0.0;
   else if (m == 2)
      b = +pow(o, 0.0);
   else if (m == 3)
      b = +2.0 * pow(o, 1.0);
   else if (m == 4)
      b = +4.0 * pow(o, 0.0) + 4.0 * pow(o, 2.0);
   else if (m == 5)
      b = +20.0 * pow(o, 1.0) + 8.0 * pow(o, 3.0);
   else if (m == 6)
      b = +32.0 * pow(o, 0.0) + 72.0 * pow(o, 2.0) + 16.0 * pow(o, 4.0);
   else if (m == 7)
      b = +264.0 * pow(o, 1.0) + 224.0 * pow(o, 3.0) + 32.0 * pow(o, 5.0);
   else if (m == 8)
      b = +384.0 * pow(o, 0.0) + 1392.0 * pow(o, 2.0) + 640.0 * pow(o, 4.0) + 64.0 * pow(o, 6.0);
   else if (m == 9)
      b = +4464.0 * pow(o, 1.0) + 5920.0 * pow(o, 3.0) + 1728.0 * pow(o, 5.0) + 128.0 * pow(o, 7.0);
   else if (m == 10)
      b = +6144.0 * pow(o, 0.0) + 31200.0 * pow(o, 2.0) + 22080.0 * pow(o, 4.0) + 4480.0 * pow(o, 6.0) +
          256.0 * pow(o, 8.0);
   else if (m == 11)
      b = +92640.0 * pow(o, 1.0) + 168960.0 * pow(o, 3.0) + 75264.0 * pow(o, 5.0) + 11264.0 * pow(o, 7.0) +
          512.0 * pow(o, 9.0);
   else if (m == 12)
      b = +122880.0 * pow(o, 0.0) + 809280.0 * pow(o, 2.0) + 779520.0 * pow(o, 4.0) + 240128.0 * pow(o, 6.0) +
          27648.0 * pow(o, 8.0) + 1024.0 * pow(o, 10.0);
   else if (m == 13)
      b = +2283840.0 * pow(o, 1.0) + 5335680.0 * pow(o, 3.0) + 3214848.0 * pow(o, 5.0) + 728064.0 * pow(o, 7.0) +
          66560.0 * pow(o, 9.0) + 2048.0 * pow(o, 11.0);
   else if (m == 14)
      b = +2949120.0 * pow(o, 0.0) + 23990400.0 * pow(o, 2.0) + 29379840.0 * pow(o, 4.0) + 12192768.0 * pow(o, 6.0) +
          2119680.0 * pow(o, 8.0) + 157696.0 * pow(o, 10.0) + 4096.0 * pow(o, 12.0);
   else if (m == 15)
      b = +65278080.0 * pow(o, 1.0) + 186708480.0 * pow(o, 3.0) + 142345728.0 * pow(o, 5.0) + 43315200.0 * pow(o, 7.0) +
          5969920.0 * pow(o, 9.0) + 368640.0 * pow(o, 11.0) + 8192.0 * pow(o, 13.0);
   else if (m == 16)
      b = +82575360.0 * pow(o, 0.0) + 802287360.0 * pow(o, 2.0) + 1196052480.0 * pow(o, 4.0) +
          626088960.0 * pow(o, 6.0) + 145981440.0 * pow(o, 8.0) + 16355328.0 * pow(o, 10.0) + 851968.0 * pow(o, 12.0) +
          16384.0 * pow(o, 14.0);
   else if (m == 17)
      b = +2123493120.0 * pow(o, 1.0) + 7205829120.0 * pow(o, 3.0) + 6662476800.0 * pow(o, 5.0) +
          2551633920.0 * pow(o, 7.0) + 471060480.0 * pow(o, 9.0) + 43769856.0 * pow(o, 11.0) +
          1949696.0 * pow(o, 13.0) + 32768.0 * pow(o, 15.0);
   else if (m == 18)
      b = +2642411520.0 * pow(o, 0.0) + 29920182272.0 * pow(o, 2.0) + 52685336576.0 * pow(o, 4.0) +
          33359800320.0 * pow(o, 6.0) + 9774673920.0 * pow(o, 8.0) + 1465491456.0 * pow(o, 10.0) +
          114802688.0 * pow(o, 12.0) + 4423680.0 * pow(o, 14.0) + 65536.0 * pow(o, 16.0);
   else if (m == 19)
      b = +77483589632.0 * pow(o, 1.0) + 304838541312.0 * pow(o, 3.0) + 331894882304.0 * pow(o, 5.0) +
          153475153920.0 * pow(o, 7.0) + 35565404160.0 * pow(o, 9.0) + 4419158016.0 * pow(o, 11.0) +
          295895040.0 * pow(o, 13.0) + 9961472.0 * pow(o, 15.0) + 131072.0 * pow(o, 17.0);
   else if (m == 20)
      b = +95126814720.0 * pow(o, 0.0) + 1232093708288.0 * pow(o, 2.0) + 2506349346816.0 * pow(o, 4.0) +
          1864742600704.0 * pow(o, 6.0) + 658838585344.0 * pow(o, 8.0) + 123888500736.0 * pow(o, 10.0) +
          12971212800.0 * pow(o, 12.0) + 751042560.0 * pow(o, 14.0) + 22282240.0 * pow(o, 16.0) +
          262144.0 * pow(o, 18.0);

   return b;
}

Double_t PMTModel::DoEvalPar(Double_t x, const Double_t *p) const
{
   switch (m_modtype) {
   case PMType::Model::SIMPLEGAUSS: return SIMPLEGAUSS(x, p);
   case PMType::Model::TRUNCGAUSS: return TRUNCGAUSS(x, p);
   case PMType::Model::ANATRUNCG: return ANATRUNCG(x, p);
   case PMType::Model::EXPTRUNCG: return EXPTRUNCG(x, p);
   }
   return 0.0;
}

Double_t PMTModel::SIMPLEGAUSS(Double_t x, const Double_t *pars) const
{
   Double_t result = 0.0;
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];
   Double_t Q = pars[4];
   Double_t s = pars[5];
   Double_t alpha = pars[6];
   Double_t w = pars[7];

   result += TMath::Poisson(0.0, mu) * TMath::Gaus(x, Q0, s0, kTRUE);

   Double_t omega = (Q0 + alpha * TMath::Sq(s0) - x) * 0.5 * TMath::Sqrt2() / s0;
   Double_t SR1 = w * alpha * 0.5 * TMath::Exp(-alpha * (x - Q0) + TMath::Sq(alpha * s0) * 0.5) * TMath::Erfc(omega);
   SR1 += (1.0 - w) * TMath::Gaus(x, Q0 + Q, TMath::Hypot(s0, s), kTRUE);
   result += TMath::Poisson(1.0, mu) * SR1;

   Double_t Qs = w / alpha + (1.0 - w) * Q;
   Double_t ss2 = w / TMath::Sq(alpha) + (1 - w) * TMath::Sq(s) + w * (1.0 - w) * TMath::Sq(Q - 1.0 / alpha);

   for (Int_t n = 2; n < 25; n++) {
      result += TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(TMath::Sq(s0) + n * ss2), kTRUE);
   }

   return result * Norm * m_wBin;
}

Double_t PMTModel::TRUNCGAUSS(Double_t x, const Double_t *pars) const
{
   Double_t result = 0.0;
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];
   Double_t Q = pars[4];
   Double_t s = pars[5];
   Double_t alpha = pars[6];
   Double_t w = pars[7];

   result += TMath::Poisson(0.0, mu) * TMath::Gaus(x, Q0 + Q, TMath::Hypot(s0, s), kTRUE);

   Double_t s1 = TMath::Hypot(s0, s);
   Double_t omega = (Q0 + alpha * TMath::Sq(s0) - x) * 0.5 * TMath::Sqrt2() / s0;
   Double_t SR1 = w * alpha / 2.0 * TMath::Exp(-alpha * (x - Q0) + TMath::Sq(alpha * s0) / 2.0) * TMath::Erfc(omega);

   Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
   Double_t A = (Q0 - x) * TMath::Sq(s) - Q * TMath::Sq(s0);
   Double_t B = TMath::Sqrt2() * s0 * s * s1;
   SR1 += (1.0 - w) * 0.5 / gn * TMath::Gaus(x, Q0 + Q, s1, kTRUE) * TMath::Erfc(A / B);
   result += TMath::Poisson(1.0, mu) * SR1;

   Double_t k = s / (gn * TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-TMath::Sq(Q) / (2.0 * TMath::Sq(s)));
   Double_t Qf = Q + k;
   Double_t sf2 = TMath::Sq(s) - (Q + k) * k;

   Double_t Qs = w / alpha + (1.0 - w) * Qf;
   Double_t ss2 = w / TMath::Sq(alpha) + (1 - w) * sf2 + w * (1.0 - w) * TMath::Sq(Qf - 1.0 / alpha);

   for (Int_t n = 2; n < 25; n++) {
      result += TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(TMath::Sq(s0) + n * ss2), kTRUE);
   }

   return result * Norm * m_wBin;
}

Double_t PMTModel::ANATRUNCG(Double_t x, const Double_t *pars) const
{
   Double_t result = 0.0;
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];
   Double_t Q = pars[4];
   Double_t s = pars[5];
   Double_t alpha = pars[6];
   Double_t w = pars[7];

   result += TMath::Poisson(0, mu) * TMath::Gaus(x, Q0, s0, kTRUE); // n = 0

   Double_t omega0 = (x - Q0 - alpha * TMath::Sq(s0)) * 0.5 * TMath::Sqrt2() / s0;
   Double_t SR1 = w * alpha * 0.5 * TMath::Exp(TMath::Sq(alpha * s0) * 0.5 - alpha * (x - Q0)) * TMath::Erfc(-omega0);

   Double_t s1 = TMath::Hypot(s0, s);
   Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
   Double_t A = (Q0 - x) * TMath::Sq(s) - Q * TMath::Sq(s0);
   Double_t B = TMath::Sqrt2() * s0 * s * s1;
   SR1 += (1.0 - w) * 0.5 / gn * TMath::Gaus(x, Q0 + Q, s1, kTRUE) * TMath::Erfc(A / B);
   SR1 *= TMath::Poisson(1, mu);
   result += SR1; // n = 1

   Double_t k = s / gn / TMath::Sqrt(TMath::TwoPi()) * TMath::Exp(-0.5 * TMath::Sq(Q) / TMath::Sq(s));
   Double_t Qg = Q + k;
   Double_t sg2 = TMath::Sq(s) - Qg * k;

   for (Int_t n = 2; n < m_nlim; n++) {
      Double_t SRn = TMath::Power(1.0 - w, n) * TMath::Gaus(x,  Q0 + n * Qg, TMath::Sqrt(TMath::Sq(s0) + n * sg2), kTRUE);

      for (Int_t m = 1; m <= n; m++) {
         Double_t Qmn = Q0 + (n - m) * Qg;
         Double_t smn2 = TMath::Sq(s0) + (n - m) * sg2;
         Double_t smn = TMath::Sqrt(smn2);

         Double_t cmn = alpha * TMath::Power(alpha * smn * TMath::Sqrt2(), m - 1.0) / TMath::Factorial(m - 1.0);

         Double_t psi = (x - Qmn) * 0.5 * TMath::Sqrt2() / smn;
         Double_t psi2 = TMath::Sq(psi);
         Double_t omega = (x - Qmn - alpha * TMath::Sq(smn)) * 0.5 * TMath::Sqrt2() / smn;
         Double_t omega2 = TMath::Sq(omega);

         Double_t A1m = m * 0.5;
         Double_t A2m = (m + 1.0) * 0.5;

         Double_t Imn = 0.0;
         Double_t hi_limit = 25.0;

         if (omega >= hi_limit) {
            Imn = TMath::Exp(omega2 - psi2 + (m - 1.0) * TMath::Log(omega));

         } else if (omega < hi_limit && omega >= 0.0) {
            Double_t t1 = TMath::Gamma(A1m) * ROOT::Math::conf_hyperg(0.5 - A1m, 0.5, -omega2);
            Double_t t2 =
               2.0 * omega * TMath::Gamma(A2m) * ROOT::Math::conf_hyperg(1.5 - A2m, 1.5, -omega2);
            Imn = 0.5 / TMath::Sqrt(TMath::Pi()) * (t1 + t2) * TMath::Exp(omega2 - psi2);

         } else if (omega < 0.0) {
            Double_t t3 = 0.5 * TMath::InvPi() * TMath::Gamma(A1m) * TMath::Gamma(A2m);
            Imn = t3 * ROOT::Math::conf_hypergU(A1m, 0.5, omega2) * TMath::Exp(-psi2);
         }

         Double_t hmnB = cmn * Imn;
         SRn += TMath::Binomial(n, m) * TMath::Power(w, m) * TMath::Power(1.0 - w, n - m) * hmnB;
      }

      result += SRn * TMath::Poisson(n, mu); // n = 2-nlim
   }

   Double_t Qs = w / alpha + (1.0 - w) * Qg;
   Double_t ss2 = w / TMath::Sq(alpha) + (1.0 - w) * sg2 + w * (1.0 - w) * TMath::Sq(Qg - 1.0 / alpha);

   for (Int_t n = m_nlim; n < 45; n++) {
      result += TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(TMath::Sq(s0) + n * ss2), kTRUE); // n >= nlim
   }

   return result * Norm * m_wBin;
}

Double_t PMTModel::EXPTRUNCG(Double_t x, const Double_t *pars) const
{
   Double_t result = 0.0;
   Double_t Norm = pars[0];
   Double_t Q0 = pars[1];
   Double_t s0 = pars[2];
   Double_t mu = pars[3];
   Double_t Q = pars[4];
   Double_t s = pars[5];
   Double_t alpha = pars[6];
   Double_t w = pars[7];

   result += TMath::Poisson(0, mu) * TMath::Gaus(x, Q0, s0, kTRUE); // n = 0

   Double_t omega0 = (x - Q0 - alpha * TMath::Sq(s0)) * 0.5 * TMath::Sqrt2() / s0;
   Double_t SR1 = 0.25 * w * alpha * TMath::Exp(TMath::Sq(alpha * s0) - alpha * (x - Q0)) * TMath::Erfc(-omega0);

   Double_t Q1 = Q0 + Q;
   Double_t s1 = TMath::Hypot(s0, s);

   Double_t arg1 = 0.0;
   if (s1 != 0.0)
      arg1 = (x - Q1) / s1;

   Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
   Double_t A = (Q0 - x) * TMath::Sq(s) - Q * TMath::Sq(s0);
   Double_t B = TMath::Sqrt2() * s0 * s * s1;
   SR1 += (1.0 - w) * 0.5 / gn * TMath::Gaus(x, Q0 + Q, s1, kTRUE) * TMath::Erfc(A / B);
   SR1 *= TMath::Poisson(1, mu);
   result += SR1; // n = 1

   Double_t k = s / (gn * TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-TMath::Sq(Q) / (2.0 * TMath::Sq(s)));
   Double_t Qg = Q + k;
   Double_t sg2 = TMath::Sq(s) - Qg * k;

   for (Int_t n = 2; n < m_nlim2; n++) {
      Double_t SRn = TMath::Power(1.0 - w, n) * TMath::Gaus(x, Q0 + n * Qg, TMath::Sqrt(TMath::Sq(s0) + n * sg2), kTRUE);

      for (Int_t m = 1; m <= n; m++) {
         Double_t Qmn = Q0 + (n - m) * Qg;
         Double_t smn2 = TMath::Sq(s0) + (n - m) * sg2;
         Double_t smn = TMath::Sqrt(smn2);

         Double_t cmn = alpha * TMath::Power(alpha * smn * TMath::Sqrt2(), m - 1.0) / TMath::Factorial(m - 1.0);

         Double_t psi = (x - Qmn) * 0.5 * TMath::Sqrt2() / smn;
         Double_t psi2 = TMath::Sq(psi);
         Double_t omega = (x - Qmn - alpha * TMath::Sq(smn)) * 0.5 * TMath::Sqrt2() / smn;
         Double_t omega2 = TMath::Sq(omega);
         Double_t Imn = TMath::Power(2.0, -m) * (am(m, omega) * TMath::Exp(omega2 - psi2) * TMath::Erfc(-omega) +
                                                 2.0 / TMath::Sqrt(TMath::Pi()) * bm(m, omega) * TMath::Exp(-psi2));

         Double_t hmnB = cmn * Imn;
         SRn += TMath::Binomial(n, m) * pow(w, m) * pow(1.0 - w, n - m) * hmnB;
      }

      SRn *= TMath::Poisson(n, mu);
      result += SRn; // n = 2-nlim
   }

   Double_t Qs = w / alpha + (1.0 - w) * Qg;
   Double_t ss2 = w / TMath::Sq(alpha) + (1.0 - w) * sg2 + w * (1.0 - w) * TMath::Sq(Qg - 1.0 / alpha);

   for (Int_t n = m_nlim2; n < 45; n++) {
      result += TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(TMath::Sq(s0) + n * ss2), kTRUE);
   }

   return result * Norm * m_wBin;
}

TGraph *PMTModel::GetGraph()
{
   TGraph *_gr = new TGraph(m_nBins);
   for (Int_t i = 0; i < m_nBins; i++) {
      Double_t x = m_xMin + m_step / 2 + i * m_step;
      Double_t y = DoEvalPar(x, Parameters());
      if (y < 1.0e-10)
         y = 1.e-3;

      _gr->SetPoint(i, x, y);
   }

   _gr->SetLineWidth(2);
   int cc = kAzure + 1;
   _gr->SetLineColor(cc);
   _gr->SetMarkerColor(cc);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

TGraph *PMTModel::GetGraphN(Int_t n)
{
   Double_t Norm = m_parameters[0];
   Double_t Q0 = m_parameters[1];
   Double_t s0 = m_parameters[2];
   Double_t mu = m_parameters[3];
   Double_t Q = m_parameters[4];
   Double_t s = m_parameters[5];
   Double_t alpha = m_parameters[6];
   Double_t w = m_parameters[7];

   TGraph *_gr = new TGraph(m_nBins);

   for (Int_t i = 0; i < m_nBins; i++) {
      Double_t x = m_xMin + m_step / 2 + i * m_step;
      Double_t y = 0.0;

      if (n == 0) {
         Double_t SR0 = TMath::Poisson(0, mu) * TMath::Gaus(x, Q0, s0, kTRUE);
         ;
         y = Norm * m_wBin * SR0; // 0
      }

      if (n == 1) {
         Double_t omega0 = 0.5 * TMath::Sqrt2() * (x - Q0 - alpha * pow(s0, 2.0)) / s0;
         Double_t SR1 =
            0.5 * w * alpha * TMath::Exp(-alpha * (x - Q0) + 0.5 * pow(alpha * s0, 2.0)) * TMath::Erfc(-omega0);

         Double_t s1 = TMath::Hypot(s0, s);
         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t A = (Q0 - x) * pow(s, 2.0) - Q * pow(s0, 2.0);
         Double_t B = TMath::Sqrt2() * s0 * s * s1;
         SR1 += (1.0 - w) * 0.5 / gn * TMath::Gaus(x, Q0 + Q, s1, kTRUE) * TMath::Erfc(A / B);
         SR1 *= TMath::Poisson(1, mu);
         y += Norm * m_wBin * SR1; // 1
      }

      if (n >= 2 && n < m_nlim) {
         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t k = s / gn / TMath::Sqrt(TMath::TwoPi()) * TMath::Exp(-pow(Q, 2.0) / (2.0 * pow(s, 2.0)));
         Double_t Qg = Q + k;
         Double_t sg2 = pow(s, 2.0) - Qg * k;
         Double_t SRn = pow(1.0 - w, n) * TMath::Gaus(x, Q0 + n * Qg, TMath::Sqrt(pow(s0, 2.0) + n * sg2), kTRUE);

         for (Int_t m = 1; m <= n; m++) {
            Double_t Qmn = Q0 + (n - m) * Qg;
            Double_t smn2 = pow(s0, 2.0) + (n - m) * sg2;
            Double_t smn = TMath::Sqrt(smn2);

            Double_t cmn = alpha * pow(alpha * smn * TMath::Sqrt2(), m - 1.0) / TMath::Factorial(m - 1.0);

            Double_t psi = (x - Qmn) * 0.5 * TMath::Sqrt2() / smn;
            Double_t psi2 = pow(psi, 2.0);
            Double_t omega = (x - Qmn - alpha * pow(smn, 2.0)) * 0.5 * TMath::Sqrt2() / smn;
            Double_t omega2 = pow(omega, 2.0);

            Double_t A1m = m / 2.0;
            Double_t A2m = (m + 1.0) / 2.0;

            Double_t Imn = 0.0;
            Double_t hi_limit = 25.0;

            if (omega >= hi_limit) {
               Imn = TMath::Exp(omega2 - psi2 + (m - 1.0) * TMath::Log(omega));

            } else if (omega < hi_limit && omega >= 0.0) {
               Double_t t1 = TMath::Gamma(A1m) * ROOT::Math::conf_hyperg(0.5 - A1m, 0.5, -omega2);
               Double_t t2 = 2.0 * omega * TMath::Gamma(A2m) * ROOT::Math::conf_hyperg(1.5 - A2m, 1.5, -omega2);
               Imn = 0.5 / TMath::Sqrt(TMath::Pi()) * (t1 + t2) * TMath::Exp(omega2 - psi2);

            } else if (omega < 0.0) {
               Double_t t3 = TMath::Gamma(A1m) * TMath::Gamma(A2m) * 0.5 * TMath::InvPi();
               Imn = t3 * ROOT::Math::conf_hypergU(A1m, 0.5, omega2) * TMath::Exp(-psi2);
            }

            Double_t hmnB = cmn * Imn;
            SRn += TMath::Binomial(n, m) * pow(w, m) * pow(1.0 - w, n - m) * hmnB;
         }

         SRn *= TMath::Poisson(n, mu);
         y += Norm * m_wBin * SRn; // n= 2-m_nlim
      }

      if (n >= m_nlim) {
         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t k = s / gn / TMath::Sqrt(TMath::TwoPi()) * TMath::Exp(-pow(Q, 2.0) / (2.0 * pow(s, 2.0)));
         Double_t Qg = Q + k;
         Double_t sg2 = pow(s, 2.0) - Qg * k;

         Double_t Qs = w / alpha + (1.0 - w) * Qg;
         Double_t ss2 = w / pow(alpha, 2.0) + (1 - w) * sg2 + w * (1.0 - w) * pow(Qg - 1.0 / alpha, 2.0);
         Double_t SRn = TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(pow(s0, 2.0) + n * ss2), kTRUE);
         y += Norm * m_wBin * SRn; // n >= nlim
      }

      if (y < 1.0e-10)
         y = 1.e-3;

      _gr->SetPoint(i, x, y);
   }

   _gr->SetLineWidth(2);
   _gr->SetLineStyle(3);
   _gr->SetLineColor(kBlack);
   _gr->SetMarkerColor(kBlack);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

TGraph *PMTModel::GetGraphN2(Int_t n)
{
   Double_t Norm = m_parameters[0];
   Double_t Q0 = m_parameters[1];
   Double_t s0 = m_parameters[2];
   Double_t mu = m_parameters[3];
   Double_t Q = m_parameters[4];
   Double_t s = m_parameters[5];
   Double_t alpha = m_parameters[6];
   Double_t w = m_parameters[7];

   TGraph *_gr = new TGraph(m_nBins);

   for (Int_t i = 0; i < m_nBins; i++) {
      Double_t x = m_xMin + m_step / 2 + i * m_step;
      Double_t y = 0.0;

      if (n == 0) {
         Double_t arg = 0.0;
         if (s0 != 0.0)
            arg = (x - Q0) / s0;
         else
            Error("GetGraphN2", "Division by zero: sigma0");

         Double_t SR0 = 1.0 / (TMath::Sqrt(TMath::TwoPi()) * s0) * TMath::Exp(-0.5 * arg * arg);
         SR0 *= TMath::Poisson(0, mu);
         y = Norm * m_wBin * SR0; // 0
      }

      if (n == 1) {
         Double_t omega0 = (x - Q0 - alpha * TMath::Sq(s0)) * 0.5 * TMath::Sqrt2() / s0;
         Double_t SR1 =
            w * alpha / 2.0 * TMath::Exp(-alpha * (x - Q0) + TMath::Sq(alpha * s0) / 2.0) * TMath::Erfc(-omega0);

         Double_t Q1 = Q0 + Q;
         Double_t s12 = TMath::Sq(s0) + TMath::Sq(s);
         Double_t s1 = TMath::Sqrt(s12);

         Double_t arg1 = 0.0;
         if (s1 != 0.0)
            arg1 = (x - Q1) / s1;
         else
            Error("F3", "Division by zero: ");

         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t A = (Q0 - x) * TMath::Sq(s) - Q * TMath::Sq(s0);
         Double_t B = TMath::Sqrt2() * s0 * s * s1;
         SR1 +=
            (1.0 - w) / 2.0 / gn / (TMath::Sqrt(TMath::TwoPi()) * s1) * TMath::Exp(-0.5 * arg1 * arg1) * TMath::Erfc(A / B);
         SR1 *= TMath::Poisson(1, mu);
         y += Norm * m_wBin * SR1; // 1
      }

      if (n >= 2 && n < m_nlim2) {
         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t k = s / gn / TMath::Sqrt(TMath::TwoPi()) * TMath::Exp(-TMath::Sq(Q) / (2.0 * TMath::Sq(s)));
         Double_t Qg = Q + k;
         Double_t sg2 = TMath::Sq(s) - Qg * k;

         Double_t SRn = 0.0;

         Double_t Qn = Q0 + n * Qg;
         Double_t sn2 = TMath::Sq(s0) + n * sg2;
         Double_t sn = TMath::Sqrt(sn2);

         Double_t argn = 0.0;
         if (sn != 0.0)
            argn = (x - Qn) / sn;
         else
            Error("GetGraphN2", "Division by zero: sigmaN");
         Double_t gnB = 1.0 / (TMath::Sqrt(TMath::TwoPi()) * sn) * TMath::Exp(-0.5 * argn * argn);
         SRn += TMath::Power(1.0 - w, n) * gnB;

         for (Int_t m = 1; m <= n; m++) {
            Double_t Qmn = Q0 + (n - m) * Qg;
            Double_t smn2 = TMath::Sq(s0) + (n - m) * sg2;
            Double_t smn = TMath::Sqrt(smn2);

            Double_t cmn = alpha * TMath::Power(alpha * smn * TMath::Sqrt2(), m - 1.0) / TMath::Factorial(m - 1.0);

            Double_t psi = (x - Qmn) * 0.5 * TMath::Sqrt2() / smn;
            Double_t psi2 = TMath::Sq(psi);
            Double_t omega = (x - Qmn - alpha * TMath::Sq(smn)) * 0.5 * TMath::Sqrt2() / smn;
            Double_t omega2 = TMath::Sq(omega);
            Double_t Imn = TMath::Power(2.0, -m) * (am(m, omega) * TMath::Exp(omega2 - psi2) * TMath::Erfc(-omega) +
                                                    2.0 / TMath::Sqrt(TMath::Pi()) * bm(m, omega) * TMath::Exp(-psi2));

            Double_t hmnB = cmn * Imn;
            SRn += TMath::Binomial(n, m) * TMath::Power(w, m) * TMath::Power(1.0 - w, n - m) * hmnB;
         }

         SRn *= TMath::Poisson(n, mu);
         y += Norm * m_wBin * SRn; // n= 2-nlim
      }

      if (n >= m_nlim2) {
         Double_t gn = 0.5 * TMath::Erfc(-Q / (TMath::Sqrt2() * s));
         Double_t k = s / gn / TMath::Sqrt(TMath::TwoPi()) * TMath::Exp(-0.5 *TMath::Sq(Q) / TMath::Sq(s));
         Double_t Qg = Q + k;
         Double_t sg2 = TMath::Sq(s) - Qg * k;

         Double_t Qs = w / alpha + (1.0 - w) * Qg;
         Double_t ss2 = w / TMath::Sq(alpha) + (1 - w) * sg2 + w * (1.0 - w) * TMath::Sq(Qg - 1.0 / alpha);
         Double_t SRn = TMath::Poisson(n, mu) * TMath::Gaus(x, Q0 + n * Qs, TMath::Sqrt(TMath::Sq(s0) + n * ss2), kTRUE);
         y += Norm * m_wBin * SRn; // n >= nlim
      }

      if (y < 1.0e-10)
         y = 1.e-3;

      _gr->SetPoint(i, x, y);
   }

   _gr->SetLineWidth(2);
   _gr->SetLineStyle(3);
   _gr->SetLineColor(kBlack);
   _gr->SetMarkerColor(kBlack);
   _gr->SetMarkerSize(0.1);

   return _gr;
}

Double_t PMTModel::Gain(const Double_t *pars) const
{
   Double_t Q = pars[4];
   Double_t alpha = pars[6];
   Double_t w = pars[7];
   return w / alpha + (1.0 - w) * Q; // From PMTCalib example2
}

Double_t PMTModel::GainError(const Double_t *pars, const Double_t *errs) const
{
   Double_t Q = pars[4];
   Double_t Qerr = errs[4];
   Double_t alpha = pars[6];
   Double_t alphaErr = errs[6];
   Double_t w = pars[7];
   Double_t wErr = errs[7];
   return TMath::Sqrt(TMath::Power((1.0 - w) * Qerr, 2) + TMath::Power((Q - alpha) * wErr, 2) +
                      TMath::Power(w / (alpha * alpha) * alphaErr, 2));
}
