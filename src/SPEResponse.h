#ifndef SPERESPONSE_H
#define SPERESPONSE_H

#include "IModel.h"
#include "PMType.h"

#include "RtypesCore.h"

class SPEResponse : public IModel {
public:
   SPEResponse();
   SPEResponse(PMType::Response resp);
   SPEResponse(const SPEResponse &other);
   virtual SPEResponse *Clone() const override { return new SPEResponse(*this); }
   ~SPEResponse() = default;

   void SetResponse(PMType::Response resp);
   Double_t gausexpfunc(Double_t x, const Double_t *par) const;
   Double_t gaus2expfunc(Double_t x, const Double_t *par) const;
   Double_t gammaexpfunc(Double_t x, const Double_t *par) const;
   Double_t gamma2expfunc(Double_t x, const Double_t *par) const;
   Double_t weibullexpfunc(Double_t x, const Double_t *par) const;
   Double_t lognormalexpfunc(Double_t x, const Double_t *par) const;
   Double_t testfunc(Double_t x, const Double_t *par) const;
   virtual Double_t Gain(const Double_t *pars) const override;
   virtual Double_t GainError(const Double_t *pars, const Double_t *errs) const override;
   
private:
   virtual Double_t DoEvalPar(Double_t x, const Double_t *p) const override;
   PMType::Response m_resp;
};

#endif
