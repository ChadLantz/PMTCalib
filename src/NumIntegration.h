#ifndef NUMINTEGRATION_H
#define NUMINTEGRATION_H

#include "IModel.h"
#include "PMType.h"
#include "SPEResponse.h"

#include "RtypesCore.h"
#include "TGraph.h"

class NumIntegration : public IModel {

public:
   NumIntegration() = default;
   NumIntegration(const NumIntegration &other);
   NumIntegration(Int_t _nbins, Double_t _wbin, Double_t _xmin, Double_t _xmax, PMType::Response _sper);
   virtual NumIntegration *Clone() const override { return new NumIntegration(*this); }
   virtual ~NumIntegration();

   virtual void SetParameters(const Double_t *pars) override;
   void CalculateValues(const Double_t *pars) const;
   TGraph *GetGraph();
   // These are probably unsafe
   virtual Double_t Gain(const Double_t *pars) const override { return m_resp.Gain(&pars[4]); }
   virtual Double_t GainError(const Double_t *pars, const Double_t *errs) const override
   {
      return m_resp.GainError(&pars[4], &errs[4]);
   }

private:
   virtual Double_t DoEvalPar(Double_t x, const Double_t *pars) const override;
   mutable UInt_t m_nCalls{0};
   mutable SPEResponse m_resp;
   TGraph *m_gr{nullptr};
};

#endif
