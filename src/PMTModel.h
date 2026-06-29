#ifndef PMTMODEL_H
#define PMTMODEL_H

#include "IModel.h"
#include "PMType.h"

#include "TGraph.h"

class PMTModel : public IModel {
public:
   PMTModel() = default;
   PMTModel(const PMTModel &other);
   PMTModel(Int_t nbins, Double_t wbin, Double_t xmin, Double_t xmax, PMType::Model modtype);
   virtual PMTModel *Clone() const override { return new PMTModel(*this); }
   virtual ~PMTModel() = default;

   void SetModel(PMType::Model model);
   Double_t am(const Int_t m, const Double_t o) const;
   Double_t bm(const Int_t m, const Double_t o) const;
   Double_t SIMPLEGAUSS(Double_t x, const Double_t *pars) const;
   Double_t TRUNCGAUSS(Double_t x, const Double_t *pars) const;
   Double_t ANATRUNCG(Double_t x, const Double_t *pars) const;
   Double_t EXPTRUNCG(Double_t x, const Double_t *pars) const;
   
   TGraph *GetGraph();
   TGraph *GetGraphN(Int_t n);
   TGraph *GetGraphN2(Int_t n);
   
private:
   Int_t m_nlim = 10; ///< TODO:: figure out what this is
   Int_t m_nlim2 = 10; ///< 
   virtual Double_t DoEvalPar(Double_t x, const Double_t *p) const override;
   PMType::Model m_modtype{PMType::Model::EXPTRUNCG};
};

#endif
