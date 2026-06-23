
#ifndef SPERESPONSE_H
#define SPERESPONSE_H

#include "TObject.h"
#include "TF1.h"

#include "PMType.h"

class SPEResponse : public TObject {
private:
   PMType::Response spetype;

public:
   SPEResponse();
   SPEResponse(const SPEResponse &);
   SPEResponse(PMType::Response _spetype, Double_t _params[]);
   virtual ~SPEResponse();

   void SetParameters(Double_t *_pars);
   PMType::Response GetResponse() { return spetype; }
   Double_t GetValue(Double_t xx) { return spefunc->Eval(xx); }
   Double_t GenQ() { return spefunc->GetRandom(); }

   inline Double_t operator()(Double_t *xx, Double_t *pars){
      SetParameters(pars);
      return GetValue(xx[0]);
   };

   Int_t nparams;
   Double_t params[10];
   TF1 *spefunc; ///< This class is a wrapper for this TF1 and so will own it

   ClassDef(SPEResponse, 1)
};

#endif
