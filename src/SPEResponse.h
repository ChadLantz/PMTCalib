
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
   SPEResponse(PMType::Response _spetype, Double_t _params[]);
   virtual ~SPEResponse();
   
   void SetParams(Double_t _params[]);
   PMType::Response GetResponse(){ return spetype; }
   TF1 *GetFitFunc();
   Double_t GetValue(Double_t xx);
   Double_t GenQ();
   
   Int_t nparams;
   Double_t params[10] = {-1.0};
   TF1 *spefunc;

   ClassDef(SPEResponse, 1)
};

#endif
