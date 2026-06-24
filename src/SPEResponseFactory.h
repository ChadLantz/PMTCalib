
#ifndef SPERESPONSEFactory_H
#define SPERESPONSEFactory_H

#include "TObject.h"
#include "TF1.h"

#include "PMType.h"

class SPEResponseFactory : public TObject {
public:
   TF1 Build(PMType::Response _spetype, Double_t _params[]);
   ClassDef(SPEResponseFactory, 1)
};

#endif
