
#ifndef PMTYPE_H
#define PMTYPE_H

#include <array>

namespace PMType {
enum Method{
   NumIntegration,
   DFTmethod,
   PMTModel,
};
static std::array<Method, 3> AllMethods{{NumIntegration, DFTmethod, PMTModel}};

enum Response {
   GAUSS,
   GAMMA,
   WEIBULL,
   LOGNORMAL,
   GAUSS2EXP,
   GAMMA2EXP,
   TEST
};
static std::array<Response, 7> AllResponseTypes{{GAUSS, GAMMA, WEIBULL, LOGNORMAL, GAUSS2EXP, GAMMA2EXP, TEST}};

enum Model {
   SIMPLEGAUSS,
   TRUNCGAUSS,
   ANATRUNCG,
   EXPTRUNCG
};
static std::array<Model, 4> AllModels{{SIMPLEGAUSS, TRUNCGAUSS, ANATRUNCG, EXPTRUNCG}};

static inline const char *ToString(Method method)
{
   switch (method) {
   case NumIntegration: return "NumIntegration";
   case DFTmethod: return "DFTmethod";
   case PMTModel: return "PMTModel";
   }
   return "Invalid Response";
}

static inline const char *ToString(Response response)
{
   switch (response) {
   case GAUSS: return "GAUSS";
   case GAMMA: return "GAMMA";
   case WEIBULL: return "WEIBULL";
   case LOGNORMAL: return "LOGNORMAL";
   case GAUSS2EXP: return "GAUSS2EXP";
   case GAMMA2EXP: return "GAMMA2EXP";
   case TEST: return "TEST";
   }
   return "Invalid Response";
}

static inline const char *ToString(Model model)
{
   switch (model) {
   case SIMPLEGAUSS: return "SIMPLEGAUSS";
   case TRUNCGAUSS: return "TRUNCGAUSS";
   case ANATRUNCG: return "ANATRUNCG";
   case EXPTRUNCG: return "EXPTRUNCG";
   }
   return "Invalid Model option";
}

} // namespace PMType

#endif
