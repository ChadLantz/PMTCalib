#ifndef IMODEL_H
#define IMODEL_H

#include "RtypesCore.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"
#include "Math/IParamFunction.h"
#include "TError.h"

#include <vector>

class IModel : public ROOT::Math::IParametricFunctionOneDim {
public:
   using BaseFunc = IBaseFunctionOneDim;
   virtual ~IModel() {}

   // Required interface
   virtual BaseFunc *Clone() const override = 0;
   virtual Double_t Gain(const Double_t *pars) const = 0;
   virtual Double_t GainError(const Double_t *pars, const Double_t *errs) const = 0;
   Double_t Gain() const { return Gain(Parameters()); } ;
   Double_t GainError() const { return GainError(Parameters(), Errors()); };


   ////////////////////////////////////////////////////////////////////////////////
   /// Set the parameter settings
   virtual void SetChiSquare(Double_t chi2) { m_chi2 = chi2;}

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the parameter settings
   virtual void SetNDF(Double_t ndf) { m_ndf = ndf;}

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the parameter settings
   virtual void SetParamsSettings(std::vector<ROOT::Fit::ParameterSettings> &settings) { m_parSettings = settings; }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the parameter errors
   virtual void SetParErrors(const Double_t *errs)
   {
      assert(errs);
      if (errs) {
         for (UInt_t i; i < NPar(); ++i) {
            m_parErrors[i] = errs[i];
            m_parSettings[i].SetStepSize(errs[i]);
         }
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the parameters with the given array (non-const)
   virtual void SetParameters(const Double_t *pars) override
   {
      assert(pars);
      if (pars) {
         for (UInt_t i; i < NPar(); ++i) {
            m_parameters[i] = pars[i];
            m_parSettings[i].SetValue(pars[i]);
         }
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the number of parameters
   virtual void SetNpar(UInt_t n)
   {
      if (n > m_maxParams) {
         Error("SetNpar", "%d > maximum number of allowed parameters %d", n, m_maxParams);
         return;
      }
      m_parSettings.resize(n);
      m_parameters.resize(n);
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set parameter names (vector overload)
   virtual void SetParNames(std::vector<std::string> names)
   {
      if (names.size() <= NPar()) {
         for (UInt_t ipar = 0; ipar < names.size(); ++ipar) {
            m_parSettings.at(ipar).SetName(names.at(ipar));
         }
      } else if (names.size() > NPar() && names.size() < m_maxParams) {
         Warning("SetParNames", "names vector is longer than the number of parameters (%ld > %d)", names.size(),
                 NPar());
      } else if (names.size() >= m_maxParams) {
         Error("SetParNames", "names vector is longer than the maximum number of parameters (%ld > %d)", names.size(),
               m_maxParams);
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the number of bins
   virtual void SetNbins(UInt_t nbins)
   {
      m_nBins = nbins;
      m_step = ((m_xMax - m_xMin) / Double_t(m_nBins));
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set bin width
   virtual void SetBinWidth(Double_t wbin) { m_wBin = wbin; }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the function range
   virtual void SetRange(Double_t xmin, Double_t xmax)
   {
      m_xMin = xmin;
      m_xMax = xmax;
      // m_step = ((xmax - xmin) / Double_t(m_nBins));
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Set the result from the fit. Similar to TF1::SetFitResult, but I keep fewer
   /// tools around
   void SetFitResult(const ROOT::Fit::FitResult &result)
   {
      if (result.IsEmpty()) {
         Warning("SetFitResult", "Empty Fit result - nothing is set in TF1");
         return;
      }
      if (NPar() != result.NPar()) {
         Error("SetFitResult", "Invalid Fit result passed - number of parameter is %d , different than IModel::GetNpar() = %d", NPar(), result.NPar());
         return;
      }
   
      // Call set parameters so the models can override
      SetParameters(result.Parameters().data());
      SetParErrors(result.Errors().data());
      SetChiSquare(result.Chi2());
      SetNDF(result.Ndf());
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Return the number of parameters
   virtual UInt_t NPar() const override { return m_parSettings.size(); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Return the pointer to an array of parameter values
   virtual const Double_t *Parameters() const override { return m_parameters.data(); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Return the pointer to an array of parameter values
   virtual const Double_t *Errors() const { return m_parErrors.data(); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get the vector of parameter settings (non-const method)
   std::vector<ROOT::Fit::ParameterSettings> ParamsSettings() { return m_parSettings; }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get the parameter settings for the i-th parameter (non-const method)
   virtual ROOT::Fit::ParameterSettings &ParSettings(UInt_t i) { return m_parSettings.at(i); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get the parameter settings for the i-th parameter (const method)
   virtual const ROOT::Fit::ParameterSettings &ParSettings(UInt_t i) const { return m_parSettings.at(i); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get the name of the i-th parameter
   virtual std::string ParameterName(unsigned int i) const override { return m_parSettings.at(i).Name(); }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get the parameter settings for the i-th parameter (const method)
   virtual void GetRange(Double_t &xmin, Double_t &xmax)
   {
      xmin = m_xMin;
      xmax = m_xMax;
   }

   ////////////////////////////////////////////////////////////////////////////////
   /// Get Chi2
   virtual Double_t GetChiSquare() const { return m_chi2;}

   ////////////////////////////////////////////////////////////////////////////////
   /// Get NDF
   virtual Double_t GetNDF() const { return m_ndf;}

protected:
   // Constructor for shared metadata
   static const UInt_t m_maxParams{11};
   UInt_t m_nBins{0};
   Double_t m_wBin{0.0};
   Double_t m_step{0.0};
   Double_t m_xMin{0.0};
   Double_t m_xMax{0.0};
   Double_t m_chi2{0.0};
   Double_t m_ndf{0.0};

   /// A note on these two. IParametricFunctionOneDim demands a DoEval const method, but DFTmethod and NumIntegration
   /// need to keep an internal buffer
   mutable std::vector<Double_t> m_parameters;       ///< Parameter values. Mutable for DFTmethod and NumIntegration
   std::vector<Double_t> m_parErrors;                ///< Parameter errors
   mutable std::vector<ROOT::Fit::ParameterSettings> ///< Parameter settings. Mutable for DFTmethod and NumIntegration
      m_parSettings; ///< Could have been an array, but I'm doing this for compatibility with FitConfig

   IModel() : m_parameters(m_maxParams) {};

   IModel(UInt_t npar, UInt_t nbins, Double_t wbin, Double_t xmin, Double_t xmax)
      : m_nBins(nbins),
        m_wBin(wbin),
        m_step((xmax - xmin) / Double_t(nbins)),
        m_xMin(xmin),
        m_xMax(xmax),
        m_parameters(m_maxParams),
        m_parErrors(m_maxParams),
        m_parSettings(npar)
   {
   }

   IModel(const IModel &other)
      : m_nBins(other.m_nBins),
        m_wBin(other.m_wBin),
        m_step(other.m_step),
        m_xMin(other.m_xMin),
        m_xMax(other.m_xMax),
        m_parameters(other.m_parameters),
        m_parErrors(other.m_parErrors),
        m_parSettings(other.m_parSettings)
   {
   }

private:
   virtual Double_t DoEvalPar(Double_t x, const Double_t *p) const override = 0;
};

#endif