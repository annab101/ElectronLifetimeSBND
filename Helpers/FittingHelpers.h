/*
* Helper functions for fitting
*/

#ifndef FITTINGHELPERS_H
#define FITTINGHELPERS_H

namespace calib {

    /**
    * Components for LG fit
    *
    * Contains parameters, errors, covariance matrix, upper and lower bound,
    * plus a function to reset parameters to default value
    */

    class fitLGParameters {
	
        public: 
            

            double fp[4] = {0.}; //fit parameters
            double ehfp[4] = {0.}; //fit parameter upper errors
            double elfp[4] = {0.}; //fit parameter lower errors
            double cov[4] = {0.}; //covariance matrix
            double lb = -1.; //fit lower bound
            double ub = -1.;

            void reset() {
                std::fill(std::begin(fp), std::end(fp), 0.);
                std::fill(std::begin(ehfp), std::end(ehfp), 0.);
                std::fill(std::begin(elfp), std::end(elfp), 0.);
                std::fill(std::begin(cov), std::end(cov), 0.);
                lb = -1.;
                ub = -1.;
            }

    };

    /**
    * Components for expo fit
    *
    * Contains parameters, errors, covariance matrix, upper and lower bound,
    * plus a function to reset parameters to default value
    */

    class fitExpoParameters {
        
        public: 
            

            double fp[2] = {0.}; //fit parameters
            double ehfp[2] = {0.}; //fit parameter errors
            double elfp[2] = {0.}; //fit parameter errors
            double cov[2] = {0.}; //covariance matrix
            double lb = -1.; //fit lower bound
            double ub = -1.;

            void reset() {
                std::fill(std::begin(fp), std::end(fp), 0.);
                std::fill(std::begin(ehfp), std::end(ehfp), 0.);
                std::fill(std::begin(elfp), std::end(elfp), 0.);
                std::fill(std::begin(cov), std::end(cov), 0.);
                lb = -1.;
                ub = -1.;
            }
            
    };

    /**
    * Components for expo fit with constant
    *
    * Contains parameters, errors, covariance matrix, upper and lower bound,
    * plus a function to reset parameters to default value
    */

    class fitExpoConstParameters {
        
        public: 
            

            double fp[3] = {0.}; //fit parameters
            double ehfp[3] = {0.}; //fit parameter errors
            double elfp[3] = {0.}; //fit parameter errors
            double cov[3] = {0.}; //covariance matrix
            double lb = -1.; //fit lower bound
            double ub = -1.;

            void reset() {
                std::fill(std::begin(fp), std::end(fp), 0.);
                std::fill(std::begin(ehfp), std::end(ehfp), 0.);
                std::fill(std::begin(elfp), std::end(elfp), 0.);
                std::fill(std::begin(cov), std::end(cov), 0.);
                lb = -1.;
                ub = -1.;
            }
            
    };

    /**
    * @brief Landau-Gaussian convolution.
    *
    * @param x x
    * @param par parameter array
    *
    */

    double_t langaufun(Double_t *x, Double_t *par);

    /**
    * @brief Exponential function in x for SBND geometry.
    *
    * @param x x
    * @param par parameter array
    *
    */

    double_t expofunX(Double_t *x, Double_t *par);

    /**
    * @brief Exponential function with normalisation and decay constant.
    *
    * @param x x
    * @param par parameter array
    *
    */

    double_t expofunT(Double_t *x, Double_t *par);

    /**
    * @brief Exponential function with normalisation, decay constant and added constant.
    *
    * @param x x
    * @param par parameter array
    *
    */

    double_t expofunConst(Double_t *x, Double_t *par);

    /**
    * @brief Initialise parameters for LG fit.
    *
    * @param h Histogram for fit
    * @param fp Fit parameters array
    * @param ehfp Upper error in fit parameters array
    * @param elfp Lower error in fit parameters array
    * @param lb Fit lower bound
    * @param ub Fit upper bound
    *
    */

    void SetLGParameters(TH1D *h, Double_t *fp, Double_t *ehfp, Double_t *elfp, Double_t &lb, Double_t &ub);
    
    /**
    * @brief Initialise parameters for exponential fit.
    *
    * @param fp Fit parameters array
    *
    */
    
    void SetExpoParameters(Double_t *fp);

    /**
    * @brief Initialise parameters for exponential fit.
    *
    * @param fp Fit parameters array
    *
    */

    void SetExpoConstParameters(Double_t *fp);

    /**
    * @brief Fit function to histogram and return fit.
    *
    * @param h Histogram for fit
    * @param lbound Fit lower bound
    * @param ubound Fit upper bound
    * @param fitparams Fit parameters array
    * @param fiterrors Array of errors in fit parameters
    * @param covmat Covariance matrix
    * @param funcName Function name 
    *
    */

    TF1 *fitter(TH1D *h, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName);
    
    /**
    * @brief Fit function to TGraphErrors and return fit.
    *
    * @param g TGraphErrors for fit
    * @param lbound Fit lower bound
    * @param ubound Fit upper bound
    * @param fitparams Fit parameters array
    * @param ehfp Array of upper errors in fit parameters
    * @param elfp Array of lower errors in fit parameters
    * @param covmat Covariance matrix
    * @param funcName Function name 
    *
    */
    
    TF1 *fitter(TGraphErrors *g, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName);

    /**
    * @brief Function to find angular bin limits.
    *
    * @param angBins number of angular bins
    *
    */

    double** getAngBinLimits(int angBins);

}

#endif
