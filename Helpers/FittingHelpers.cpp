#include "FittingHelpers.h"

using namespace constants;

namespace calib {

    double_t langaufun(Double_t *x, Double_t *par) {

        //Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
        //Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
        //Markus Friedl (Markus.Friedl@cern.ch) "langaus.C"

        //Fit parameters:
        //par[0]=Width (scale) parameter of Landau density
        //par[1]=Most Probable (MP, location) parameter of Landau density
        //par[2]=Total area (integral -inf to inf, normalization constant)
        //par[3]=Width (sigma) of convoluted Gaussian function
        //
        //In the Landau distribution (represented by the CERNLIB approximation), 
        //the maximum is located at x=-0.22278298 with the location parameter=0.
        //This shift is corrected within this function, so that the actual
        //maximum is identical to the MP parameter.

        // Numeric constants
        Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
        Double_t mpshift  = -0.22278298;       // Landau maximum location

        // Control constants
        Double_t np = 100.0;      // number of convolution steps
        Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

        // Variables
        Double_t xx;
        Double_t mpc;
        Double_t fland;
        Double_t sum = 0.0;
        Double_t xlow,xupp;
        Double_t step;
        Double_t i;


        // MP shift correction
        mpc = par[1] - mpshift * par[0]; 

        // Range of convolution integral
        xlow = x[0] - sc * par[3];
        xupp = x[0] + sc * par[3];

        step = (xupp-xlow) / np;

        // Convolution integral of Landau and Gaussian by sum
        for(i=1.0; i<=np/2; i++) {
            xx = xlow + (i-.5) * step;
            fland = TMath::Landau(xx,mpc,par[0]) / par[0];
            sum += fland * TMath::Gaus(x[0],xx,par[3]);

            xx = xupp - (i-.5) * step;
            fland = TMath::Landau(xx,mpc,par[0]) / par[0];
            sum += fland * TMath::Gaus(x[0],xx,par[3]);
        }

        return (par[2] * step * sum * invsq2pi / par[3]);
    }

    double_t expofunX(Double_t *x, Double_t *par){

        Double_t xx = x[0];
        Double_t f = par[0]*exp(-((200-TMath::Abs(xx))/((1/par[1])*vDrift)));
        return f;
    }

    double_t expofunT(Double_t *x, Double_t *par){

        Double_t xx = x[0];
        Double_t f = par[0]*exp(-(xx*par[1]));
        return f;
    }

    double_t expofunConst(Double_t *x, Double_t *par){

        Double_t xx = x[0];
        Double_t f = par[0]*exp(-(xx*par[1])) + par[2];
        return f;
    }

    void SetLGParameters(TH1D *h, Double_t *fp, Double_t *ehfp, Double_t *elfp, Double_t &lb, Double_t &ub){
        //Adapted from Lan Nguyen's code which was adapted from Dom Barker's code etc.

        // Set fit parameters and errors  
        double nentries = h->GetEntries();
        double max = h->GetBinCenter(h->GetMaximumBin()); //x position of max dQdx value 
        double area = h->GetEntries()*h->GetBinWidth(1)*0.8 ; //normalisation constant //starting value for the fit - fit doesn't include whole hist
        double rms = h->GetRMS();
        
        //Starting guesses for parameters
        fp[0] = max / 100 * 4; //scale parameter of Landau density, starting guess about 5% of MPV x value
        fp[1] = max; //MPV: max dQdx value 		(position)
        fp[2] = area; // total area or integral //NORM (bin width) getNentries*Getbinwidth
        fp[3] = rms / 10;  //GSigma typically smaller than distribution rms
        
        //Starting guesses for errors
        ehfp[0] = max / 100 * 0.5; //0.5% of max dQdx value  	
        ehfp[1] = max / 100 * 8;	// 8% of max dQdx value
        ehfp[2] = area * 0.1;	
        ehfp[3] = rms * 0.01;
        elfp[0] = - max / 100 * 0.5; //0.5% of max dQdx value  	
        elfp[1] = - max / 100 * 8;	// 8% of max dQdx value
        elfp[2] = - area * 0.1;	
        elfp[3] = - rms * 0.01;
    
        //Lower and upper bounds for fit: sufficient to find the peak of the distribution precisely,
        //whilst staying in high stats bins to reduce statistical fluctuations
        double dQdxpeak = h->GetMaximumBin(); 
        lb = h->GetBinCenter(dQdxpeak-8);	
        ub = h->GetBinCenter(dQdxpeak+12);

    }

    void SetExpoParameters(Double_t *fp){
        fp[0] = 1000.;
        fp[1] = 0.1;

    }

    void SetExpoConstParameters(Double_t *fp){
        fp[0] = 3.44;
        fp[1] = 1.0;
        fp[2] = 11.;

    }

    fitResult fitter(TH1D *h, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName){

        fitResult postFit;

        Double_t(*func)(Double_t *,Double_t *);
        int func_index;
        int nParams;

        if(!strcmp(funcName.c_str(), "expoX")){
            func = expofunX;
            func_index = 1;
            nParams = 2;
        }
        else if(!strcmp(funcName.c_str(), "expoT")){
            func = expofunT;
            func_index = 2;
            nParams = 2;
        }
        else if(!strcmp(funcName.c_str(), "LG")){
            func = langaufun;
            func_index = 3;
            nParams = 4;
        }
        else{
            std::cout << "UNKNOWN FUNCTION" << std::endl;
        }

        Char_t FitFuncName[100]; 
        sprintf(FitFuncName,"Fitfcn_%s",h->GetName());

        TF1 *fitfunc = new TF1(FitFuncName,func,lbound,ubound, nParams);
        
        if(func_index == 1 || func_index == 2){
            fitfunc->SetParameters(fitparams[0], fitparams[1]);
            fitfunc->SetParError(0,ehfp[0]);
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParNames("Norm","1/Lifetime");
        }
        else if(func_index == 3){
            fitfunc->SetParameters(fitparams[0], fitparams[1], fitparams[2], fitparams[3]);
            fitfunc->SetParError(0,ehfp[0]);
            if (ehfp[0]==0) fitfunc->FixParameter(0,fitparams[0]); //if scale parameter error is 0 scale parameter is fixed
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParError(2,ehfp[2]);
            fitfunc->SetParError(3,ehfp[3]);
            fitfunc->SetParNames("Width","MPV","TotalArea","GSigma"); 
        }

        TFitResultPtr r = h->Fit(FitFuncName,"LREQS");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode
        //Other fitting options https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html (7.1.1)
    
        TString fitOutcome = gMinuit->fCstatu.Data();

        for(int i = 0; i < nParams; i++){
            fitparams[i] = h->GetFunction(FitFuncName)->GetParameter(i);	
            ehfp[i] = r->UpperError(i);
            elfp[i] = r->LowerError(i);
        }

        if (!fitOutcome.BeginsWith("SUCC")) { 
            for(int i = 0; i < nParams; i++){
                fitparams[i] = -1000.0;	
                ehfp[i] = -1000.0;
                elfp[i] = 1000.0;
                covmat[i] = -1000.0;
            } 
        }

        postFit.myFitFunc = fitfunc;
        postFit.myFitResult = r;

        return(postFit);
        
    }

    fitResult fitter(TGraphAsymmErrors *g, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName){

        fitResult postFit;

        Double_t(*func)(Double_t *,Double_t *);
        int func_index;
        int nParams;

        if(!strcmp(funcName.c_str(), "expoX")){
            func = expofunX;
            func_index = 1;
            nParams = 2;
        }
        else if(!strcmp(funcName.c_str(), "expoT")){
            func = expofunT;
            func_index = 2;
            nParams = 2;
        }
        else if(!strcmp(funcName.c_str(), "LG")){
            func = langaufun;
            func_index = 3;
            nParams = 4;
        }
        else if(!strcmp(funcName.c_str(), "expoConst")){
            func = expofunConst;
            func_index = 4;
            nParams = 3;
        }
        else{
            std::cout << "UNKNOWN FUNCTION" << std::endl;
        }

        Char_t FitFuncName[100]; 
        sprintf(FitFuncName,"Fitfcn_%s",g->GetName());

        TF1 *fitfunc = new TF1(FitFuncName,func,lbound,ubound, nParams);
        
        if(func_index == 1 || func_index == 2){
            fitfunc->SetParameters(fitparams[0], fitparams[1]);
            fitfunc->SetParError(0,ehfp[0]);
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParNames("Norm","1/Lifetime");
        }
        else if(func_index == 3){
            fitfunc->SetParameters(fitparams[0], fitparams[1], fitparams[2], fitparams[3]);
            fitfunc->SetParError(0,ehfp[0]);
            if (ehfp[0]==0) fitfunc->FixParameter(0,fitparams[0]); //if scale parameter error is 0 scale parameter is fixed
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParError(2,ehfp[2]);
            fitfunc->SetParError(3,ehfp[3]);
            fitfunc->SetParNames("Width","MPV","TotalArea","GSigma"); 
        }
        else if(func_index == 4){
            fitfunc->SetParameters(fitparams[0], fitparams[1], fitparams[2]);
            fitfunc->SetParError(0,ehfp[0]);
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParError(2,ehfp[2]);
            fitfunc->SetParNames("Norm","Decay","Const");
        }

        TFitResultPtr r = g->Fit(FitFuncName,"LREQS");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode
        //Other fitting options https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html (7.1.1)
    
        TString fitOutcome = gMinuit->fCstatu.Data();

        for(int i = 0; i < nParams; i++){
            fitparams[i] = g->GetFunction(FitFuncName)->GetParameter(i);	
            ehfp[i] = r->UpperError(i);
            elfp[i] = r->LowerError(i);
        }

        if (!fitOutcome.BeginsWith("SUCC")) { 
            for(int i = 0; i < nParams; i++){
                fitparams[i] = -1000.0;	
                ehfp[i] = -1000.0;
                elfp[i] = 1000.0;
                covmat[i] = -1000.0;
            } 
        }

        postFit.myFitFunc = fitfunc;
        postFit.myFitResult = r;

        return(postFit);
        
    }

    double** getAngBinLimits(int angBins){

        double angWindow = 136.4;
		double angLowLimit1 = -68.2;
		double angUpLimit1 = 68.2;
		double angLowLimit2 = 111.8;
		double angUpLimit2 = 248.2;

		double angStep = angWindow*2/(double)angBins;

        double** angLimits = new double*[2];
        angLimits[0] = new double[angBins/2+1];
        angLimits[1] = new double[angBins/2+1];
		
        angLimits[0][0] = -68.2;
        angLimits[1][0] = 111.8;

		for(int i = 1; i <= angBins/2; i++){

			angLimits[0][i] = angLimits[0][i-1] + angStep;
			angLimits[1][i] = angLimits[1][i-1] + angStep;

		}

        return angLimits;

    }

}
