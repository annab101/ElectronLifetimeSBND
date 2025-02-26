// Author: Anna Beever
// Date:   January 2025
// Code for monitoring the electron lifetime in the nearline. The lifetime values are for 
// purity monitoring, not for use in calibration. Produces 1D histograms of dQ/dx in each time
// bin for both TPCs, with the LG fit also shown; a plot of the MPV distribution and exponential
// fit for both TPCs side by side; and finally a text file with the fit results for 1/lifetime,
// and these converted to lifetime.

//C++ includes
#include<fstream>
#include<filesystem>
#include<valarray>

//ROOT includes

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TError.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

namespace constants {

    //Sheffield colour scheme
	Int_t deepViolet = TColor::GetColor("#440099");
	Int_t coral = TColor::GetColor("#E7004C");
	Int_t powderBlue = TColor::GetColor("#9ADBE8");
	Int_t flamingo = TColor::GetColor("#FF6371");
	Int_t mintGreen = TColor::GetColor("#00CE7C");
	Int_t purple = TColor::GetColor("#981F92");
	Int_t peach= TColor::GetColor("#FF9664");
	Int_t lavender= TColor::GetColor("#DAA8E2");
	Int_t mauve= TColor::GetColor("#663DB3");
	double vDrift = 156.267;

	std::map<int, std::string> configMap{{0, "classic"}, {1, "multiwire"}, {2, "multihit"}, {3, "angle"}, {4, "wireAndAngle"}, {5, "octant"}};
	std::map<int, std::string> octantName{{1, "NW Upper"}, {2, "NE Upper"}, {3, "SW Upper"}, {4, "SE Upper"}, {5, "NW Lower"}, {6, "NE Lower"}, {7, "SW Lower"}, {8, "SE Lower"}};
	std::map<int, std::string> plotFitName{{0, "noFit"}, {1, "withFit"}};
	std::map<std::string, std::pair<std::string, std::pair<double, double>>> whichSetup{{"XL", std::make_pair("xDriftL", std::make_pair(-200.,0.))}, {"XR", std::make_pair("xDriftR", std::make_pair(0.,200.))}, {"TL", std::make_pair("tDriftL", std::make_pair(0., 1.3))}, {"TR", std::make_pair("tDriftR", std::make_pair(0., 1.3))}};

}

namespace calib {

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

    double_t expofunT(Double_t *x, Double_t *par){

        Double_t xx = x[0];
        Double_t f = par[0]*exp(-(xx*par[1]));
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
        ub = h->GetBinCenter(dQdxpeak+16);

    }

    void SetExpoParameters(Double_t *fp){
        fp[0] = 1000.;
        fp[1] = 0.1;

    }

    TF1 *fitter(TH1D *h, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName){

        Double_t(*func)(Double_t *,Double_t *);
        int func_index;
        int nParams;

        if(!strcmp(funcName.c_str(), "LG")){
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
        
        if(func_index == 3){
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

        return(fitfunc);
        
    }

    TF1 *fitter(TGraphAsymmErrors *g, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *ehfp, Double_t *elfp, Double_t *covmat, std::string funcName){

        Double_t(*func)(Double_t *,Double_t *);
        int func_index;
        int nParams;

        if(!strcmp(funcName.c_str(), "expoT")){
            func = expofunT;
            func_index = 2;
            nParams = 2;
        }
        else{
            std::cout << "UNKNOWN FUNCTION" << std::endl;
        }

        Char_t FitFuncName[100]; 
        sprintf(FitFuncName,"Fitfcn_%s",g->GetName());

        TF1 *fitfunc = new TF1(FitFuncName,func,lbound,ubound, nParams);
        
        if(func_index == 2){
            fitfunc->SetParameters(fitparams[0], fitparams[1]);
            fitfunc->SetParError(0,ehfp[0]);
            fitfunc->SetParError(1,ehfp[1]);
            fitfunc->SetParNames("Norm","DecayConst");
        }

        TFitResultPtr r = g->Fit(FitFuncName,"LREQS0");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode, 0 = don't plot fit
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

        return(fitfunc);
        
    }

    template<class T> void setMarker(T *t, Color_t color, Style_t style, Size_t size){
        
        t->SetMarkerColor(color);
        t->SetMarkerStyle(style);
        t->SetMarkerSize(size);
    }

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, fitExpoParameters expoParams){
        
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f1->GetChisquare(),f1->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("1/Lifetime: %.5f^{+%.5f}_{%.5f}",f1->GetParameter(1),expoParams.ehfp[1],expoParams.elfp[1]));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        
        return pt;
    }

    template<class T> void setFontSize(T *t, int font, int size){
	
        t->GetXaxis()->SetTitleFont(font);
        t->GetXaxis()->SetTitleSize(size);
        t->GetXaxis()->SetLabelFont(font);
        t->GetXaxis()->SetLabelSize(size);

        t->GetYaxis()->SetTitleFont(font);
        t->GetYaxis()->SetTitleSize(size);
        t->GetYaxis()->SetLabelFont(font);
        t->GetYaxis()->SetLabelSize(size);

    }

}

using namespace constants;
using namespace calib;

int main(int argc, char**argv) {

    // Input file should be a list of all the calib ntuple .root files.

    //=================================================================
    // READING IN DATA FILES & SETTING PARAMETERS
    //=================================================================

	std::string filename = "noFile";

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--filename")){
			filename = argv[i+1];
		}
	}

	std::cout << "-----------------------------------------------------------" << std::endl;

	std::string tag = "nearline";

    std::string textFile = "LifetimeValues_" + tag + ".txt";
	
	int nGroupedWires = 1;
	int NBinsT = 10;
	int NBinsdQdx = 50;
	double mindQdx = 600.;
	double maxdQdx = 2400.;
	double minT = 0.;
	double maxT = 1.3;
    int writeProjY = 1;
    double expo_tDriftE_lb = 0.2;
    double expo_tDriftW_lb = 0.2;
    double expo_tDriftE_ub = 1.0;
    double expo_tDriftW_ub = 1.0;

	TChain *chain = new TChain;

	std::fstream file_list;
	file_list.open(filename.c_str(),std::ios::in);
	
	if (file_list.is_open()){
		std::string fileUrl;
		while(getline(file_list, fileUrl)){

			chain->Add((fileUrl + "/caloskim/TrackCaloSkim").c_str());
			
		}
		file_list.close();
	}

	//read in variables of interest
	TTreeReader treereader(chain);
	TTreeReaderArray<Float_t> read_x(treereader, "trk.hits2.h.sp.x");
	TTreeReaderArray<Float_t> read_y(treereader, "trk.hits2.h.sp.y");
	TTreeReaderArray<Float_t> read_z(treereader, "trk.hits2.h.sp.z");
	TTreeReaderArray<Float_t> read_dqdx(treereader, "trk.hits2.dqdx");
	TTreeReaderArray<Float_t> read_T(treereader, "trk.hits2.h.time");
	TTreeReaderValue<Float_t> read_t0(treereader, "trk.t0");
	TTreeReaderValue<Float_t> read_xi(treereader, "trk.start.x");
	TTreeReaderValue<Float_t> read_yi(treereader, "trk.start.y");
	TTreeReaderValue<Float_t> read_zi(treereader, "trk.start.z");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<Float_t> read_yf(treereader, "trk.end.y");
	TTreeReaderValue<Float_t> read_zf(treereader, "trk.end.z");
	TTreeReaderValue<int> read_event(treereader, "trk.meta.evt");
	TTreeReaderValue<int> read_selected(treereader, "trk.selected");
	TTreeReaderArray<uint16_t> read_wire(treereader, "trk.hits2.h.wire");

    //=================================================================
    // 2D HISTOGRAM OF dQ/dx vs t
    //=================================================================

	int track_count = 0;

	TH2D* h_dQdx_tDriftE = new TH2D(TString::Format("h_dQdx_tDriftE_%dwires", nGroupedWires), "dQ/dx vs t East TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	TH2D* h_dQdx_tDriftW = new TH2D(TString::Format("h_dQdx_tDriftW_%dwires", nGroupedWires), "dQ/dx vs t West TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

    std::cout << "Reading in trees" << std::endl;

	treereader.Restart();
    while(treereader.Next()){

        if(*read_selected != 1){
            continue;
        }

        if(track_count % 10 == 0){std::cout << "track_count: " << track_count << std::endl;}
        track_count ++;

        auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
        int start_index = std::distance(read_x.begin(), it_start_x);

        if(start_index == read_x.GetSize()){continue;}
        auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
        int end_index = std::distance(read_x.begin(), it_end_x);
        if(end_index < 0){continue;}

        
        for(int i = start_index; i <= end_index; i++){            
                                    
            if(-200. < read_x[i] && read_x[i] < 0.){
                h_dQdx_tDriftE->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
            }

            if(0. < read_x[i] && read_x[i] < 200.){
                h_dQdx_tDriftW->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
            }

        }

    }

    //=================================================================
    // 1D PROJECTION IN EACH TIME BIN + LG FIT + PLOT FOR EACH BIN
    //=================================================================
	
	TH1D* projY_tDriftE = new TH1D();
	TH1D* projY_tDriftW = new TH1D();

	TF1* LGfit_tDriftE = new TF1();
	TF1* LGfit_tDriftW = new TF1();

	fitLGParameters fp_tDriftE;
	fitLGParameters fp_tDriftW;

	TGraphAsymmErrors *MPVplot_tDriftE = new TGraphAsymmErrors();
	TGraphAsymmErrors *MPVplot_tDriftW = new TGraphAsymmErrors();

	for(int i = 1; i <= NBinsT; i++){

		std::cout << "--------- Fitting bin " + std::to_string(i) + "----------" << std::endl;

		projY_tDriftE = h_dQdx_tDriftE->ProjectionY(TString::Format("projY_tDriftE_bin%d", i), i, i);
		projY_tDriftW = h_dQdx_tDriftW->ProjectionY(TString::Format("projY_tDriftW_bin%d", i), i, i);

		SetLGParameters(projY_tDriftE, fp_tDriftE.fp, fp_tDriftE.ehfp, fp_tDriftE.elfp, fp_tDriftE.lb, fp_tDriftE.ub);
		SetLGParameters(projY_tDriftW, fp_tDriftW.fp, fp_tDriftW.ehfp, fp_tDriftW.elfp, fp_tDriftW.lb, fp_tDriftW.ub);

		LGfit_tDriftE = fitter(projY_tDriftE, fp_tDriftE.lb, fp_tDriftE.ub, fp_tDriftE.fp, fp_tDriftE.ehfp, fp_tDriftE.elfp, fp_tDriftE.cov, "LG");
		LGfit_tDriftW = fitter(projY_tDriftW, fp_tDriftW.lb, fp_tDriftW.ub, fp_tDriftW.fp, fp_tDriftW.ehfp, fp_tDriftW.elfp, fp_tDriftW.cov, "LG");

		if(fp_tDriftE.fp[1] >= 0.){

			MPVplot_tDriftE->SetPoint(MPVplot_tDriftE->GetN(), h_dQdx_tDriftE->GetXaxis()->GetBinCenter(i), fp_tDriftE.fp[1]);
			MPVplot_tDriftE->SetPointError(MPVplot_tDriftE->GetN() - 1, 0., 0., -1*fp_tDriftE.elfp[1], fp_tDriftE.ehfp[1]);
		}

		if(fp_tDriftW.fp[1] >= 0.){

			MPVplot_tDriftW->SetPoint(MPVplot_tDriftW->GetN(), h_dQdx_tDriftW->GetXaxis()->GetBinCenter(i), fp_tDriftW.fp[1]);
			MPVplot_tDriftW->SetPointError(MPVplot_tDriftW->GetN() - 1, 0., 0., -1*fp_tDriftW.elfp[1], fp_tDriftW.ehfp[1]);
		}

		if((bool)writeProjY){
	
            std::stringstream s1_tDriftE;
            s1_tDriftE << std::setprecision(4) << h_dQdx_tDriftE->GetXaxis()->GetBinLowEdge(i);
            std::string str1_tDriftE = s1_tDriftE.str();
            std::stringstream s2_tDriftE;
            s2_tDriftE << std::setprecision(4) << (h_dQdx_tDriftE->GetXaxis()->GetBinLowEdge(i) + h_dQdx_tDriftE->GetXaxis()->GetBinWidth(i));
            std::string str2_tDriftE = s2_tDriftE.str();
            std::string projYtitle_tDriftE = "dQ/dx for " + str1_tDriftE + " ms < t < " + str2_tDriftE + " ms, East TPC;dQ/dx (Arb. Units);Entries";
            projY_tDriftE->SetTitle(projYtitle_tDriftE.c_str());
            projY_tDriftE->SetTitleFont(133);
            projY_tDriftE->SetTitleSize(28);

            projY_tDriftE->GetXaxis()->SetLabelOffset(0.01);
            projY_tDriftE->GetXaxis()->SetNdivisions(505);
            setFontSize<TH1>(projY_tDriftE, 133, 30);
            projY_tDriftE->GetYaxis()->SetLabelOffset(0.01);
            projY_tDriftE->SetFillColor(flamingo);
            projY_tDriftE->SetLineColor(flamingo);
            projY_tDriftE->SetStats(0);
            LGfit_tDriftE->SetLineColor(deepViolet);
            LGfit_tDriftE->SetLineWidth(2.3);
        
            TCanvas *c_east = new TCanvas();
            c_east->SetLeftMargin(0.15);
            c_east->SetBottomMargin(0.14);
            projY_tDriftE->Draw();
            LGfit_tDriftE->Draw("same");
            c_east->SaveAs(("nearlinedQdxProjection_tDriftE_bin" + std::to_string(i) + ".png").c_str());

            std::stringstream s1_tDriftW;
            s1_tDriftW << std::setprecision(4) << h_dQdx_tDriftW->GetXaxis()->GetBinLowEdge(i);
            std::string str1_tDriftW = s1_tDriftW.str();
            std::stringstream s2_tDriftW;
            s2_tDriftW << std::setprecision(4) << (h_dQdx_tDriftW->GetXaxis()->GetBinLowEdge(i) + h_dQdx_tDriftW->GetXaxis()->GetBinWidth(i));
            std::string str2_tDriftW = s2_tDriftW.str();
            std::string projYtitle_tDriftW = "dQ/dx for " + str1_tDriftW + " ms < t < " + str2_tDriftW + " ms, West TPC;dQ/dx (Arb. Units);Entries";
            projY_tDriftW->SetTitle(projYtitle_tDriftW.c_str());
            projY_tDriftW->SetTitleFont(133);
            projY_tDriftW->SetTitleSize(28);

            projY_tDriftW->GetXaxis()->SetLabelOffset(0.01);
            projY_tDriftW->GetXaxis()->SetNdivisions(505);
            setFontSize<TH1>(projY_tDriftW, 133, 30);
            projY_tDriftW->GetYaxis()->SetLabelOffset(0.01);
            projY_tDriftW->SetFillColor(flamingo);
            projY_tDriftW->SetLineColor(flamingo);
            projY_tDriftW->SetStats(0);
            LGfit_tDriftW->SetLineColor(deepViolet);
            LGfit_tDriftW->SetLineWidth(2.3);
        
            TCanvas *c_west = new TCanvas();
            c_west->SetLeftMargin(0.15);
            c_west->SetBottomMargin(0.14);
            projY_tDriftW->Draw();
            LGfit_tDriftW->Draw("same");
            c_west->SaveAs(("nearlinedQdxProjection_tDriftW_bin" + std::to_string(i) + ".png").c_str());
		}

	}

    //=================================================================
    // EXPO FIT TO MPVS FROM EACH LG FIT + MPV PLOT + LIFETIME VALUES
    //=================================================================

	std::cout << "-----Initialising expo fit parameters-----" << std::endl;

	fitExpoParameters expoParams_tDriftE;
	fitExpoParameters expoParams_tDriftW;

	TF1 *expoFit_tDriftE = new TF1();
	TF1 *expoFit_tDriftW = new TF1();

	expoParams_tDriftE.reset();
	expoParams_tDriftW.reset();

	SetExpoParameters(expoParams_tDriftE.fp);
	SetExpoParameters(expoParams_tDriftW.fp);

	expoParams_tDriftE.lb = expo_tDriftE_lb;
	expoParams_tDriftE.ub = expo_tDriftE_ub;
	expoParams_tDriftW.lb = expo_tDriftW_lb;
	expoParams_tDriftW.ub = expo_tDriftW_ub;

	std::cout << "-----Fitting exponential-----" << std::endl;

	expoFit_tDriftE = fitter(MPVplot_tDriftE, expoParams_tDriftE.lb, expoParams_tDriftE.ub, expoParams_tDriftE.fp, expoParams_tDriftE.ehfp, expoParams_tDriftE.elfp, expoParams_tDriftE.cov, "expoT");
	expoFit_tDriftW = fitter(MPVplot_tDriftW, expoParams_tDriftW.lb, expoParams_tDriftW.ub, expoParams_tDriftW.fp, expoParams_tDriftW.ehfp, expoParams_tDriftW.elfp, expoParams_tDriftW.cov, "expoT");

	expoFit_tDriftE->SetName("MPVfit_tDriftE");
	expoFit_tDriftW->SetName("MPVfit_tDriftW");

	std::cout << "-----Saving plots and writing lifetime values to file-----" << std::endl;

    //plotting MPVs and expo fit goes here

	int MPVmin = 600;
	int MPVmax = 2400;

    setMarker<TGraphAsymmErrors>(MPVplot_tDriftE, 1, 20, 0.7);
    MPVplot_tDriftE->SetTitle("East TPC;t_{#scale[1.2]{drift}} (ms); dQ/dx MPV (Arb. Units)");
    MPVplot_tDriftE->SetMinimum(MPVmin);
    MPVplot_tDriftE->SetMaximum(MPVmax);
    expoFit_tDriftE->SetLineColor(coral);
    expoFit_tDriftE->SetLineWidth(2.0);
    setFontSize<TGraphAsymmErrors>(MPVplot_tDriftE, 133, 30);

    setMarker<TGraphAsymmErrors>(MPVplot_tDriftW, 1, 20, 0.7);
    MPVplot_tDriftW->SetTitle("West TPC;t_{#scale[1.2]{drift}} (ms); dQ/dx MPV (Arb. Units)");
    MPVplot_tDriftW->SetMinimum(MPVmin);
    MPVplot_tDriftW->SetMaximum(MPVmax);
    expoFit_tDriftW->SetLineColor(kAzure - 3);
    expoFit_tDriftW->SetLineWidth(2.0);
    setFontSize<TGraphAsymmErrors>(MPVplot_tDriftW, 133, 30);

    TCanvas *c_split = new TCanvas();
    c_split->SetWindowSize(2000,500);
    
    TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
    TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
    c_split->cd();
    pad1->SetRightMargin(0.15);
    pad1->SetLeftMargin(0.18);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.15);
    pad2->SetRightMargin(0.15);
    pad2->SetLeftMargin(0.18);
    pad2->SetTopMargin(0.08);
    pad2->SetBottomMargin(0.15);    
    c_split->cd();		
    pad1->Draw();		
    pad1->cd();	
    MPVplot_tDriftE->Draw("AP");	
    expoFit_tDriftE->Draw("SAME");
    TPaveText *statsE = statsBox({.45,.65,.80,.88}, track_count, expoFit_tDriftE, expoParams_tDriftE);
    statsE->Draw();
    c_split->cd();
    pad2->Draw();
    pad2->cd();
    MPVplot_tDriftW->Draw("AP");
    expoFit_tDriftW->Draw("SAME");
    TPaveText *statsW = statsBox({.45,.65,.80,.88}, track_count, expoFit_tDriftW, expoParams_tDriftW);
    statsW->Draw();

    c_split->SaveAs("nearlineMPVs.png");

	std::ofstream outFile;
	outFile.open(textFile.c_str());
	outFile << std::setw(20) << "TPC  |  "
			<< std::setw(20) << "1/etime  |  "
			<< std::setw(20) << "error low  |  "
			<< std::setw(20) << "error up  | \n"
			<< std::setw(20) << "E  |  "
			<< std::setw(15) << expoParams_tDriftE.fp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftE.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftE.ehfp[1] << "  | \n"
			<< std::setw(20) << "W  |  "
			<< std::setw(15) << expoParams_tDriftW.fp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftW.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftW.ehfp[1] << "  | \n"
            << std::setw(20) << "TPC  |  "
			<< std::setw(20) << "etime  |  "
			<< std::setw(20) << "error low  |  "
			<< std::setw(20) << "error up  | \n"
			<< std::setw(20) << "E  |  "
			<< std::setw(15) << 1/expoParams_tDriftE.fp[1] << "  |  "
			<< std::setw(15) << -((1/expoParams_tDriftE.fp[1]) - (1/(expoParams_tDriftE.fp[1] + expoParams_tDriftE.ehfp[1]))) << "  |  "
			<< std::setw(15) << ((1/(expoParams_tDriftE.fp[1]+expoParams_tDriftE.elfp[1])) - (1/expoParams_tDriftE.fp[1])) << "  | \n"
			<< std::setw(20) << "W  |  "
			<< std::setw(15) << 1/expoParams_tDriftW.fp[1] << "  |  "
			<< std::setw(15) << -((1/expoParams_tDriftW.fp[1]) - (1/(expoParams_tDriftW.fp[1] + expoParams_tDriftW.ehfp[1]))) << "  |  "
			<< std::setw(15) << ((1/(expoParams_tDriftW.fp[1]+expoParams_tDriftW.elfp[1])) - (1/expoParams_tDriftW.fp[1])) << "  | \n"
            << "Note: 1/etime of -1000 means the fit did not work!\n";
	outFile.close();

	return 0;

}