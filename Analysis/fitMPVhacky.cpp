//Author: Anna Beever
//Date:   April 2024

//C++ includes
#include<fstream>
#include<filesystem>
#include<valarray>

//ROOT includes
#include "../Utilities/ROOTincludes.h"

//Local includes
#include "../Helpers/Constants.h"
#include "../Helpers/PlottingHelpers.h"
#include "../Helpers/PlottingHelpers.cpp"
#include "../Helpers/FittingHelpers.h"
#include "../Helpers/FittingHelpers.cpp"
#include "../Utilities/ConfigReader.h"
#include "../Utilities/ConfigReader.cpp"

using namespace calib;
using namespace constants;
using namespace cppsecrets;


int main(int argc, char**argv) {

    std::string filename = "noFile";

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
	}

	ConfigReader* p = ConfigReader::getInstance();
	p->parseFile(filename);
	std::cout << " Variables from configuration file: " << std::endl;
	p->dumpFileValues();
	std::cout << "-----------------------------------------------------------" << std::endl;

	std::string inputData = "noInputData";
	std::string tag = "noTag";
	std::string dataset = "noDataSet";
	std::string saveLoc = "noSaveLoc";
	int codeConfig = 0;
	int maxWireGroup = 20;
	int angBins = 6;
    double expo_xDriftL_lb = -160.;
    double expo_xDriftR_lb = 40.;
    double expo_xDriftL_ub = -40.;
    double expo_xDriftR_ub = 160.;
    double expo_tDriftL_lb = 0.2;
    double expo_tDriftR_lb = 0.2;
    double expo_tDriftL_ub = 1.0;
    double expo_tDriftR_ub = 1.0;
	std::string configLabel = "noConfigLabel";

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);
	p->getValue("maxWireGroup", maxWireGroup);
	p->getValue("angBins", angBins);
    p->getValue("expo_xDriftL_lb", expo_xDriftL_lb);
	p->getValue("expo_xDriftR_lb", expo_xDriftR_lb);
	p->getValue("expo_xDriftL_ub", expo_xDriftL_ub);
	p->getValue("expo_xDriftR_ub", expo_xDriftR_ub);
	p->getValue("expo_tDriftL_lb", expo_tDriftL_lb);
	p->getValue("expo_tDriftR_lb", expo_tDriftR_lb);
	p->getValue("expo_tDriftL_ub", expo_tDriftL_ub);
	p->getValue("expo_tDriftR_ub", expo_tDriftR_ub);
	
	configLabel = configMap[codeConfig];

	TFile f((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());
	TFile ff((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");

    if(codeConfig == 0){

		TGraphAsymmErrors *MPVplot_xDrift_basic = (TGraphAsymmErrors*)f.Get("MPVplot_xDrift_basic");
		TGraphAsymmErrors *MPVplot_tDriftL_basic = (TGraphAsymmErrors*)f.Get("MPVplot_tDriftL_basic");
		TGraphAsymmErrors *MPVplot_tDriftR_basic = (TGraphAsymmErrors*)f.Get("MPVplot_tDriftR_basic");

		ff.cd();

		std::cout << "-----Initialising parameters-----" << std::endl;

		fitExpoParameters expoParams_xDriftL;
		fitExpoParameters expoParams_xDriftR;
		fitExpoParameters expoParams_tDriftL;
		fitExpoParameters expoParams_tDriftR;

		TF1 *expoFit_xDriftL = new TF1();
		TF1 *expoFit_xDriftR = new TF1();
		TF1 *expoFit_tDriftL = new TF1();
		TF1 *expoFit_tDriftR = new TF1();

		expoParams_xDriftL.reset();
		expoParams_xDriftR.reset();
		expoParams_tDriftL.reset();
		expoParams_tDriftR.reset();

		SetExpoParameters(expoParams_xDriftL.fp);
		SetExpoParameters(expoParams_xDriftR.fp);
		SetExpoParameters(expoParams_tDriftL.fp);
		SetExpoParameters(expoParams_tDriftR.fp);

		expoParams_xDriftL.lb = expo_xDriftL_lb;
		expoParams_xDriftL.ub = expo_xDriftL_ub;
		expoParams_xDriftR.lb = expo_xDriftR_lb;
		expoParams_xDriftR.ub = expo_xDriftR_ub;
		expoParams_tDriftL.lb = expo_tDriftL_lb;
		expoParams_tDriftL.ub = expo_tDriftL_ub;
		expoParams_tDriftR.lb = expo_tDriftR_lb;
		expoParams_tDriftR.ub = expo_tDriftR_ub;

		std::cout << "-----Fitting exponential-----" << std::endl;

		expoFit_xDriftL = fitter(MPVplot_xDrift_basic, expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.ehfp, expoParams_xDriftL.elfp, expoParams_xDriftL.cov, "expoX");
		expoFit_xDriftR = fitter(MPVplot_xDrift_basic, expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.ehfp, expoParams_xDriftR.elfp, expoParams_xDriftR.cov, "expoX");
		std::cout << "fit tDrift L: " << std::endl;
		expoFit_tDriftL = fitter(MPVplot_tDriftL_basic, expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.ehfp, expoParams_tDriftL.elfp, expoParams_tDriftL.cov, "expoT");
		std::cout << "fit tDrift R: " << std::endl;
		expoFit_tDriftR = fitter(MPVplot_tDriftR_basic, expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.ehfp, expoParams_tDriftR.elfp, expoParams_tDriftR.cov, "expoT");

		expoFit_xDriftL->SetName("MPVfit_xDriftL_basic");
		expoFit_xDriftR->SetName("MPVfit_xDriftR_basic");
		expoFit_tDriftL->SetName("MPVfit_tDriftL_basic");
        expoFit_tDriftR->SetName("MPVfit_tDriftR_basic");

		std::cout << "-----Writing to file-----" << std::endl;
        
		expoFit_xDriftL->Write();
		expoFit_xDriftR->Write();
		expoFit_tDriftL->Write();
        expoFit_tDriftR->Write();

	}

    if(codeConfig == 1){

		fitExpoParameters expoParams_xDriftL;
		fitExpoParameters expoParams_xDriftR;
		fitExpoParameters expoParams_tDriftL;
		fitExpoParameters expoParams_tDriftR;

        TF1 *expoFit_xDriftL = new TF1();
		TF1 *expoFit_xDriftR = new TF1();
		TF1 *expoFit_tDriftL = new TF1();
		TF1 *expoFit_tDriftR = new TF1();

        ff.cd();
        
        int i = 10;

            TGraphAsymmErrors *MPVplot_xDrift_wires = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_%dwires", i));
            TGraphAsymmErrors *MPVplot_tDriftL_wires = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_%dwires", i));
            TGraphAsymmErrors *MPVplot_tDriftR_wires = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_%dwires", i));

            expoParams_xDriftL.reset();
            expoParams_xDriftR.reset();
            expoParams_tDriftL.reset();
            expoParams_tDriftR.reset();

            SetExpoParameters(expoParams_xDriftL.fp);
            SetExpoParameters(expoParams_xDriftR.fp);
            SetExpoParameters(expoParams_tDriftL.fp);
            SetExpoParameters(expoParams_tDriftR.fp);

            expoParams_xDriftL.lb = expo_xDriftL_lb;
            expoParams_xDriftL.ub = expo_xDriftL_ub;
            expoParams_xDriftR.lb = expo_xDriftR_lb;
            expoParams_xDriftR.ub = expo_xDriftR_ub;
            expoParams_tDriftL.lb = expo_tDriftL_lb;
            expoParams_tDriftL.ub = expo_tDriftL_ub;
            expoParams_tDriftR.lb = expo_tDriftR_lb;
            expoParams_tDriftR.ub = expo_tDriftR_ub;

            expoFit_xDriftL = fitter(MPVplot_xDrift_wires, expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.ehfp, expoParams_xDriftL.elfp, expoParams_xDriftL.cov, "expoX");
			expoFit_xDriftR = fitter(MPVplot_xDrift_wires, expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.ehfp, expoParams_xDriftR.elfp, expoParams_xDriftR.cov, "expoX");
			std::cout << "tL" << std::endl;
			expoFit_tDriftL = fitter(MPVplot_tDriftL_wires, expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.ehfp, expoParams_tDriftL.elfp, expoParams_tDriftL.cov, "expoT");
			std::cout << "tR" << std::endl;
			expoFit_tDriftR = fitter(MPVplot_tDriftR_wires, expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.ehfp, expoParams_tDriftR.elfp, expoParams_tDriftR.cov, "expoT");

            expoFit_xDriftL->SetName(TString::Format("MPVfit_xDriftL_%dwires", i));
            expoFit_xDriftR->SetName(TString::Format("MPVfit_xDriftR_%dwires", i));
            expoFit_tDriftL->SetName(TString::Format("MPVfit_tDriftL_%dwires", i));
            expoFit_tDriftR->SetName(TString::Format("MPVfit_tDriftR_%dwires", i));
            
            expoFit_xDriftL->Write();
            expoFit_xDriftR->Write();
            expoFit_tDriftL->Write();
            expoFit_tDriftR->Write();

        

	}

    if(codeConfig == 2){

		fitExpoParameters expoParams_xDriftL;
		fitExpoParameters expoParams_xDriftR;
		fitExpoParameters expoParams_tDriftL;
		fitExpoParameters expoParams_tDriftR;

        TF1 *expoFit_xDriftL = new TF1();
		TF1 *expoFit_xDriftR = new TF1();
		TF1 *expoFit_tDriftL = new TF1();
		TF1 *expoFit_tDriftR = new TF1();

        ff.cd();
        
        for(int i = 1; i <= maxWireGroup; i++){

            TGraphAsymmErrors *MPVplot_xDrift_hits = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_%dhits", i));
            TGraphAsymmErrors *MPVplot_tDriftL_hits = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_%dhits", i));
            TGraphAsymmErrors *MPVplot_tDriftR_hits = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_%dhits", i));

            expoParams_xDriftL.reset();
            expoParams_xDriftR.reset();
            expoParams_tDriftL.reset();
            expoParams_tDriftR.reset();

            SetExpoParameters(expoParams_xDriftL.fp);
            SetExpoParameters(expoParams_xDriftR.fp);
            SetExpoParameters(expoParams_tDriftL.fp);
            SetExpoParameters(expoParams_tDriftR.fp);

            expoParams_xDriftL.lb = expo_xDriftL_lb;
            expoParams_xDriftL.ub = expo_xDriftL_ub;
            expoParams_xDriftR.lb = expo_xDriftR_lb;
            expoParams_xDriftR.ub = expo_xDriftR_ub;
            expoParams_tDriftL.lb = expo_tDriftL_lb;
            expoParams_tDriftL.ub = expo_tDriftL_ub;
            expoParams_tDriftR.lb = expo_tDriftR_lb;
            expoParams_tDriftR.ub = expo_tDriftR_ub;

            expoFit_xDriftL = fitter(MPVplot_xDrift_hits, expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.ehfp, expoParams_xDriftL.elfp, expoParams_xDriftL.cov, "expoX");
			expoFit_xDriftR = fitter(MPVplot_xDrift_hits, expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.ehfp, expoParams_xDriftR.elfp, expoParams_xDriftR.cov, "expoX");
			expoFit_tDriftL = fitter(MPVplot_tDriftL_hits, expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.ehfp, expoParams_tDriftL.elfp, expoParams_tDriftL.cov, "expoT");
			expoFit_tDriftR = fitter(MPVplot_tDriftR_hits, expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.ehfp, expoParams_tDriftR.elfp, expoParams_tDriftR.cov, "expoT");

            expoFit_xDriftL->SetName(TString::Format("MPVfit_xDriftL_%dhits", i));
            expoFit_xDriftR->SetName(TString::Format("MPVfit_xDriftR_%dhits", i));
            expoFit_tDriftL->SetName(TString::Format("MPVfit_tDriftL_%dhits", i));
            expoFit_tDriftR->SetName(TString::Format("MPVfit_tDriftR_%dhits", i));
            
            expoFit_xDriftL->Write();
            expoFit_xDriftR->Write();
            expoFit_tDriftL->Write();
            expoFit_tDriftR->Write();

        }

	}

    if(codeConfig == 3){

		fitExpoParameters expoParams_xDriftL_toW;
		fitExpoParameters expoParams_xDriftR_toW;
		fitExpoParameters expoParams_tDriftL_toW;
		fitExpoParameters expoParams_tDriftR_toW;

		fitExpoParameters expoParams_xDriftL_toE;
		fitExpoParameters expoParams_xDriftR_toE;
		fitExpoParameters expoParams_tDriftL_toE;
		fitExpoParameters expoParams_tDriftR_toE;

        TF1 *expoFit_xDriftL_toW = new TF1();
		TF1 *expoFit_xDriftR_toW = new TF1();
		TF1 *expoFit_tDriftL_toW = new TF1();
		TF1 *expoFit_tDriftR_toW = new TF1();

		TF1 *expoFit_xDriftL_toE = new TF1();
		TF1 *expoFit_xDriftR_toE = new TF1();
		TF1 *expoFit_tDriftL_toE = new TF1();
		TF1 *expoFit_tDriftR_toE = new TF1();

        double** angBinLimits = getAngBinLimits(angBins);

        ff.cd();
        
        for(int i = 0; i < angBins/2; i++){

            TGraphAsymmErrors *MPVplot_xDrift_toW= (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            TGraphAsymmErrors *MPVplot_tDriftL_toW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            TGraphAsymmErrors *MPVplot_tDriftR_toW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));

			TGraphAsymmErrors *MPVplot_xDrift_toE= (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1])); 
			TGraphAsymmErrors *MPVplot_tDriftL_toE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
			TGraphAsymmErrors *MPVplot_tDriftR_toE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));

            expoParams_xDriftL_toW.reset();
            expoParams_xDriftR_toW.reset();
            expoParams_tDriftL_toW.reset();
            expoParams_tDriftR_toW.reset();

			expoParams_xDriftL_toE.reset();
            expoParams_xDriftR_toE.reset();
            expoParams_tDriftL_toE.reset();
            expoParams_tDriftR_toE.reset();

            SetExpoParameters(expoParams_xDriftL_toW.fp);
            SetExpoParameters(expoParams_xDriftR_toW.fp);
            SetExpoParameters(expoParams_tDriftL_toW.fp);
            SetExpoParameters(expoParams_tDriftR_toW.fp);

			SetExpoParameters(expoParams_xDriftL_toE.fp);
            SetExpoParameters(expoParams_xDriftR_toE.fp);
            SetExpoParameters(expoParams_tDriftL_toE.fp);
            SetExpoParameters(expoParams_tDriftR_toE.fp);

            expoParams_xDriftL_toW.lb = expo_xDriftL_lb;
            expoParams_xDriftL_toW.ub = expo_xDriftL_ub;
            expoParams_xDriftR_toW.lb = expo_xDriftR_lb;
            expoParams_xDriftR_toW.ub = expo_xDriftR_ub;
            expoParams_tDriftL_toW.lb = expo_tDriftL_lb;
            expoParams_tDriftL_toW.ub = expo_tDriftL_ub;
            expoParams_tDriftR_toW.lb = expo_tDriftR_lb;
            expoParams_tDriftR_toW.ub = expo_tDriftR_ub;

			expoParams_xDriftL_toE.lb = expo_xDriftL_lb;
            expoParams_xDriftL_toE.ub = expo_xDriftL_ub;
            expoParams_xDriftR_toE.lb = expo_xDriftR_lb;
            expoParams_xDriftR_toE.ub = expo_xDriftR_ub;
            expoParams_tDriftL_toE.lb = expo_tDriftL_lb;
            expoParams_tDriftL_toE.ub = expo_tDriftL_ub;
            expoParams_tDriftR_toE.lb = expo_tDriftR_lb;
            expoParams_tDriftR_toE.ub = expo_tDriftR_ub;

            expoFit_xDriftL_toW = fitter(MPVplot_xDrift_toW, expoParams_xDriftL_toW.lb, expoParams_xDriftL_toW.ub, expoParams_xDriftL_toW.fp, expoParams_xDriftL_toW.ehfp, expoParams_xDriftL_toW.elfp, expoParams_xDriftL_toW.cov, "expoX");
			expoFit_xDriftR_toW = fitter(MPVplot_xDrift_toW, expoParams_xDriftR_toW.lb, expoParams_xDriftR_toW.ub, expoParams_xDriftR_toW.fp, expoParams_xDriftR_toW.ehfp, expoParams_xDriftR_toW.elfp, expoParams_xDriftR_toW.cov, "expoX");
			expoFit_tDriftL_toW = fitter(MPVplot_tDriftL_toW, expoParams_tDriftL_toW.lb, expoParams_tDriftL_toW.ub, expoParams_tDriftL_toW.fp, expoParams_tDriftL_toW.ehfp, expoParams_tDriftL_toW.elfp, expoParams_tDriftL_toW.cov, "expoT");
			expoFit_tDriftR_toW = fitter(MPVplot_tDriftR_toW, expoParams_tDriftR_toW.lb, expoParams_tDriftR_toW.ub, expoParams_tDriftR_toW.fp, expoParams_tDriftR_toW.ehfp, expoParams_tDriftR_toW.elfp, expoParams_tDriftR_toW.cov, "expoT");

			expoFit_xDriftL_toE = fitter(MPVplot_xDrift_toE, expoParams_xDriftL_toE.lb, expoParams_xDriftL_toE.ub, expoParams_xDriftL_toE.fp, expoParams_xDriftL_toE.ehfp, expoParams_xDriftL_toE.elfp, expoParams_xDriftL_toE.cov, "expoX");
			expoFit_xDriftR_toE = fitter(MPVplot_xDrift_toE, expoParams_xDriftR_toE.lb, expoParams_xDriftR_toE.ub, expoParams_xDriftR_toE.fp, expoParams_xDriftR_toE.ehfp, expoParams_xDriftR_toE.elfp, expoParams_xDriftR_toE.cov, "expoX");
			expoFit_tDriftL_toE = fitter(MPVplot_tDriftL_toE, expoParams_tDriftL_toE.lb, expoParams_tDriftL_toE.ub, expoParams_tDriftL_toE.fp, expoParams_tDriftL_toE.ehfp, expoParams_tDriftL_toE.elfp, expoParams_tDriftL_toE.cov, "expoT");
			expoFit_tDriftR_toE = fitter(MPVplot_tDriftR_toE, expoParams_tDriftR_toE.lb, expoParams_tDriftR_toE.ub, expoParams_tDriftR_toE.fp, expoParams_tDriftR_toE.ehfp, expoParams_tDriftR_toE.elfp, expoParams_tDriftR_toE.cov, "expoT");

            expoFit_xDriftL_toW->SetName(TString::Format("MPVfit_xDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            expoFit_xDriftR_toW->SetName(TString::Format("MPVfit_xDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            expoFit_tDriftL_toW->SetName(TString::Format("MPVfit_tDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            expoFit_tDriftR_toW->SetName(TString::Format("MPVfit_tDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
            
			expoFit_xDriftL_toE->SetName(TString::Format("MPVfit_xDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
            expoFit_xDriftR_toE->SetName(TString::Format("MPVfit_xDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
            expoFit_tDriftL_toE->SetName(TString::Format("MPVfit_tDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
            expoFit_tDriftR_toE->SetName(TString::Format("MPVfit_tDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));

            expoFit_xDriftL_toW->Write();
            expoFit_xDriftR_toW->Write();
            expoFit_tDriftL_toW->Write();
            expoFit_tDriftR_toW->Write();

			expoFit_xDriftL_toE->Write();
            expoFit_xDriftR_toE->Write();
            expoFit_tDriftL_toE->Write();
            expoFit_tDriftR_toE->Write();

        }

	}

    if(codeConfig == 4){

		fitExpoParameters expoParams_xDriftL_toW;
		fitExpoParameters expoParams_xDriftR_toW;
		fitExpoParameters expoParams_tDriftL_toW;
		fitExpoParameters expoParams_tDriftR_toW;

		fitExpoParameters expoParams_xDriftL_toE;
		fitExpoParameters expoParams_xDriftR_toE;
		fitExpoParameters expoParams_tDriftL_toE;
		fitExpoParameters expoParams_tDriftR_toE;

        TF1 *expoFit_xDriftL_toW = new TF1();
		TF1 *expoFit_xDriftR_toW = new TF1();
		TF1 *expoFit_tDriftL_toW = new TF1();
		TF1 *expoFit_tDriftR_toW = new TF1();

		TF1 *expoFit_xDriftL_toE = new TF1();
		TF1 *expoFit_xDriftR_toE = new TF1();
		TF1 *expoFit_tDriftL_toE = new TF1();
		TF1 *expoFit_tDriftR_toE = new TF1();

        double** angBinLimits = getAngBinLimits(angBins);

        ff.cd();
        
        for(int i = 0; i < angBins/2; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				TGraphAsymmErrors *MPVplot_xDrift_toW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				TGraphAsymmErrors *MPVplot_tDriftL_toW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				TGraphAsymmErrors *MPVplot_tDriftR_toW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));

				TGraphAsymmErrors *MPVplot_xDrift_toE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j)); 
				TGraphAsymmErrors *MPVplot_tDriftL_toE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));
				TGraphAsymmErrors *MPVplot_tDriftR_toE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));

				expoParams_xDriftL_toW.reset();
				expoParams_xDriftR_toW.reset();
				expoParams_tDriftL_toW.reset();
				expoParams_tDriftR_toW.reset();

				expoParams_xDriftL_toE.reset();
				expoParams_xDriftR_toE.reset();
				expoParams_tDriftL_toE.reset();
				expoParams_tDriftR_toE.reset();

				SetExpoParameters(expoParams_xDriftL_toW.fp);
				SetExpoParameters(expoParams_xDriftR_toW.fp);
				SetExpoParameters(expoParams_tDriftL_toW.fp);
				SetExpoParameters(expoParams_tDriftR_toW.fp);

				SetExpoParameters(expoParams_xDriftL_toE.fp);
				SetExpoParameters(expoParams_xDriftR_toE.fp);
				SetExpoParameters(expoParams_tDriftL_toE.fp);
				SetExpoParameters(expoParams_tDriftR_toE.fp);

				expoParams_xDriftL_toW.lb = expo_xDriftL_lb;
				expoParams_xDriftL_toW.ub = expo_xDriftL_ub;
				expoParams_xDriftR_toW.lb = expo_xDriftR_lb;
				expoParams_xDriftR_toW.ub = expo_xDriftR_ub;
				expoParams_tDriftL_toW.lb = expo_tDriftL_lb;
				expoParams_tDriftL_toW.ub = expo_tDriftL_ub;
				expoParams_tDriftR_toW.lb = expo_tDriftR_lb;
				expoParams_tDriftR_toW.ub = expo_tDriftR_ub;

				expoParams_xDriftL_toE.lb = expo_xDriftL_lb;
				expoParams_xDriftL_toE.ub = expo_xDriftL_ub;
				expoParams_xDriftR_toE.lb = expo_xDriftR_lb;
				expoParams_xDriftR_toE.ub = expo_xDriftR_ub;
				expoParams_tDriftL_toE.lb = expo_tDriftL_lb;
				expoParams_tDriftL_toE.ub = expo_tDriftL_ub;
				expoParams_tDriftR_toE.lb = expo_tDriftR_lb;
				expoParams_tDriftR_toE.ub = expo_tDriftR_ub;

				expoFit_xDriftL_toW = fitter(MPVplot_xDrift_toW, expoParams_xDriftL_toW.lb, expoParams_xDriftL_toW.ub, expoParams_xDriftL_toW.fp, expoParams_xDriftL_toW.ehfp, expoParams_xDriftL_toW.elfp, expoParams_xDriftL_toW.cov, "expoX");
				expoFit_xDriftR_toW = fitter(MPVplot_xDrift_toW, expoParams_xDriftR_toW.lb, expoParams_xDriftR_toW.ub, expoParams_xDriftR_toW.fp, expoParams_xDriftR_toW.ehfp, expoParams_xDriftR_toW.elfp, expoParams_xDriftR_toW.cov, "expoX");
				expoFit_tDriftL_toW = fitter(MPVplot_tDriftL_toW, expoParams_tDriftL_toW.lb, expoParams_tDriftL_toW.ub, expoParams_tDriftL_toW.fp, expoParams_tDriftL_toW.ehfp, expoParams_tDriftL_toW.elfp, expoParams_tDriftL_toW.cov, "expoT");
				expoFit_tDriftR_toW = fitter(MPVplot_tDriftR_toW, expoParams_tDriftR_toW.lb, expoParams_tDriftR_toW.ub, expoParams_tDriftR_toW.fp, expoParams_tDriftR_toW.ehfp, expoParams_tDriftR_toW.elfp, expoParams_tDriftR_toW.cov, "expoT");

				expoFit_xDriftL_toE = fitter(MPVplot_xDrift_toE, expoParams_xDriftL_toE.lb, expoParams_xDriftL_toE.ub, expoParams_xDriftL_toE.fp, expoParams_xDriftL_toE.ehfp, expoParams_xDriftL_toE.elfp, expoParams_xDriftL_toE.cov, "expoX");
				expoFit_xDriftR_toE = fitter(MPVplot_xDrift_toE, expoParams_xDriftR_toE.lb, expoParams_xDriftR_toE.ub, expoParams_xDriftR_toE.fp, expoParams_xDriftR_toE.ehfp, expoParams_xDriftR_toE.elfp, expoParams_xDriftR_toE.cov, "expoX");
				expoFit_tDriftL_toE = fitter(MPVplot_tDriftL_toE, expoParams_tDriftL_toE.lb, expoParams_tDriftL_toE.ub, expoParams_tDriftL_toE.fp, expoParams_tDriftL_toE.ehfp, expoParams_tDriftL_toE.elfp, expoParams_tDriftL_toE.cov, "expoT");
				expoFit_tDriftR_toE = fitter(MPVplot_tDriftR_toE, expoParams_tDriftR_toE.lb, expoParams_tDriftR_toE.ub, expoParams_tDriftR_toE.fp, expoParams_tDriftR_toE.ehfp, expoParams_tDriftR_toE.elfp, expoParams_tDriftR_toE.cov, "expoT");

				expoFit_xDriftL_toW->SetName(TString::Format("MPVfit_xDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				expoFit_xDriftR_toW->SetName(TString::Format("MPVfit_xDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				expoFit_tDriftL_toW->SetName(TString::Format("MPVfit_tDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				expoFit_tDriftR_toW->SetName(TString::Format("MPVfit_tDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j));
				
				expoFit_xDriftL_toE->SetName(TString::Format("MPVfit_xDriftL_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));
				expoFit_xDriftR_toE->SetName(TString::Format("MPVfit_xDriftR_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));
				expoFit_tDriftL_toE->SetName(TString::Format("MPVfit_tDriftL_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));
				expoFit_tDriftR_toE->SetName(TString::Format("MPVfit_tDriftR_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j));

				expoFit_xDriftL_toW->Write();
				expoFit_xDriftR_toW->Write();
				expoFit_tDriftL_toW->Write();
				expoFit_tDriftR_toW->Write();

				expoFit_xDriftL_toE->Write();
				expoFit_xDriftR_toE->Write();
				expoFit_tDriftL_toE->Write();
				expoFit_tDriftR_toE->Write();

			}

        }

	}

	if(codeConfig == 5){

		fitExpoParameters expoParams_xDriftL;
		fitExpoParameters expoParams_xDriftR;
		fitExpoParameters expoParams_tDriftL;
		fitExpoParameters expoParams_tDriftR;

		TF1 *expoFit_xDriftL = new TF1();
		TF1 *expoFit_xDriftR = new TF1();
		TF1 *expoFit_tDriftL = new TF1();
		TF1 *expoFit_tDriftR = new TF1();

		ff.cd();
        
        for(int i = 1; i <= 4; i++){

            TGraphAsymmErrors *MPVplot_xDriftL = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDriftL_oct%i",2*i));
			TGraphAsymmErrors *MPVplot_xDriftR = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDriftR_oct%i",(2*i)-1));
			TGraphAsymmErrors *MPVplot_tDriftL = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftL_oct%i",2*i));
            TGraphAsymmErrors *MPVplot_tDriftR = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftR_oct%i",(2*i)-1));
            
			expoParams_xDriftL.reset();
            expoParams_xDriftR.reset();
			expoParams_tDriftL.reset();
            expoParams_tDriftR.reset();
            
			SetExpoParameters(expoParams_xDriftL.fp);
            SetExpoParameters(expoParams_xDriftR.fp);
			SetExpoParameters(expoParams_tDriftL.fp);
            SetExpoParameters(expoParams_tDriftR.fp);
            
			expoParams_xDriftL.lb = expo_xDriftL_lb;
            expoParams_xDriftL.ub = expo_xDriftL_ub;
            expoParams_xDriftR.lb = expo_xDriftR_lb;
            expoParams_xDriftR.ub = expo_xDriftR_ub;
			expoParams_tDriftL.lb = expo_tDriftL_lb;
            expoParams_tDriftL.ub = expo_tDriftL_ub;
            expoParams_tDriftR.lb = expo_tDriftR_lb;
            expoParams_tDriftR.ub = expo_tDriftR_ub;
            
			expoFit_xDriftL = fitter(MPVplot_xDriftL, expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.ehfp, expoParams_xDriftL.elfp, expoParams_xDriftL.cov, "expoX");
			expoFit_xDriftR = fitter(MPVplot_xDriftR, expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.ehfp, expoParams_xDriftR.elfp, expoParams_xDriftR.cov, "expoX");
			expoFit_tDriftL = fitter(MPVplot_tDriftL, expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.ehfp, expoParams_tDriftL.elfp, expoParams_tDriftL.cov, "expoT");
			expoFit_tDriftR = fitter(MPVplot_tDriftR, expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.ehfp, expoParams_tDriftR.elfp, expoParams_tDriftR.cov, "expoT");

			
			expoFit_xDriftL->SetName(TString::Format("MPVfit_xDriftL_oct%i",2*i));
            expoFit_xDriftR->SetName(TString::Format("MPVfit_xDriftR_oct%i",(2*i)-1));
			expoFit_tDriftL->SetName(TString::Format("MPVfit_tDriftL_oct%i",2*i));
            expoFit_tDriftR->SetName(TString::Format("MPVfit_tDriftR_oct%i",(2*i)-1));
            
			expoFit_xDriftL->Write();
            expoFit_xDriftR->Write();
			expoFit_tDriftL->Write();
            expoFit_tDriftR->Write();
            
        }

	}

	return 0;

}