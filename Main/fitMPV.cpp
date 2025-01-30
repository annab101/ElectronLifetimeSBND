//Author: Anna Beever
//Date:   April 2024

//C++ includes
#include<fstream>
#include<iostream>
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
	int nGroupedWires = 1;
	double expo_xDriftE_lb = -160.;
    double expo_xDriftW_lb = 40.;
    double expo_xDriftE_ub = -40.;
    double expo_xDriftW_ub = 160.;
    double expo_tDriftE_lb = 0.2;
    double expo_tDriftW_lb = 0.2;
    double expo_tDriftE_ub = 1.0;
    double expo_tDriftW_ub = 1.0;

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
    p->getValue("nGroupedWires", nGroupedWires);
    p->getValue("expo_xDriftE_lb", expo_xDriftE_lb);
	p->getValue("expo_xDriftW_lb", expo_xDriftW_lb);
	p->getValue("expo_xDriftE_ub", expo_xDriftE_ub);
	p->getValue("expo_xDriftW_ub", expo_xDriftW_ub);
	p->getValue("expo_tDriftE_lb", expo_tDriftE_lb);
	p->getValue("expo_tDriftW_lb", expo_tDriftW_lb);
	p->getValue("expo_tDriftE_ub", expo_tDriftE_ub);
	p->getValue("expo_tDriftW_ub", expo_tDriftW_ub);

	TFile f((saveLoc + dataset  + "/MPV_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str());
	TFile ff((saveLoc + dataset  + "/MPVfit_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str(), "new");
	std::string textFile = saveLoc + dataset  + "/LifetimeValues_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".txt";

	TGraphAsymmErrors *MPVplot_xDrift = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_xDrift_%dwires",nGroupedWires));
	TGraphAsymmErrors *MPVplot_tDriftE = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftE_%dwires",nGroupedWires));
	TGraphAsymmErrors *MPVplot_tDriftW = (TGraphAsymmErrors*)f.Get(TString::Format("MPVplot_tDriftW_%dwires",nGroupedWires));

	ff.cd();

	std::cout << "-----Initialising parameters-----" << std::endl;

	fitExpoParameters expoParams_xDriftE;
	fitExpoParameters expoParams_xDriftW;
	fitExpoParameters expoParams_tDriftE;
	fitExpoParameters expoParams_tDriftW;

	TF1 *expoFit_xDriftE = new TF1();
	TF1 *expoFit_xDriftW = new TF1();
	TF1 *expoFit_tDriftE = new TF1();
	TF1 *expoFit_tDriftW = new TF1();

	expoParams_xDriftE.reset();
	expoParams_xDriftW.reset();
	expoParams_tDriftE.reset();
	expoParams_tDriftW.reset();

	SetExpoParameters(expoParams_xDriftE.fp);
	SetExpoParameters(expoParams_xDriftW.fp);
	SetExpoParameters(expoParams_tDriftE.fp);
	SetExpoParameters(expoParams_tDriftW.fp);

	expoParams_xDriftE.lb = expo_xDriftE_lb;
	expoParams_xDriftE.ub = expo_xDriftE_ub;
	expoParams_xDriftW.lb = expo_xDriftW_lb;
	expoParams_xDriftW.ub = expo_xDriftW_ub;
	expoParams_tDriftE.lb = expo_tDriftE_lb;
	expoParams_tDriftE.ub = expo_tDriftE_ub;
	expoParams_tDriftW.lb = expo_tDriftW_lb;
	expoParams_tDriftW.ub = expo_tDriftW_ub;

	std::cout << "-----Fitting exponential-----" << std::endl;

	expoFit_xDriftE = fitter(MPVplot_xDrift, expoParams_xDriftE.lb, expoParams_xDriftE.ub, expoParams_xDriftE.fp, expoParams_xDriftE.ehfp, expoParams_xDriftE.elfp, expoParams_xDriftE.cov, "expoX");
	expoFit_xDriftW = fitter(MPVplot_xDrift, expoParams_xDriftW.lb, expoParams_xDriftW.ub, expoParams_xDriftW.fp, expoParams_xDriftW.ehfp, expoParams_xDriftW.elfp, expoParams_xDriftW.cov, "expoX");
	expoFit_tDriftE = fitter(MPVplot_tDriftE, expoParams_tDriftE.lb, expoParams_tDriftE.ub, expoParams_tDriftE.fp, expoParams_tDriftE.ehfp, expoParams_tDriftE.elfp, expoParams_tDriftE.cov, "expoT");
	expoFit_tDriftW = fitter(MPVplot_tDriftW, expoParams_tDriftW.lb, expoParams_tDriftW.ub, expoParams_tDriftW.fp, expoParams_tDriftW.ehfp, expoParams_tDriftW.elfp, expoParams_tDriftW.cov, "expoT");

	expoFit_xDriftE->SetName("MPVfit_xDriftE");
	expoFit_xDriftW->SetName("MPVfit_xDriftW");
	expoFit_tDriftE->SetName("MPVfit_tDriftE");
	expoFit_tDriftW->SetName("MPVfit_tDriftW");

	std::cout << "-----Writing to file-----" << std::endl;
	
	expoFit_xDriftE->Write();
	expoFit_xDriftW->Write();
	expoFit_tDriftE->Write();
	expoFit_tDriftW->Write();

	ofstream outFile;
	outFile.open(textFile.c_str());
	outFile << std::setw(20) << "TPC, dist/time  |  "
			<< std::setw(20) << "etime  |  "
			<< std::setw(20) << "error low  |  "
			<< std::setw(20) << "error up  | \n"
			<< std::setw(20) << "E, t  |  "
			<< std::setw(15) << expoParams_tDriftE.fp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftE.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftE.ehfp[1] << "  | \n"
			<< std::setw(20) << "W, t  |  "
			<< std::setw(15) << expoParams_tDriftW.fp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftW.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_tDriftW.ehfp[1] << "  | \n"
			<< std::setw(20) << "E, x  |  "
			<< std::setw(15) << expoParams_xDriftE.fp[1] << "  |  "
			<< std::setw(15) << expoParams_xDriftE.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_xDriftE.ehfp[1] << "  | \n"
			<< std::setw(20) << "W, x  |  "
			<< std::setw(15) << expoParams_xDriftW.fp[1] << "  |  "
			<< std::setw(15) << expoParams_xDriftW.elfp[1] << "  |  "
			<< std::setw(15) << expoParams_xDriftW.ehfp[1] << "  | \n";
	outFile.close();

	return 0;

}