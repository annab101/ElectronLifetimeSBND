//Author: Anna Beever
//Date:   April 2024

//C++ includes

//ROOT includes
#include "Utilities/ROOTincludes.h"

//Local includes
#include "Helpers/Constants.h"
#include "Helpers/PlottingHelpers.h"
#include "Helpers/PlottingHelpers.cpp"
#include "Utilities/ConfigReader.h"
#include "Utilities/ConfigReader.cpp"
#include "Helpers/FittingHelpers.h"
#include "Helpers/FittingHelpers.cpp"
#include "Helpers/StatsHelpers.h"
#include "Helpers/StatsHelpers.cpp"

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

    int codeConfig = 0;
    std::string tag = "noTag";
	std::string dataset = "noDataSet";
	std::string saveLoc = "noSaveLoc";
    std::string configLabel = "noConfigLabel";
	int maxWireGroup = 20;

	p->getValue("maxWireGroup", maxWireGroup);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);

    configLabel = configMap[codeConfig];

    std::cout << configLabel << std::endl;
    
    TFile f_fit((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

	TGraphErrors *g_lifeVwiresXL = new TGraphErrors;
	TGraphErrors *g_lifeVwiresXR = new TGraphErrors;
	TGraphErrors *g_lifeVwiresTL = new TGraphErrors;
	TGraphErrors *g_lifeVwiresTR = new TGraphErrors;

	for(int i = 1; i <= maxWireGroup; i++){

		TF1* expoFit_xDriftL = (TF1*)f_fit.Get(("MPVfit_xDriftL_" + std::to_string(i) + "wires").c_str());
		TF1* expoFit_xDriftR = (TF1*)f_fit.Get(("MPVfit_xDriftR_" + std::to_string(i) + "wires").c_str());
		TF1* expoFit_tDriftL = (TF1*)f_fit.Get(("MPVfit_tDriftL_" + std::to_string(i) + "wires").c_str());
		TF1* expoFit_tDriftR = (TF1*)f_fit.Get(("MPVfit_tDriftR_" + std::to_string(i) + "wires").c_str());

		if(expoFit_xDriftL->GetParameter(1) >= 0. && expoFit_xDriftL->GetParameter(1) <= 50. && expoFit_xDriftL->GetParError(1) <= 50){
			g_lifeVwiresXL->SetPoint(g_lifeVwiresXL->GetN(), i, expoFit_xDriftL->GetParameter(1));
			g_lifeVwiresXL->SetPointError(g_lifeVwiresXL->GetN()-1, 0., expoFit_xDriftL->GetParError(1));
		}
		if(expoFit_xDriftR->GetParameter(1) >= 0. && expoFit_xDriftR->GetParameter(1) <= 50. && expoFit_xDriftR->GetParError(1) <= 50){
			g_lifeVwiresXR->SetPoint(g_lifeVwiresXR->GetN(), i, expoFit_xDriftR->GetParameter(1));
			g_lifeVwiresXR->SetPointError(g_lifeVwiresXR->GetN()-1, 0., expoFit_xDriftR->GetParError(1));
		}
		if(expoFit_tDriftL->GetParameter(1) >= 0. && expoFit_tDriftL->GetParameter(1) <= 50. && expoFit_tDriftL->GetParError(1) <= 50){
			g_lifeVwiresTL->SetPoint(g_lifeVwiresTL->GetN(), i, expoFit_tDriftL->GetParameter(1));
			g_lifeVwiresTL->SetPointError(g_lifeVwiresTL->GetN()-1, 0., expoFit_tDriftL->GetParError(1));
		}
		if(expoFit_tDriftR->GetParameter(1) >= 0. && expoFit_tDriftR->GetParameter(1) <= 50. && expoFit_tDriftR->GetParError(1) <= 50){
			g_lifeVwiresTR->SetPoint(g_lifeVwiresTR->GetN(), i, expoFit_tDriftR->GetParameter(1));
			g_lifeVwiresTR->SetPointError(g_lifeVwiresTR->GetN()-1, 0., expoFit_tDriftR->GetParError(1));
		}

	}

	TCanvas *c_plain = new TCanvas();
	c_plain->SetLeftMargin(0.15);
	c_plain->SetRightMargin(0.1);
	c_plain->SetBottomMargin(0.15);
	TMultiGraph *mg = new TMultiGraph();
	setMarker(g_lifeVwiresXL, coral, 21, 1.0);
	setMarker(g_lifeVwiresXR, kAzure - 3, 21, 1.0);
	setMarker(g_lifeVwiresTL, coral, 47, 1.0);
	setMarker(g_lifeVwiresTR, kAzure - 3, 47, 1.0);

	g_lifeVwiresXL->SetTitle("x, East TPC");
	g_lifeVwiresXR->SetTitle("x, West TPC");
	g_lifeVwiresTL->SetTitle("t, East TPC");
	g_lifeVwiresTR->SetTitle("t, West TPC");

	mg->SetTitle(";N; #tau (ms)");
	mg->Add(g_lifeVwiresXL);
	mg->Add(g_lifeVwiresXR);
	mg->Add(g_lifeVwiresTL);
	mg->Add(g_lifeVwiresTR);
	setFontSize<TMultiGraph>(mg, 133, 25);
	mg->Draw("AP");
	c_plain->BuildLegend(0.4,0.7,0.6,0.88);
	saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/plot_multiwire_" + tag);

	TFile f00((saveLoc + dataset  + "_" + configLabel + "/multiwireExpos_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
	f00.cd();

	fitExpoConstParameters expoXL;
	fitExpoConstParameters expoXR;
	fitExpoConstParameters expoTL;
	fitExpoConstParameters expoTR;

	SetExpoConstParameters(expoXL.fp);
	SetExpoConstParameters(expoXR.fp);
	SetExpoConstParameters(expoTL.fp);
	SetExpoConstParameters(expoTR.fp);

	TF1 *expoFitXL = new TF1();
	TF1 *expoFitXR = new TF1();
	TF1 *expoFitTL = new TF1();
	TF1 *expoFitTR = new TF1();

	expoFitXL = fitter(g_lifeVwiresXL, expoXL.lb, expoXL.ub, expoXL.fp, expoXL.efp, expoXL.cov, "expoConst");
	expoFitXR = fitter(g_lifeVwiresXR, expoXR.lb, expoXR.ub, expoXR.fp, expoXR.efp, expoXR.cov, "expoConst");
	expoFitTL = fitter(g_lifeVwiresTL, expoTL.lb, expoTL.ub, expoTL.fp, expoTL.efp, expoTL.cov, "expoConst");
	expoFitTR = fitter(g_lifeVwiresTR, expoTR.lb, expoTR.ub, expoTR.fp, expoTR.efp, expoTR.cov, "expoConst");

	g_lifeVwiresXL->SetName("g_lifeVwiresXL");
	g_lifeVwiresXR->SetName("g_lifeVwiresXR");
	g_lifeVwiresTL->SetName("g_lifeVwiresTL");
	g_lifeVwiresTR->SetName("g_lifeVwiresTR");

	g_lifeVwiresXL->Write();
	g_lifeVwiresXR->Write();
	g_lifeVwiresTL->Write();
	g_lifeVwiresTR->Write();

	expoFitXL->SetName("expoFit_XL");
	expoFitXR->SetName("expoFit_XR");
	expoFitTL->SetName("expoFit_TL");
	expoFitTR->SetName("expoFit_TR");

	expoFitXL->Write();
	expoFitXR->Write();
	expoFitTL->Write();
	expoFitTR->Write();

    return 0;

}