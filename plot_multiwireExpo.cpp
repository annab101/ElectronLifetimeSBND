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
	std::string id = "noID"; //should be X or T then L or R e.g XR for distance and TPC West

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
		if(!strcmp(argv[i], "--id")){
			id = argv[i+1];
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

	std::string TPCLabel = "";
	int color = 1;
	int marker;
	if (id.find("L") != std::string::npos) {
    	TPCLabel = "East TPC";
		color = coral;
	}
	else if (id.find("R") != std::string::npos) {
    	TPCLabel = "West TPC";
		color = kAzure - 3;
	}
	else {
    	std::cout << "id doesn't contain L or R - won't be able to find the graphs in the root file...";
	}

	if (id.find("X") != std::string::npos) {
    	marker = 21;
	}
	else if (id.find("T") != std::string::npos) {
    	marker = 47;
	}
    
	TFile f((saveLoc + dataset  + "_" + configLabel + "/multiwireExpos_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

	TGraphErrors *g = (TGraphErrors*)f.Get(("g_lifeVwires" + id).c_str());

	TF1* expoFit = (TF1*)f.Get(("expoFit_" + id).c_str());
	
	setMarker(g, color, marker, 1.0);

	g->SetTitle(";N; #tau (ms)");
	
	setFontSize<TGraphErrors>(g, 133, 25);

	TPaveText *pt = new TPaveText(.45,.65,.80,.88, "blNDC");

	pt->SetBorderSize(1);
	pt->SetFillColor(0);
	pt->AddText(TPCLabel.c_str());
    ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
	pt->AddText(Form("#chi^{2} / DoF: %g / %i",expoFit->GetChisquare(),expoFit->GetNDF()));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(color);
	pt->AddText(Form("N: %g#pm%g",expoFit->GetParameter(0),expoFit->GetParError(0)));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(color);
	pt->AddText(Form("#lambda: %g#pm%g",expoFit->GetParameter(1),expoFit->GetParError(1)));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(color);
	pt->AddText(Form("#tau_{0}: %g#pm%g",expoFit->GetParameter(2),expoFit->GetParError(2)));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(color);


	TCanvas *c = new TCanvas();
	c->SetLeftMargin(0.15);
	c->SetRightMargin(0.18);
	c->SetBottomMargin(0.12);
	g->Draw("AP");
	pt->Draw();

	gPad->Update();
    TF1 *f1 = (TF1*)g->GetListOfFunctions()->FindObject("Fitfcn_");
    f1->SetLineColor(1);

	c->DrawClone();

	saveFig(c, saveLoc + dataset  + "_" + configLabel + "/plot_multiwireExpo" + id + "_" + tag);

    return 0;

}