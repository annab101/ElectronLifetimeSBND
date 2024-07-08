//Author: Anna Beever
//Date:   April 2024

//C++ includes
#include <iomanip>

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
    std::string histName = "noHist";
	int binNum = -1;

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
        if(!strcmp(argv[i], "--histName")){
			histName = argv[i+1];
		}
		if(!strcmp(argv[i], "--binNum")){
			binNum = std::stoi(argv[i+1]);
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

    p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);

    configLabel = configMap[codeConfig];

    std::cout << configLabel << std::endl;
    
	TFile f_stats((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());
    TFile f_MPV((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

	TH2D* h_stats = (TH2D*)f_stats.Get("h_stats");
    TH1D* projY = (TH1D*)f_MPV.Get(("projY_xDrift_" + histName + "_bin" + std::to_string(binNum)).c_str());
	TF1* LGfit = (TF1*)f_MPV.Get(("LGfit_xDrift_" + histName + "_bin" + std::to_string(binNum)).c_str());

	std::string projYtitle;
	std::stringstream s1;
	s1 << std::setprecision(4) << h_stats->GetXaxis()->GetBinLowEdge(binNum);
	std::string str1 = s1.str();
	std::stringstream s2;
	s2 << std::setprecision(4) << (h_stats->GetXaxis()->GetBinLowEdge(binNum) + h_stats->GetXaxis()->GetBinWidth(binNum));
	std::string str2 = s2.str();
	projYtitle = "dQ/dx for " + str1 + " cm < x < " + str2 + " cm;dQ/dx (ADC/cm);Entries";
	projY->SetTitle(projYtitle.c_str());
	projY->SetTitleFont(133);
	projY->SetTitleSize(28);

	projY->GetXaxis()->SetLabelOffset(0.01);
	projY->GetXaxis()->SetNdivisions(505);
	setFontSize<TH1>(projY, 133, 25);
	projY->GetYaxis()->SetLabelOffset(0.01);

	projY->SetFillColor(flamingo);
	projY->SetLineColor(flamingo);
	projY->SetStats(0);
	LGfit->SetLineColor(deepViolet);
	LGfit->SetLineWidth(2.3);
 
	TCanvas *c_plain = new TCanvas();
	c_plain->SetLeftMargin(0.15);
	c_plain->SetBottomMargin(0.12);
	projY->Draw();
	LGfit->Draw("same");
	saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/plot_projY_" + histName + tag);

	//Cheating for now
	fitLGParameters fp_xDrift;
	SetLGParameters(projY, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
	LGfit = fitter(projY, fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
	std::cout << "MPV: " << fp_xDrift.fp[1] << std::endl;
	std::cout << "eMPV: " << pointError(projY, fp_xDrift) << std::endl;
	std::cout << "scale: " << fp_xDrift.fp[0] << std::endl;
	std::cout << "escale: " << fp_xDrift.efp[0] << std::endl;
	std::cout << "norm: " << fp_xDrift.fp[2] << std::endl;
	std::cout << "enorm: " << fp_xDrift.efp[2] << std::endl;
	std::cout << "sigma: " << fp_xDrift.fp[3] << std::endl;
	std::cout << "esigma: " << fp_xDrift.efp[3] << std::endl;

    return 0;

}