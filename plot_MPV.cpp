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

using namespace calib;
using namespace constants;
using namespace cppsecrets;

int main(int argc, char**argv) {
	
	
	std::string filename = "noFile";
    std::string histName = "noHist";
    int distOrTime = 0; //0 = distance, 1 = time

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
        if(!strcmp(argv[i], "--histName")){
			histName = argv[i+1];
		}
        if(!strcmp(argv[i], "--distOrTime")){
			distOrTime = std::stoi(argv[i+1]);
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
    TFile f_fit((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

    TH2D* h_stats = (TH2D*)f_stats.Get("h_stats");
    int track_count = h_stats->GetEntries();

    if(distOrTime == 0){

        TGraphErrors* MPVplot = (TGraphErrors*)f_MPV.Get(("MPVplot_xDrift_" + histName).c_str());
        TF1* expoFitL = (TF1*)f_fit.Get(("MPVfit_xDriftL_" + histName).c_str());
        TF1* expoFitR = (TF1*)f_fit.Get(("MPVfit_xDriftR_" + histName).c_str());

        TCanvas *c_plain = new TCanvas();
        c_plain->SetLeftMargin(0.15);
        c_plain->SetRightMargin(0.18);
        c_plain->SetBottomMargin(0.12);

        setMarker<TGraphErrors>(MPVplot, kAzure - 3, 21, 0.6);
		MPVplot->SetTitle(";x (cm);dQ/dx MPV (ADC/cm)");
		MPVplot->GetXaxis()->SetNdivisions(505);
		setFontSize<TGraphErrors>(MPVplot, 133, 25);
		MPVplot->Draw("AP");
		expoFitL->SetLineColor(coral);
		expoFitL->SetLineWidth(2.3);
		expoFitL->Draw("SAME");
		expoFitR->SetLineColor(deepViolet);
		expoFitR->SetLineWidth(2.3);
		expoFitR->Draw("SAME");
		TPaveText *stats = statsBox({.34,.62,.63,.88}, track_count, expoFitL, expoFitR);
		stats->Draw();

		saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/plot_MPV_xDrift_" + histName + tag);

    }

    if(distOrTime == 1){

        TGraphErrors* MPVplotL = (TGraphErrors*)f_MPV.Get(("MPVplot_tDriftL_" + histName).c_str());
        TGraphErrors* MPVplotR = (TGraphErrors*)f_MPV.Get(("MPVplot_tDriftR_" + histName).c_str());
        TF1* expoFitL = (TF1*)f_fit.Get(("MPVfit_tDriftL_" + histName).c_str());
        TF1* expoFitR = (TF1*)f_fit.Get(("MPVfit_tDriftR_" + histName).c_str());

        TCanvas *c_split = new TCanvas();
        c_split->SetWindowSize(2000,500);
        TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
        TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
        c_split->cd();
        pad1->SetRightMargin(0.20);
        pad1->SetLeftMargin(0.18);
        pad2->SetRightMargin(0.20);
        pad2->SetLeftMargin(0.18);

        setMarker<TGraphErrors>(MPVplotL, kAzure - 3, 21, 0.6);
		MPVplotL->SetTitle("dQ/dx MPV vs t Left TPC;t (ms); dQ/dx MPV (ADC/cm)");
		expoFitL->SetLineColor(coral);
		setFontSize<TGraphErrors>(MPVplotL, 133, 25);
		TPaveText *statsL = statsBox({.45,.65,.80,.88}, track_count, expoFitL);
		setMarker<TGraphErrors>(MPVplotR, kAzure - 3, 21, 0.5);
		MPVplotR->SetTitle("dQ/dx MPV vs t Right TPC;t (ms); dQ/dx MPV (ADC/cm)");
		expoFitR->SetLineColor(deepViolet);
		setFontSize<TGraphErrors>(MPVplotR, 133, 25);
		TPaveText *statsR = statsBox({.45,.65,.80,.88}, track_count, expoFitR);
		
		c_split->cd();		
		pad1->Draw();		
		pad1->cd();	
		MPVplotL->Draw("AP");	
		expoFitL->Draw("SAME");
		statsL->Draw();
		c_split->cd();
		pad2->Draw();
		pad2->cd();
		MPVplotR->Draw("AP");
		expoFitR->Draw("SAME");
		statsR->Draw();
        saveFig(c_split, saveLoc + dataset  + "_" + configLabel + "/plot_MPV_tDrift_" + histName + tag);

    }

    return 0;

}