//Author: Anna Beever
//Date:   April 2024

//C++ includes

//ROOT includes
#include "../Utilities/ROOTincludes.h"

//Local includes
#include "../Helpers/Constants.h"
#include "../Helpers/PlottingHelpers.h"
#include "../Helpers/PlottingHelpers.cpp"
#include "../Utilities/ConfigReader.h"
#include "../Utilities/ConfigReader.cpp"

using namespace calib;
using namespace constants;
using namespace cppsecrets;

int main(int argc, char**argv) {
	
	std::string filename = "noFile";
    std::string histName = "noHist";
    int distOrTime = 0; //0 = distance, 1 = time
	int plotFit = 0;
	int plotStats = 0;
	int MPVmin = 940;
	int MPVmax = 1060;
	int col = 0;

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
		if(!strcmp(argv[i], "--plotFit")){
			plotFit = std::stoi(argv[i+1]);
		}
		if(!strcmp(argv[i], "--plotStats")){
			plotStats = std::stoi(argv[i+1]);
		}
		if(!strcmp(argv[i], "--MPVmin")){
			MPVmin = std::stoi(argv[i+1]);
		}
		if(!strcmp(argv[i], "--MPVmax")){
			MPVmax = std::stoi(argv[i+1]);
		}
		if(!strcmp(argv[i], "--col")){
			col = std::stoi(argv[i+1]);
		}
	}

	ConfigReader* p = ConfigReader::getInstance();
	p->parseFile(filename);
	std::cout << " Variables from configuration file: " << std::endl;
	p->dumpFileValues();
	std::cout << "-----------------------------------------------------------" << std::endl;

    int codeConfig = 0;
	int nGroupedWires = 0;
    std::string tag = "noTag";
	std::string dataset = "noDataSet";
	std::string saveLoc = "noSaveLoc";
    std::string configLabel = "noConfigLabel";
	std::string plotFitLabel = "noFitLabel";

    p->getValue("tag", tag);
	p->getValue("nGroupedWires", nGroupedWires);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);

    configLabel = configMap[codeConfig];
	plotFitLabel = plotFitName[plotFit];
    
    TFile f_stats((saveLoc + dataset + "/dQdx_hist_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str());
    TFile f_MPV((saveLoc + dataset + "/MPV_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str());
    TFile f_fit((saveLoc + dataset + "/MPVfit_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str());

    TH2D* h_stats = (TH2D*)f_stats.Get("h_stats");
    int track_count = h_stats->GetEntries();

	if(distOrTime == 0){

        TGraphAsymmErrors* MPVplot = (TGraphAsymmErrors*)f_MPV.Get(("MPVplot_xDrift_" + histName).c_str());
        TF1* expoFitL = (TF1*)f_fit.Get("MPVfit_xDriftE");
        TF1* expoFitR = (TF1*)f_fit.Get("MPVfit_xDriftW");

		TFitResult* expoFitLResult = (TFitResult*)f_fit.Get("MPVfit_xDriftE_fitResult");
        TFitResult* expoFitRResult = (TFitResult*)f_fit.Get("MPVfit_xDriftW_fitResult");

		TCanvas *c_plain = new TCanvas();
        c_plain->SetLeftMargin(0.15);
        c_plain->SetRightMargin(0.18);
        c_plain->SetBottomMargin(0.12);
		if(col){
			c_plain->SetFillColor(powderBlue);
			c_plain->SetFillStyle(1001);
			c_plain->SetFrameLineColor(deepViolet);
		}

		setMarker<TGraphAsymmErrors>(MPVplot, 1, 20, 0.7);
		MPVplot->SetTitle(";x (cm);dQ/dx MPV (Arb. Units)");
		MPVplot->SetMinimum(MPVmin);
    	MPVplot->SetMaximum(MPVmax);
		MPVplot->GetXaxis()->SetNdivisions(505);
		if(col){
			MPVplot->GetXaxis()->SetAxisColor(deepViolet);
			MPVplot->GetXaxis()->SetLabelColor(deepViolet);
			MPVplot->GetXaxis()->SetTitleColor(deepViolet);
			MPVplot->GetYaxis()->SetAxisColor(deepViolet);
			MPVplot->GetYaxis()->SetLabelColor(deepViolet);
			MPVplot->GetYaxis()->SetTitleColor(deepViolet);
		}
		setFontSize<TGraphAsymmErrors>(MPVplot, 133, 30);
		MPVplot->Draw("AP");
		expoFitL->SetLineColor(coral);
		expoFitL->SetLineWidth(4.0);
		if(plotFit){expoFitL->Draw("SAME");}
		expoFitR->SetLineColor(kAzure - 3);
		expoFitR->SetLineWidth(4.0);
		if(plotFit){expoFitR->Draw("SAME");}
		if(plotStats){
			if(plotFit){TPaveText *stats = statsBox({.34,.62,.63,.88}, track_count, expoFitL, expoFitR, expoFitLResult, expoFitRResult);stats->Draw();}
			if(!plotFit){TPaveText *stats = statsBox({.34,.8,.63,.88}, track_count);stats->Draw();}
		}

		if(col){saveFig(c_plain, saveLoc + dataset + "/plot_MPV_colour_xDrift_" + histName + "_" + plotFitLabel + "_" + tag);}
        if(!col){saveFig(c_plain, saveLoc + dataset + "/plot_MPV_xDrift_" + histName + "_" + plotFitLabel + "_" + tag);}

    }

    if(distOrTime == 1){

        TGraphAsymmErrors* MPVplotL = (TGraphAsymmErrors*)f_MPV.Get(("MPVplot_tDriftE_" + histName).c_str());
        TGraphAsymmErrors* MPVplotR = (TGraphAsymmErrors*)f_MPV.Get(("MPVplot_tDriftW_" + histName).c_str());
        TF1* expoFitL = (TF1*)f_fit.Get("MPVfit_tDriftE");
        TF1* expoFitR = (TF1*)f_fit.Get("MPVfit_tDriftW");
		TFitResult* expoFitLResult = (TFitResult*)f_fit.Get("MPVfit_tDriftE_fitResult");
        TFitResult* expoFitRResult = (TFitResult*)f_fit.Get("MPVfit_tDriftW_fitResult");

		TCanvas *c_split = new TCanvas();
        c_split->SetWindowSize(2000,500);
		if(col){
			c_split->SetFillColor(powderBlue);
			c_split->SetFillStyle(1001);
			c_split->SetFrameLineColor(deepViolet);
		}
        TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
        TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
        c_split->cd();
        pad1->SetRightMargin(0.15);
        pad1->SetLeftMargin(0.18);
		pad1->SetTopMargin(0.05);
		pad1->SetBottomMargin(0.15);
        pad2->SetRightMargin(0.15);
        pad2->SetLeftMargin(0.18);
		pad2->SetTopMargin(0.05);
		pad2->SetBottomMargin(0.15);
		if(col){
			pad1->SetFillColor(powderBlue);
			pad1->SetFillStyle(4100);
			pad1->SetFrameFillColor(powderBlue);
			pad1->SetFrameLineColor(deepViolet);
			pad1->SetFrameFillStyle(4100);
			pad2->SetFillColor(powderBlue);
			pad2->SetFillStyle(4100);
			pad2->SetFrameFillColor(powderBlue);
			pad2->SetFrameLineColor(deepViolet);
			pad2->SetFrameFillStyle(4100);
		}

		setMarker<TGraphAsymmErrors>(MPVplotL, 1, 20, 0.7);
		MPVplotL->SetTitle(";t_{#scale[1.2]{drift}} (ms); dQ/dx MPV (Arb. Units)");
		MPVplotL->SetMinimum(MPVmin);
    	MPVplotL->SetMaximum(MPVmax);

		if(col){
			MPVplotL->GetXaxis()->SetAxisColor(deepViolet);
			MPVplotL->GetXaxis()->SetLabelColor(deepViolet);
			MPVplotL->GetXaxis()->SetTitleColor(deepViolet);
			MPVplotL->GetYaxis()->SetAxisColor(deepViolet);
			MPVplotL->GetYaxis()->SetLabelColor(deepViolet);
			MPVplotL->GetYaxis()->SetTitleColor(deepViolet);
		}
		expoFitL->SetLineColor(coral);
		expoFitL->SetLineWidth(4.0);

		setFontSize<TGraphAsymmErrors>(MPVplotL, 133, 30);

		setMarker<TGraphAsymmErrors>(MPVplotR, 1, 20, 0.7);
		MPVplotR->SetTitle(";t_{#scale[1.2]{drift}} (ms); dQ/dx MPV (Arb. Units)");
		MPVplotR->SetMinimum(MPVmin);
    	MPVplotR->SetMaximum(MPVmax);

		if(col){
			MPVplotR->GetXaxis()->SetAxisColor(deepViolet);
			MPVplotR->GetXaxis()->SetLabelColor(deepViolet);
			MPVplotR->GetXaxis()->SetTitleColor(deepViolet);
			MPVplotR->GetYaxis()->SetAxisColor(deepViolet);
			MPVplotR->GetYaxis()->SetLabelColor(deepViolet);
			MPVplotR->GetYaxis()->SetTitleColor(deepViolet);
		}

		expoFitR->SetLineColor(kAzure - 3);
		expoFitR->SetLineWidth(4.0);

		setFontSize<TGraphAsymmErrors>(MPVplotR, 133, 30);

		
		c_split->cd();		
		pad1->Draw();		
		pad1->cd();	
		MPVplotL->Draw("AP");	
		if(plotFit){expoFitL->Draw("SAME");}
		if(plotStats){
			if(plotFit){TPaveText *statsL = statsBox({.45,.65,.80,.88}, track_count, expoFitL, expoFitLResult);statsL->Draw();}
			if(!plotFit){TPaveText *statsL = statsBox({.45,.82,.80,.88}, track_count);statsL->Draw();}
		}
		c_split->cd();
		pad2->Draw();
		pad2->cd();
		MPVplotR->Draw("AP");
		if(plotFit){expoFitR->Draw("SAME");}
		if(plotStats){
		if(plotFit){TPaveText *statsR = statsBox({.45,.65,.80,.88}, track_count, expoFitR, expoFitRResult);statsR->Draw();}
		if(!plotFit){TPaveText *statsR = statsBox({.45,.82,.80,.88}, track_count);statsR->Draw();}
		}
		
        if(col){saveFig(c_split, saveLoc + dataset + "/plot_MPV_colour_tDrift_" + histName + "_" + plotFitLabel + "_" + tag);}
        if(!col){saveFig(c_split, saveLoc + dataset  + "/plot_MPV_tDrift_" + histName + "_" + plotFitLabel + "_" + tag);}

		std::cout << track_count << std::endl;

    }

    return 0;

}