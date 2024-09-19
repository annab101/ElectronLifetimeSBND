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

    Color_t actualWhite = TColor::GetFreeColorIndex();
    TColor *ci = new TColor(actualWhite, 1., 1., 1.);
	
	std::string filename = "noFile";
    std::string histName = "noHist";
    int distOrTime = 0; //0 = distance, 1 = time
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

    TFile f((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

    if(distOrTime == 0){

        TH2D *h;

        h = (TH2D*)f.Get(("h_dQdx_xDrift_" + histName).c_str());

        TCanvas *c_plain = new TCanvas();
        c_plain->SetLeftMargin(0.15);
        c_plain->SetRightMargin(0.18);
        c_plain->SetBottomMargin(0.12);

        h->SetStats(0);
        h->SetTitle(";x_{#scale[1.2]{drift}} (cm);Collection Plane dQ/dx (ADC/cm);Entries");
        h->GetXaxis()->SetLabelOffset(0.01);
        h->GetXaxis()->SetNdivisions(505);
        h->GetYaxis()->SetLabelOffset(0.01);
        h->GetYaxis()->SetNdivisions(505);
        h->GetYaxis()->SetTitleOffset(1.8);
        setFontSize<TH2>(h, 133, 25);
        setFontSizeZ(h, 133, 25);
        h->GetZaxis()->CenterTitle(true);
        h->GetZaxis()->SetTitleOffset(2.0);
        h->Draw("COLZ");

        std::vector<std::string> SBNDlabels = {
            "SBND: Data Run 14608",
            "Anode-cathode crossing tracks",
            "Preliminary"
        };

        TPaveText *labels = plotLabels({.18,.70,.46,.85}, SBNDlabels, actualWhite, 133,22);
        labels->Draw();
        saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_xDrift_" + histName + tag);

    }

    if(distOrTime == 1){

        TH2D *hL, *hR;

        hL = (TH2D*)f.Get(("h_dQdx_tDriftL_" + histName).c_str());
        hR = (TH2D*)f.Get(("h_dQdx_tDriftR_" + histName).c_str());

        TCanvas *c_split = new TCanvas();
        c_split->SetWindowSize(2000,500);
        if(col){
            c_split->SetFillColor(powderBlue);
            c_split->SetFillStyle(1001);
        }
        TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
        TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
        c_split->cd();
        pad1->SetRightMargin(0.22);
        pad1->SetLeftMargin(0.18);
        pad1->SetBottomMargin(0.15);
        pad2->SetRightMargin(0.22);
        pad2->SetLeftMargin(0.18);
        pad2->SetBottomMargin(0.15);
        if(col){
            pad1->SetFillColor(powderBlue);
            pad1->SetFillStyle(4100);
            pad1->SetFrameFillColor(powderBlue);
            pad1->SetFrameFillStyle(4100);
            pad2->SetFillColor(powderBlue);
            pad2->SetFillStyle(4100);
            pad2->SetFrameFillColor(powderBlue);
            pad2->SetFrameFillStyle(4100);
        }

        hL->SetStats(0);
        hL->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (ADC/cm);Entries");
        hL->GetXaxis()->SetLabelOffset(0.01);
        hL->GetXaxis()->SetNdivisions(505);
        hL->GetYaxis()->SetLabelOffset(0.01);
        hL->GetYaxis()->SetNdivisions(505);
        hL->GetYaxis()->SetTitleOffset(1.8);
        if(col){
            hL->GetXaxis()->SetAxisColor(deepViolet);
            hL->GetXaxis()->SetLabelColor(deepViolet);
            hL->GetXaxis()->SetTitleColor(deepViolet);
            hL->GetYaxis()->SetAxisColor(deepViolet);
            hL->GetYaxis()->SetLabelColor(deepViolet);
            hL->GetYaxis()->SetTitleColor(deepViolet);
        }
        setFontSize<TH2>(hL, 133, 30);
        setFontSizeZ(hL, 133, 30);
        hL->GetZaxis()->CenterTitle(true);
        hL->GetZaxis()->SetTitleOffset(2.0);
        if(col){
            hL->GetZaxis()->SetAxisColor(deepViolet);
            hL->GetZaxis()->SetLabelColor(deepViolet);
            hL->GetZaxis()->SetTitleColor(deepViolet);
        }

        hR->SetStats(0);
        hR->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (ADC/cm);Entries");
        hR->GetXaxis()->SetLabelOffset(0.01);
        hR->GetXaxis()->SetNdivisions(505);
        hR->GetYaxis()->SetLabelOffset(0.01);
        hR->GetYaxis()->SetNdivisions(505);
        hR->GetYaxis()->SetTitleOffset(1.8);
        if(col){
            hR->GetXaxis()->SetAxisColor(deepViolet);
            hR->GetXaxis()->SetLabelColor(deepViolet);
            hR->GetXaxis()->SetTitleColor(deepViolet);
            hR->GetYaxis()->SetAxisColor(deepViolet);
            hR->GetYaxis()->SetLabelColor(deepViolet);
            hR->GetYaxis()->SetTitleColor(deepViolet);
        }
        setFontSize<TH2>(hR, 133, 30);
        setFontSizeZ(hR, 133, 30);
        hR->GetZaxis()->CenterTitle(true);
        hR->GetZaxis()->SetTitleOffset(2.0);
        if(col){
            hR->GetZaxis()->SetAxisColor(deepViolet);
            hR->GetZaxis()->SetLabelColor(deepViolet);
            hR->GetZaxis()->SetTitleColor(deepViolet);
        }

        c_split->cd();		
        pad1->Draw();		
        pad1->cd();	
        hL->Draw("COLZ");
        std::vector<std::string> SBNDlabelsE = {
            "SBND: Data Run 14608",
            "Anode-cathode crossing tracks",
            "Preliminary",
            "East TPC"
        };
        TPaveText *labelsE = plotLabels({.21,.65,.49,.85}, SBNDlabelsE, actualWhite, 133,24);
        labelsE->Draw();
        c_split->cd();
        pad2->Draw();
        pad2->cd();
        hR->Draw("COLZ");
        std::vector<std::string> SBNDlabelsW = {
            "SBND: Data Run 14608",
            "Anode-cathode crossing tracks",
            "Preliminary",
            "West TPC"
        };
        TPaveText *labelsW = plotLabels({.21,.65,.49,.85}, SBNDlabelsW, actualWhite, 133,24);
        labelsW->Draw();
        if(col){saveFig(c_split, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_colour_tDrift_" + histName + tag);}
        if(!col){saveFig(c_split, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_tDrift_" + histName + tag);}

    }

    if(distOrTime == 2){ //separate plots for E and W

        TH2D *hL, *hR;

        hL = (TH2D*)f.Get(("h_dQdx_tDriftL_" + histName).c_str());
        hR = (TH2D*)f.Get(("h_dQdx_tDriftR_" + histName).c_str());

        TCanvas *c = new TCanvas();
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.12);

        hL->SetStats(0);
        hL->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (ADC/cm);Entries");
        hL->GetXaxis()->SetLabelOffset(0.01);
        hL->GetXaxis()->SetNdivisions(505);
        hL->GetYaxis()->SetLabelOffset(0.01);
        hL->GetYaxis()->SetNdivisions(505);
        hL->GetYaxis()->SetTitleOffset(1.8);
        setFontSize<TH2>(hL, 133, 25);
        setFontSizeZ(hL, 133, 25);
        hL->GetZaxis()->CenterTitle(true);
        hL->GetZaxis()->SetTitleOffset(2.0);
        
        hR->SetStats(0);
        hR->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (ADC/cm);Entries");
        hR->GetXaxis()->SetLabelOffset(0.01);
        hR->GetXaxis()->SetNdivisions(505);
        hR->GetYaxis()->SetLabelOffset(0.01);
        hR->GetYaxis()->SetNdivisions(505);
        hR->GetYaxis()->SetTitleOffset(1.8);
        setFontSize<TH2>(hR, 133, 25);
        setFontSizeZ(hR, 133, 25);
        hR->GetZaxis()->CenterTitle(true);
        hR->GetZaxis()->SetTitleOffset(2.0);

        c->cd();		
        hL->Draw("COLZ");
        std::vector<std::string> SBNDlabelsE = {
            "SBND: Data Run 14608",
            "Anode-cathode crossing tracks",
            "Preliminary",
            "East TPC"
        };
        TPaveText *labelsE = plotLabels({.18,.65,.46,.85}, SBNDlabelsE, actualWhite, 133,22);
        labelsE->Draw();
        saveFig(c, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_colour_tDriftE_" + histName + tag);

        openAndClear(c);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.12);
        c->cd();
        hR->Draw("COLZ");
        std::vector<std::string> SBNDlabelsW = {
            "SBND: Data Run 14608",
            "Anode-cathode crossing tracks",
            "Preliminary",
            "West TPC"
        };
        TPaveText *labelsW = plotLabels({.18,.65,.46,.85}, SBNDlabelsW, actualWhite, 133,22);
        labelsW->Draw();
        saveFig(c, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_colour_tDriftW_" + histName + tag);

    }

    return 0;

}