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

    bool hideWirePlanes = 0;

    Color_t actualWhite = TColor::GetFreeColorIndex();
    TColor *ci = new TColor(actualWhite, 1., 1., 1.);
	
	std::string filename = "noFile";
    std::string histName = "noHist";

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
        if(!strcmp(argv[i], "--histName")){
			histName = argv[i+1];
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

    /*==============================================================
                    HISTOGRAMS FOR DRIFT DISTANCE
    ==============================================================*/
    
    TH2D *h;
    h = (TH2D*)f.Get(("h_dQdx_xDrift_" + histName).c_str());

    TCanvas *c_plain = new TCanvas();
    c_plain->SetLeftMargin(0.15);
    c_plain->SetRightMargin(0.18);
    c_plain->SetBottomMargin(0.13);

    h->SetStats(0);
    h->SetTitle(";x_{#scale[1.2]{drift}} (cm);Collection Plane dQ/dx (Arb. Units);Entries");
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetNdivisions(505);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitleOffset(1.8);
    setFontSize<TH2>(h, 133, 25);
    setFontSizeZ(h, 133, 25);
    h->GetZaxis()->CenterTitle(true);
    h->GetZaxis()->SetTitleOffset(1.6);
    h->Draw("COLZ");

    TBox *anode1 = new TBox(-200, 200, -190, 1800);
    TBox *anode2 = new TBox(190, 200, 200, 1800); 
    TBox *cathode_xDrift = new TBox(-5, 200, 5, 1800);

    anode1->SetFillColor(kGray + 3);
    anode2->SetFillColor(kGray + 3);
    cathode_xDrift->SetFillColor(kGray + 3);

    if(!hideWirePlanes){
        anode1->Draw();
        anode2->Draw();
        cathode_xDrift->Draw();
        gPad->RedrawAxis();
    }

    TPaveText *pt_anode1 = new TPaveText(0.1075, 0.45, 0.2075, 0.6, "blNDC");

    pt_anode1->SetBorderSize(0);
    pt_anode1->SetFillColorAlpha(0,0.0);

    pt_anode1->AddText("ANODE");
    ((TText*)pt_anode1->GetListOfLines()->Last())->SetTextColor(actualWhite);
    ((TText*)pt_anode1->GetListOfLines()->Last())->SetTextAlign(22);
    ((TText*)pt_anode1->GetListOfLines()->Last())->SetTextAngle(90);
    ((TText*)pt_anode1->GetListOfLines()->Last())->SetTextFont(133);
    ((TText*)pt_anode1->GetListOfLines()->Last())->SetTextSize(15);

    if(!hideWirePlanes){
        pt_anode1->Draw();
    }

    TPaveText *pt_anode2 = new TPaveText(0.763, 0.45, 0.863, 0.6, "blNDC");

    pt_anode2->SetBorderSize(0);
    pt_anode2->SetFillColorAlpha(0,0.0);

    pt_anode2->AddText("ANODE");
    ((TText*)pt_anode2->GetListOfLines()->Last())->SetTextColor(actualWhite);
    ((TText*)pt_anode2->GetListOfLines()->Last())->SetTextAlign(22);
    ((TText*)pt_anode2->GetListOfLines()->Last())->SetTextAngle(90);
    ((TText*)pt_anode2->GetListOfLines()->Last())->SetTextFont(133);
    ((TText*)pt_anode2->GetListOfLines()->Last())->SetTextSize(15);

    if(!hideWirePlanes){    
        pt_anode2->Draw();
    }

    TPaveText *pt_cathode_xDrift = new TPaveText(0.43525, 0.45, 0.53525, 0.6, "blNDC");

    pt_cathode_xDrift->SetBorderSize(0);
    pt_cathode_xDrift->SetFillColorAlpha(0,0.0);

    pt_cathode_xDrift->AddText("CATHODE");
    ((TText*)pt_cathode_xDrift->GetListOfLines()->Last())->SetTextColor(actualWhite);
    ((TText*)pt_cathode_xDrift->GetListOfLines()->Last())->SetTextAlign(22);
    ((TText*)pt_cathode_xDrift->GetListOfLines()->Last())->SetTextAngle(90);
    ((TText*)pt_cathode_xDrift->GetListOfLines()->Last())->SetTextFont(133);
    ((TText*)pt_cathode_xDrift->GetListOfLines()->Last())->SetTextSize(15);
    
    if(!hideWirePlanes){
        pt_cathode_xDrift->Draw();
    }

    std::vector<std::string> SBNDlabels = {
        "SBND " + dataset,
        "Anode-Cathode-Crossing Tracks"
    };

    TPaveText *labels = plotLabels({.18,.75,.46,.85}, SBNDlabels, actualWhite, 133,20);
    
    labels->Draw();
    
    saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_xDrift_" + dataset + " " + histName + " " + tag);

    /*==============================================================
                        HISTOGRAMS FOR DRIFT TIME
    ==============================================================*/

    TH2D *hL, *hR;

    hL = (TH2D*)f.Get(("h_dQdx_tDriftL_" + histName).c_str());
    hR = (TH2D*)f.Get(("h_dQdx_tDriftR_" + histName).c_str());

    TBox *anode = new TBox(0.0, 200, 0.06, 1800);
    TBox *cathode_tDrift = new TBox(1.24, 200, 1.3, 1800);

    anode->SetFillColor(kGray + 3);
    cathode_tDrift->SetFillColor(kGray + 3);

    TPaveText *pt_anode = new TPaveText(0.115, 0.45, 0.215, 0.6, "blNDC");

    pt_anode->SetBorderSize(0);
    pt_anode->SetFillColorAlpha(0,0.0);

    pt_anode->AddText("ANODE");
    ((TText*)pt_anode->GetListOfLines()->Last())->SetTextColor(actualWhite);
    ((TText*)pt_anode->GetListOfLines()->Last())->SetTextAlign(22);
    ((TText*)pt_anode->GetListOfLines()->Last())->SetTextAngle(90);
    ((TText*)pt_anode->GetListOfLines()->Last())->SetTextFont(133);
    ((TText*)pt_anode->GetListOfLines()->Last())->SetTextSize(15);

    TPaveText *pt_cathode_tDrift = new TPaveText(0.755, 0.45, 0.855, 0.6, "blNDC");

    pt_cathode_tDrift->SetBorderSize(0);
    pt_cathode_tDrift->SetFillColorAlpha(0,0.0);

    pt_cathode_tDrift->AddText("CATHODE");
    ((TText*)pt_cathode_tDrift->GetListOfLines()->Last())->SetTextColor(actualWhite);
    ((TText*)pt_cathode_tDrift->GetListOfLines()->Last())->SetTextAlign(22);
    ((TText*)pt_cathode_tDrift->GetListOfLines()->Last())->SetTextAngle(90);
    ((TText*)pt_cathode_tDrift->GetListOfLines()->Last())->SetTextFont(133);
    ((TText*)pt_cathode_tDrift->GetListOfLines()->Last())->SetTextSize(15);

    TCanvas *c = new TCanvas();
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.18);
    c->SetBottomMargin(0.12);

    hL->SetStats(0);
    hL->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (Arb. Units);Entries");
    hL->GetXaxis()->SetLabelOffset(0.01);
    hL->GetXaxis()->SetNdivisions(505);
    hL->GetYaxis()->SetLabelOffset(0.01);
    hL->GetYaxis()->SetNdivisions(505);
    hL->GetYaxis()->SetTitleOffset(1.8);
    setFontSize<TH2>(hL, 133, 25);
    setFontSizeZ(hL, 133, 25);
    hL->GetZaxis()->CenterTitle(true);
    hL->GetZaxis()->SetTitleOffset(1.6);
    
    hR->SetStats(0);
    hR->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (Arb. Units);Entries");
    hR->GetXaxis()->SetLabelOffset(0.01);
    hR->GetXaxis()->SetNdivisions(505);
    hR->GetYaxis()->SetLabelOffset(0.01);
    hR->GetYaxis()->SetNdivisions(505);
    hR->GetYaxis()->SetTitleOffset(1.8);
    setFontSize<TH2>(hR, 133, 25);
    setFontSizeZ(hR, 133, 25);
    hR->GetZaxis()->CenterTitle(true);
    hR->GetZaxis()->SetTitleOffset(1.6);

    c->cd();		
    hL->Draw("COLZ");

    if(!hideWirePlanes){
        anode->Draw();
        cathode_tDrift->Draw();
        gPad->RedrawAxis();

        pt_anode->Draw();
        pt_cathode_tDrift->Draw();
    }

    std::vector<std::string> SBNDlabelsE = {
        "SBND " + dataset,
        "Anode-Cathode-Crossing Tracks",
        "East TPC"
    };

    TPaveText *labelsE = plotLabels({.18,.70,.46,.85}, SBNDlabelsE, actualWhite, 133,20);
    labelsE->Draw();

    saveFig(c, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_tDriftE_" + dataset + " " + histName + " " + tag);

    openAndClear(c);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.18);
    c->SetBottomMargin(0.12);
    c->cd();
    hR->Draw("COLZ");

    if(!hideWirePlanes){
        anode->Draw();
        cathode_tDrift->Draw();
        gPad->RedrawAxis();

        pt_anode->Draw();
        pt_cathode_tDrift->Draw();
    }

    std::vector<std::string> SBNDlabelsW = {
        "SBND " + dataset,
        "Anode-Cathode-Crossing Tracks",
        "West TPC"
    };

    TPaveText *labelsW = plotLabels({.18,.70,.46,.85}, SBNDlabelsW, actualWhite, 133,20);
    labelsW->Draw();
    
    saveFig(c, saveLoc + dataset  + "_" + configLabel + "/plot_dQdx_tDriftW_" + dataset + " " + histName + " " + tag);
    

    return 0;

}