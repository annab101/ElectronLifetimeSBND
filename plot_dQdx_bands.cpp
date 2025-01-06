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

    bool fancyPlot = 1;

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
    std::string MPVfileName = saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root";

    TFile f_bands("/exp/sbnd/data/users/abeever/cosmics_analysis/Lifetime/run14608choppyCut_classic/lifetimestudy.root");

    std::vector<double> lifetimes {3,10};
    std::vector<int> colours {peach, coral};

    if(distOrTime == 0){

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
        h->SetMinimum(0);
        h->Draw("COLZ");

        /*if(fancyPlot){
            plotExp(MPVfileName, lifetimes, colours, "XL", histName);
            plotExp(MPVfileName, lifetimes, colours, "XR", histName);
        }*/

        //Finding the norms for the exponentials for the bands
        TFile f_MPVs((MPVfileName).c_str());
        TGraphAsymmErrors *MPVs_xDrift = (TGraphAsymmErrors*)f_MPVs.Get(("MPVplot_xDrift_" + histName).c_str());
        
        double normL, normR;
        N_MPVs_xDrift = MPVs_xDrift->GetN();

        for(int i=4; i < 8; i++){
            normL+=MPVs_xDrift->GetPointY(i)/4;
            normR+=MPVs_xDrift->GetPointY(N_MPVs_xDrift-i)/4;
        }

        TGraph *sfGraph = (TGraph*)f_bands.Get("lifetimegraph_total_SF");

        int Npoints =  sfGraph->GetN();

        TGraphAsymmErrors *eTimeBand3ms = new TGraphAsymmErrors();
        TGraph *upperLine3ms = new TGraph();
        TGraphAsymmErrors *eTimeBand10ms = new TGraphAsymmErrors();
        TGraph *upperLine10ms = new TGraph();

        std::cout << Npoints << std::endl;

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(i);
            double sf = sfGraph->GetPointY(i);
            double sign = sf-1;
            if(x < 0){
                eTimeBand3ms->SetPoint(i,x,normL*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                upperLine3ms->SetPoint(i,x,sf*normL*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                if(sign>0){
                    eTimeBand3ms->SetPointError(i,0,0,0,(sf-1)*normL*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                }
                else if (sign<0){
                    eTimeBand3ms->SetPointError(i,0,0,(1-sf)*normL*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)),0);
                }
            }
            else if(x > 0){
                eTimeBand3ms->SetPoint(i,x,normR*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                upperLine3ms->SetPoint(i,x,sf*normR*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                if(sign>0){
                    eTimeBand3ms->SetPointError(i,0,0,0,(sf-1)*normR*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)));
                }
                else if (sign<0){
                    eTimeBand3ms->SetPointError(i,0,0,(1-sf)*normR*TMath::Exp(-(200-TMath::Abs(x))/(3*vDrift)),0);
                }
            }
        }

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(i);
            double sf = sfGraph->GetPointY(i);
            double sign = sf-1;
            if(x < 0){
                eTimeBand10ms->SetPoint(i,x,normL*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                upperLine10ms->SetPoint(i,x,sf*normL*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                if(sign>0){
                    eTimeBand10ms->SetPointError(i,0,0,0,(sf-1)*normL*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                }
                else if (sign<0){
                    eTimeBand10ms->SetPointError(i,0,0,(1-sf)*normL*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)),0);
                }
            }
            else if(x > 0){
                eTimeBand10ms->SetPoint(i,x,normR*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                upperLine10ms->SetPoint(i,x,sf*normR*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                if(sign>0){
                    eTimeBand10ms->SetPointError(i,0,0,0,(sf-1)*normR*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)));
                }
                else if (sign<0){
                    eTimeBand10ms->SetPointError(i,0,0,(1-sf)*normR*TMath::Exp(-(200-TMath::Abs(x))/(10*vDrift)),0);
                }
            }
        }

        eTimeBand3ms->SetFillColor(peach);
        eTimeBand3ms->SetLineColor(peach);
        eTimeBand3ms->SetLineWidth(2);
        upperLine3ms->SetLineColor(peach);
        upperLine3ms->SetLineWidth(2);
        eTimeBand3ms->Draw("3L");
        upperLine3ms->Draw("L");

        eTimeBand10ms->SetFillColor(coral);
        eTimeBand10ms->SetLineColor(coral);
        eTimeBand10ms->SetLineWidth(2);
        upperLine10ms->SetLineColor(coral);
        upperLine10ms->SetLineWidth(2);
        eTimeBand10ms->Draw("3L");
        upperLine10ms->Draw("L");


        TBox *anode1 = new TBox(-200, 200, -190, 1800);
        TBox *anode2 = new TBox(190, 200, 200, 1800); 
        TBox *cathode = new TBox(-5, 200, 5, 1800);
        TBox *key = new TBox(120, 1410, 180, 1640);

        anode1->SetFillColor(kGray + 3);
        anode2->SetFillColor(kGray + 3);
        cathode->SetFillColor(kGray + 3);
        key->SetFillColorAlpha(0, 0.3);

        if(fancyPlot){
            anode1->Draw();
            anode2->Draw();
            cathode->Draw();
            key->Draw();
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

        if(fancyPlot){
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

        if(fancyPlot){    
            pt_anode2->Draw();
        }

        TPaveText *pt_cathode = new TPaveText(0.43525, 0.45, 0.53525, 0.6, "blNDC");

        pt_cathode->SetBorderSize(0);
        pt_cathode->SetFillColorAlpha(0,0.0);

        pt_cathode->AddText("CATHODE");
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextColor(actualWhite);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextAlign(22);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextAngle(90);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextFont(133);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextSize(15);
        
        if(fancyPlot){
            pt_cathode->Draw();
        }

        std::vector<std::string> SBNDlabels = {
            "SBND: Data Run 14608",
            "Anode-Cathode-Crossing Tracks",
            "Preliminary"
        };

        TPaveText *labels = plotLabels({.18,.70,.46,.85}, SBNDlabels, actualWhite, 133,20);
        
        if(fancyPlot){
            labels->Draw();
        }

        TPaveText *pt_colour = new TPaveText(.67, .72, .80, .87, "blNDC");

        pt_colour->SetBorderSize(0);
        pt_colour->SetFillColorAlpha(0,0.0);
        pt_colour->SetTextAlign(22);
        pt_colour->SetTextFont(133);
        pt_colour->SetTextSize(20);

        pt_colour->AddText("Lifetime:");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(actualWhite);
        //pt_colour->AddText("1 ms");
        //((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(peach);
        pt_colour->AddText("3 ms");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(peach);
        pt_colour->AddText("10 ms");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(coral);
        //pt_colour->AddText("100 ms");
        //((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(mauve);

        if(fancyPlot){
            pt_colour->Draw();
        }
        saveFig(c_plain, saveLoc + dataset  + "_" + configLabel + "/BANDS_plot_dQdx_xDrift_" + histName + tag);

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
        hL->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (Arb. Units);Entries");
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
        hL->GetZaxis()->SetTitleOffset(1.6);
        if(col){
            hL->GetZaxis()->SetAxisColor(deepViolet);
            hL->GetZaxis()->SetLabelColor(deepViolet);
            hL->GetZaxis()->SetTitleColor(deepViolet);
        }

        hR->SetStats(0);
        hR->SetTitle(";t_{#scale[1.2]{drift}} (ms);Collection Plane dQ/dx (Arb. Units);Entries");
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
        hR->GetZaxis()->SetTitleOffset(1.6);
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
            "Anode-Cathode-Crossing Tracks",
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
            "Anode-Cathode-Crossing Tracks",
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

        TBox *anode = new TBox(0.0, 200, 0.06, 1800);
        TBox *cathode = new TBox(1.24, 200, 1.3, 1800);
        TBox *key = new TBox(1.0, 1410, 1.19, 1640);

        anode->SetFillColor(kGray + 3);
        cathode->SetFillColor(kGray + 3);
        key->SetFillColorAlpha(0, 0.3);

        TPaveText *pt_anode = new TPaveText(0.115, 0.45, 0.215, 0.6, "blNDC");

        pt_anode->SetBorderSize(0);
        pt_anode->SetFillColorAlpha(0,0.0);

        pt_anode->AddText("ANODE");
        ((TText*)pt_anode->GetListOfLines()->Last())->SetTextColor(actualWhite);
        ((TText*)pt_anode->GetListOfLines()->Last())->SetTextAlign(22);
        ((TText*)pt_anode->GetListOfLines()->Last())->SetTextAngle(90);
        ((TText*)pt_anode->GetListOfLines()->Last())->SetTextFont(133);
        ((TText*)pt_anode->GetListOfLines()->Last())->SetTextSize(15);

        TPaveText *pt_cathode = new TPaveText(0.755, 0.45, 0.855, 0.6, "blNDC");

        pt_cathode->SetBorderSize(0);
        pt_cathode->SetFillColorAlpha(0,0.0);

        pt_cathode->AddText("CATHODE");
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextColor(actualWhite);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextAlign(22);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextAngle(90);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextFont(133);
        ((TText*)pt_cathode->GetListOfLines()->Last())->SetTextSize(15);

        TPaveText *pt_colour = new TPaveText(.65, .72, .78, .87, "blNDC");

        pt_colour->SetBorderSize(0);
        pt_colour->SetFillColorAlpha(0,0.0);
        pt_colour->SetTextAlign(22);
        pt_colour->SetTextFont(133);
        pt_colour->SetTextSize(20);

        pt_colour->AddText("Lifetime:");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(actualWhite);
        //pt_colour->AddText("1 ms");
        //((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(peach);
        pt_colour->AddText("3 ms");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(peach);
        pt_colour->AddText("10 ms");
        ((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(coral);
       // pt_colour->AddText("100 ms");
        //((TText*)pt_colour->GetListOfLines()->Last())->SetTextColor(mauve);

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
        hL->GetZaxis()->SetLimits(0,600);
        hL->SetMinimum(0);
        hL->SetMaximum(600);
        
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
        hR->GetZaxis()->SetLimits(0,600);
        hR->SetMinimum(0);
        hR->SetMaximum(600);

        c->cd();		
        hL->Draw("COLZ");

        double normL = 1100;

        TGraph *sfGraph = (TGraph*)f_bands.Get("lifetimegraph_total_SF");

        int Npoints =  sfGraph->GetN();

        TGraphAsymmErrors *eTimeBand3msL = new TGraphAsymmErrors();
        TGraph *upperLine3msL = new TGraph();
        TGraphAsymmErrors *eTimeBand10msL = new TGraphAsymmErrors();
        TGraph *upperLine10msL = new TGraph();

        std::cout << Npoints << std::endl;

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(i);
            double t = (200-TMath::Abs(x))/vDrift;
            double sf = sfGraph->GetPointY(i);
            double sign = sf-1;
            if(x < 0){
                eTimeBand3msL->SetPoint(i,t,normL*TMath::Exp(-t/3));
                upperLine3msL->SetPoint(i,t,sf*normL*TMath::Exp(-t/3));
                if(sign>0){
                    eTimeBand3msL->SetPointError(i,0,0,0,(sf-1)*normL*TMath::Exp(-t/3));
                }
                else if (sign<0){
                    eTimeBand3msL->SetPointError(i,0,0,(1-sf)*normL*TMath::Exp(-t/3),0);
                }
            }
        }

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(i);
            double t = (200-TMath::Abs(x))/vDrift;
            double sf = sfGraph->GetPointY(i);
            double sign = sf-1;
            if(x < 0){
                eTimeBand10msL->SetPoint(i,t,normL*TMath::Exp(-t/10));
                upperLine10msL->SetPoint(i,t,sf*normL*TMath::Exp(-t/10));
                if(sign>0){
                    eTimeBand10msL->SetPointError(i,0,0,0,(sf-1)*normL*TMath::Exp(-t/10));
                }
                else if (sign<0){
                    eTimeBand10msL->SetPointError(i,0,0,(1-sf)*normL*TMath::Exp(-t/10),0);
                }
            }
        }

        eTimeBand3msL->SetFillColor(peach);
        eTimeBand3msL->SetLineColor(peach);
        eTimeBand3msL->SetLineWidth(2);
        upperLine3msL->SetLineColor(peach);
        upperLine3msL->SetLineWidth(2);
        eTimeBand3msL->Draw("3L");
        upperLine3msL->Draw("L");

        eTimeBand10msL->SetFillColor(coral);
        eTimeBand10msL->SetLineColor(coral);
        eTimeBand10msL->SetLineWidth(2);
        upperLine10msL->SetLineColor(coral);
        upperLine10msL->SetLineWidth(2);
        eTimeBand10msL->Draw("3L");
        upperLine10msL->Draw("L");

        /*if(fancyPlot){
            plotExp(MPVfileName, lifetimes, colours, "TL", histName);
        }*/

        if(fancyPlot){
            anode->Draw();
            cathode->Draw();
            key->Draw();
            gPad->RedrawAxis();

            pt_anode->Draw();
            pt_cathode->Draw();
        }

        std::vector<std::string> SBNDlabelsE = {
            "SBND: Data Run 14608",
            "Anode-Cathode-Crossing Tracks",
            "Preliminary",
            "East TPC"
        };
        TPaveText *labelsE = plotLabels({.18,.65,.46,.85}, SBNDlabelsE, actualWhite, 133,20);
        if(fancyPlot){
            //labelsE->Draw();
            //pt_colour->Draw();
        }

        saveFig(c, saveLoc + dataset  + "_" + configLabel + "/BANDS_FINAL_plot_dQdx_tDriftE_" + histName + tag);

        openAndClear(c);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.12);
        c->cd();
        hR->Draw("COLZ");

        double normR = 1100;

        TGraphAsymmErrors *eTimeBand3msR = new TGraphAsymmErrors();
        TGraph *upperLine3msR = new TGraph();
        TGraphAsymmErrors *eTimeBand10msR = new TGraphAsymmErrors();
        TGraph *upperLine10msR = new TGraph();

        std::cout << Npoints << std::endl;

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(Npoints-i-1);
            double t = (200-TMath::Abs(x))/vDrift;
            double sf = sfGraph->GetPointY(Npoints-i-1);
            double sign = sf-1;
            if(x > 0){
                std::cout << t << std::endl;
                std::cout << normR*TMath::Exp(-t/3) << std::endl;
                std::cout << sf*normR*TMath::Exp(-t/3) << std::endl;

                eTimeBand3msR->SetPoint(i,t,normR*TMath::Exp(-t/3));
                upperLine3msR->SetPoint(i,t,sf*normR*TMath::Exp(-t/3));
                if(sign>0){
                    eTimeBand3msR->SetPointError(i,0,0,0,(sf-1)*normR*TMath::Exp(-t/3));
                }
                else if (sign<0){
                    eTimeBand3msR->SetPointError(i,0,0,(1-sf)*normR*TMath::Exp(-t/3),0);
                }
            }
        }

        for(int i = 0; i < Npoints; i++){
            double x = sfGraph->GetPointX(Npoints-i-1);
            double t = (200-TMath::Abs(x))/vDrift;
            double sf = sfGraph->GetPointY(Npoints-i-1);
            double sign = sf-1;
            if(x > 0){
                eTimeBand10msR->SetPoint(i,t,normR*TMath::Exp(-t/10));
                upperLine10msR->SetPoint(i,t,sf*normR*TMath::Exp(-t/10));
                if(sign>0){
                    eTimeBand10msR->SetPointError(i,0,0,0,(sf-1)*normR*TMath::Exp(-t/10));
                }
                else if (sign<0){
                    eTimeBand10msR->SetPointError(i,0,0,(1-sf)*normR*TMath::Exp(-t/10),0);
                }
            }
        }

        eTimeBand3msR->SetFillColor(peach);
        eTimeBand3msR->SetLineColor(peach);
        eTimeBand3msR->SetLineWidth(2);
        upperLine3msR->SetLineColor(peach);
        upperLine3msR->SetLineWidth(2);
        eTimeBand3msR->Draw("3L");
        upperLine3msR->Draw("L");

        eTimeBand10msR->SetFillColor(coral);
        eTimeBand10msR->SetLineColor(coral);
        eTimeBand10msR->SetLineWidth(2);
        upperLine10msR->SetLineColor(coral);
        upperLine10msR->SetLineWidth(2);
        eTimeBand10msR->Draw("3L");
        upperLine10msR->Draw("L");

        /*if(fancyPlot){
            plotExp(MPVfileName, lifetimes, colours, "TR", histName);
        }*/

        if(fancyPlot){
            anode->Draw();
            cathode->Draw();
            key->Draw();
            gPad->RedrawAxis();

            pt_anode->Draw();
            pt_cathode->Draw();
        }

        std::vector<std::string> SBNDlabelsW = {
            "SBND: Data Run 14608",
            "Anode-Cathode-Crossing Tracks",
            "Preliminary",
            "West TPC"
        };
        TPaveText *labelsW = plotLabels({.18,.65,.46,.85}, SBNDlabelsW, actualWhite, 133,20);
        if(fancyPlot){
            //labelsW->Draw();
            //pt_colour->Draw();
        }
        
        saveFig(c, saveLoc + dataset  + "_" + configLabel + "/BANDS_FINAL_plot_dQdx_tDriftW_" + histName + tag);

    }

    return 0;

}