#include "PlottingHelpers.h"

using namespace constants;

namespace calib {
    
    template<class T> void setMarker(T *t, Color_t color, Style_t style, Size_t size){
        
        t->SetMarkerColor(color);
        t->SetMarkerStyle(style);
        t->SetMarkerSize(size);
    }

    void saveFig(TCanvas *c, std::string plotName){

        std::string saveLocTempPNG = plotName + ".png";
        std::string saveLocTempPDF = plotName + ".pdf";
        std::string saveLocTempROOT = plotName + ".root";
        
        c->SaveAs(saveLocTempPNG.c_str());
        c->SaveAs(saveLocTempPDF.c_str());
        c->SaveAs(saveLocTempROOT.c_str());	

    }

    void openAndClear(TCanvas *c){
	
        c->cd();
        c->Clear();
        
    }

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TF1* f2, TFitResult* r1, TFitResult* r2){
	
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        //calculate lifetime from 1/lifetime

        double tau1 = 1/f1->GetParameter(1);
        double tau2 = 1/f2->GetParameter(1);

        double elow1 = -((1/f1->GetParameter(1)) - (1/(f1->GetParameter(1) + r1->UpperError(1))));
        double ehigh1 = ((1/(f1->GetParameter(1)+r1->LowerError(1))) - (1/f1->GetParameter(1)));

        double elow2 = -((1/f2->GetParameter(1)) - (1/(f2->GetParameter(1) + r2->UpperError(1))));
        double ehigh2 = ((1/(f2->GetParameter(1)+r2->LowerError(1))) - (1/f2->GetParameter(1)));

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f1->GetChisquare(),f1->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("#tau: %.5f^{+%.5f}_{%.5f}",tau1,ehigh1,elow1));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f2->GetChisquare(),f2->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
        pt->AddText(Form("#tau: %.5f^{+%.5f}_{%.5f}",tau2,ehigh2,elow2));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
        
        return pt;
    }

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TFitResult* r1){
        
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        //calculate lifetime from 1/lifetime

        double tau1 = 1/f1->GetParameter(1);

        double elow1 = -((1/f1->GetParameter(1)) - (1/(f1->GetParameter(1) + r1->UpperError(1))));
        double ehigh1 = ((1/(f1->GetParameter(1)+r1->LowerError(1))) - (1/f1->GetParameter(1)));

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f1->GetChisquare(),f1->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("#tau: %.5f^{+%.5f}_{%.5f}",tau1,ehigh1,elow1));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        
        return pt;
    }

    TPaveText *statsBox(std::vector<double> pos, int track_count){
        
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        
        return pt;
    }

    TPaveText *plotLabels(std::vector<double> pos, std::vector<std::string> lines, Color_t textColor, int font, int size){
        
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        pt->SetBorderSize(0);
        pt->SetFillColorAlpha(0,0.0);
        pt->SetTextAlign(11);
        pt->SetTextFont(font);
        pt->SetTextSize(size);

        for(int i=0; i < lines.size(); i++){
            pt->AddText(lines[i].c_str());
            ((TText*)pt->GetListOfLines()->Last())->SetTextColor(textColor);
        }
        
        return pt;
    }

    template<class T> void setFontSize(T *t, int font, int size){
	
        t->GetXaxis()->SetTitleFont(font);
        t->GetXaxis()->SetTitleSize(size);
        t->GetXaxis()->SetLabelFont(font);
        t->GetXaxis()->SetLabelSize(size);

        t->GetYaxis()->SetTitleFont(font);
        t->GetYaxis()->SetTitleSize(size);
        t->GetYaxis()->SetLabelFont(font);
        t->GetYaxis()->SetLabelSize(size);

    }

    void setFontSizeZ(TH2 *h, int font, int size){
        
        h->GetZaxis()->SetTitleFont(font);
        h->GetZaxis()->SetTitleSize(size);
        h->GetZaxis()->SetLabelFont(font);
        h->GetZaxis()->SetLabelSize(size);

    }

    void plotExp(std::string MPVfileName, std::vector<double> lifetimes, std::vector<int> colours, std::string plotSetup, std::string histName){
        
        TFile f((MPVfileName).c_str());
        std::pair<std::string, std::pair<double,double>> thisSetup = whichSetup[plotSetup];
        TF1 *expoFit = (TF1*)f.Get(("MPVfit_" + thisSetup.first + "_" + histName).c_str());
        
        if(strcmp(plotSetup.c_str(), "XL") == 0 || strcmp(plotSetup.c_str(), "XR") == 0){

            for(int i=0; i < lifetimes.size(); i++){

                auto expoFuncX = new TF1("f1",("[0]*exp(-((200-TMath::Abs(x))/([1]*" + std::to_string(vDrift) + ")))").c_str(), thisSetup.second.first, thisSetup.second.second);
                expoFuncX->SetParameter(0, expoFit->GetParameter(0));
                expoFuncX->SetParameter(1,lifetimes[i]);
                expoFuncX->SetLineColor(colours[i]);
                expoFuncX->SetLineWidth(3);
                gStyle->SetLineStyleString(11,"20 20");
                expoFuncX->SetLineStyle(11);
                expoFuncX->Draw("same");

            }

        }
        if(strcmp(plotSetup.c_str(), "TL") == 0 || strcmp(plotSetup.c_str(), "TR") == 0){

            for(int i=0; i < lifetimes.size(); i++){

                auto expoFuncT = new TF1("f1","[0]*exp(-(x/[1]))", thisSetup.second.first, thisSetup.second.second);
                expoFuncT->SetParameter(0, expoFit->GetParameter(0));
                expoFuncT->SetParameter(1,lifetimes[i]);
                expoFuncT->SetLineColor(colours[i]);
                expoFuncT->SetLineWidth(3);
                gStyle->SetLineStyleString(11,"20 20");
                expoFuncT->SetLineStyle(11);
                expoFuncT->Draw("same");

            }

        }
        

    }
}