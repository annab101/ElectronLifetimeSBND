#include "PlottingHelpers.h"

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

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TF1* f2){
	
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f1->GetChisquare(),f1->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("#tau: %g#pm%g",f1->GetParameter(1),f1->GetParError(1)));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f2->GetChisquare(),f2->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
        pt->AddText(Form("#tau: %g#pm%g",f2->GetParameter(1),f2->GetParError(1)));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
        
        return pt;
    }

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1){
        
        TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->AddText(Form("AC cosmics: %d",track_count));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
        pt->AddText(Form("#chi^{2} / DoF: %g / %i",f1->GetChisquare(),f1->GetNDF()));
        ((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
        pt->AddText(Form("Lifetime: %g#pm%g",f1->GetParameter(1),f1->GetParError(1)));
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
}