//Author: Anna Beever
//Date:   April 2024

//C++ includes

//ROOT includes
#include "Utilities/ROOTincludes.h"

//Local includes
#include "Helpers/Constants.h"
#include "Helpers/PlottingHelpers.h"
#include "Helpers/PlottingHelpers.cpp"
#include "Helpers/FittingHelpers.h"
#include "Helpers/FittingHelpers.cpp"
#include "Utilities/ConfigReader.h"
#include "Utilities/ConfigReader.cpp"

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

	std::string inputData = "noInputData";
	std::string tag = "noTag";
	std::string dataset = "noDataSet";
	std::string saveLoc = "noSaveLoc";
	int codeConfig = 0;
	int maxWireGroup = 20;
	int angBins = 6;
    std::string configLabel = "noConfigLabel";

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);
	p->getValue("maxWireGroup", maxWireGroup);
	p->getValue("angBins", angBins);
    
    configLabel = configMap[codeConfig];

    if(codeConfig == 0 || codeConfig == 3){

        std::cout << "Lifetime v wire plot doesn't work for this configuration!";

    }

    if(codeConfig == 1){

        TFile f1((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

        TFile f11((saveLoc + dataset  + "_" + configLabel + "/multiwire_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
        f11.cd();

        TGraphErrors *g_lifeVwiresXL = new TGraphErrors();
        TGraphErrors *g_lifeVwiresXR = new TGraphErrors();
        TGraphErrors *g_lifeVwiresTL = new TGraphErrors();
        TGraphErrors *g_lifeVwiresTR = new TGraphErrors();

        for(int i = 1; i <= maxWireGroup; i++){

            TF1 *expoFit_xDriftL = (TF1*)f1.Get(TString::Format("MPVfit_xDriftL_%dwires", i));
            TF1 *expoFit_xDriftR = (TF1*)f1.Get(TString::Format("MPVfit_xDriftR_%dwires", i));
            TF1 *expoFit_tDriftL = (TF1*)f1.Get(TString::Format("MPVfit_tDriftL_%dwires", i));
            TF1 *expoFit_tDriftR = (TF1*)f1.Get(TString::Format("MPVfit_tDriftR_%dwires", i));

            double tau_XL = expoFit_xDriftL->GetParameter(1);
            double etau_XL = expoFit_xDriftL->GetParError(1);
            double tau_XR = expoFit_xDriftR->GetParameter(1);
            double etau_XR = expoFit_xDriftR->GetParError(1);
            double tau_TL = expoFit_tDriftL->GetParameter(1);
            double etau_TL = expoFit_tDriftL->GetParError(1);
            double tau_TR = expoFit_tDriftR->GetParameter(1);
            double etau_TR = expoFit_tDriftR->GetParError(1);

            if(tau_XL >= 0. && tau_XL <= 50. && etau_XL <= 50){
				g_lifeVwiresXL->SetPoint(g_lifeVwiresXL->GetN(), i, tau_XL);
				g_lifeVwiresXL->SetPointError(g_lifeVwiresXL->GetN()-1, 0., etau_XL);
			}
			if(tau_XR >= 0. && tau_XR <= 50. && etau_XR <= 50){
				g_lifeVwiresXR->SetPoint(g_lifeVwiresXR->GetN(), i, tau_XR);
				g_lifeVwiresXR->SetPointError(g_lifeVwiresXR->GetN()-1, 0., etau_XR);
			}
			if(tau_TL >= 0. && tau_TL <= 50. && etau_TL <= 50){
				g_lifeVwiresTL->SetPoint(g_lifeVwiresTL->GetN(), i, tau_TL);
				g_lifeVwiresTL->SetPointError(g_lifeVwiresTL->GetN()-1, 0., etau_TL);
			}
			if(tau_TR >= 0. && tau_TR <= 50. && etau_TR <= 50){
				g_lifeVwiresTR->SetPoint(g_lifeVwiresTR->GetN(), i, tau_TR);
				g_lifeVwiresTR->SetPointError(g_lifeVwiresTR->GetN()-1, 0., etau_TR);
			}

        }

        g_lifeVwiresXL->SetTitle("x, left TPC");
        g_lifeVwiresXR->SetTitle("x, right TPC");
        g_lifeVwiresTL->SetTitle("t, left TPC");
        g_lifeVwiresTR->SetTitle("t, right TPC");

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(";N;#tau (ms)");
        mg->Add(g_lifeVwiresXL);
        mg->Add(g_lifeVwiresXR);
        mg->Add(g_lifeVwiresTL);
        mg->Add(g_lifeVwiresTR);

        mg->SetName("lifeVwires");
        
        mg->Write();

    }

    if(codeConfig == 2){

        TFile f2((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

        TFile f22((saveLoc + dataset  + "_" + configLabel + "/multihit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
        f22.cd();

        TGraphErrors *g_lifeVhitsXL = new TGraphErrors();
        TGraphErrors *g_lifeVhitsXR = new TGraphErrors();
        TGraphErrors *g_lifeVhitsTL = new TGraphErrors();
        TGraphErrors *g_lifeVhitsTR = new TGraphErrors();

        for(int i = 1; i <= maxWireGroup; i++){

            TF1 *expoFit_xDriftL = (TF1*)f2.Get(TString::Format("MPVfit_xDriftL_%dhits", i));
            TF1 *expoFit_xDriftR = (TF1*)f2.Get(TString::Format("MPVfit_xDriftR_%dhits", i));
            TF1 *expoFit_tDriftL = (TF1*)f2.Get(TString::Format("MPVfit_tDriftL_%dhits", i));
            TF1 *expoFit_tDriftR = (TF1*)f2.Get(TString::Format("MPVfit_tDriftR_%dhits", i));

            double tau_XL = expoFit_xDriftL->GetParameter(1);
            double etau_XL = expoFit_xDriftL->GetParError(1);
            double tau_XR = expoFit_xDriftR->GetParameter(1);
            double etau_XR = expoFit_xDriftR->GetParError(1);
            double tau_TL = expoFit_tDriftL->GetParameter(1);
            double etau_TL = expoFit_tDriftL->GetParError(1);
            double tau_TR = expoFit_tDriftR->GetParameter(1);
            double etau_TR = expoFit_tDriftR->GetParError(1);

            if(tau_XL >= 0. && tau_XL <= 50. && etau_XL <= 50){
				g_lifeVhitsXL->SetPoint(g_lifeVhitsXL->GetN(), i, tau_XL);
				g_lifeVhitsXL->SetPointError(g_lifeVhitsXL->GetN()-1, 0., etau_XL);
			}
			if(tau_XR >= 0. && tau_XR <= 50. && etau_XR <= 50){
				g_lifeVhitsXR->SetPoint(g_lifeVhitsXR->GetN(), i, tau_XR);
				g_lifeVhitsXR->SetPointError(g_lifeVhitsXR->GetN()-1, 0., etau_XR);
			}
			if(tau_TL >= 0. && tau_TL <= 50. && etau_TL <= 50){
				g_lifeVhitsTL->SetPoint(g_lifeVhitsTL->GetN(), i, tau_TL);
				g_lifeVhitsTL->SetPointError(g_lifeVhitsTL->GetN()-1, 0., etau_TL);
			}
			if(tau_TR >= 0. && tau_TR <= 50. && etau_TR <= 50){
				g_lifeVhitsTR->SetPoint(g_lifeVhitsTR->GetN(), i, tau_TR);
				g_lifeVhitsTR->SetPointError(g_lifeVhitsTR->GetN()-1, 0., etau_TR);
			}

        }

        g_lifeVhitsXL->SetTitle("x, left TPC");
        g_lifeVhitsXR->SetTitle("x, right TPC");
        g_lifeVhitsTL->SetTitle("t, left TPC");
        g_lifeVhitsTR->SetTitle("t, right TPC");

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(";N;#tau (ms)");
        mg->Add(g_lifeVhitsXL);
        mg->Add(g_lifeVhitsXR);
        mg->Add(g_lifeVhitsTL);
        mg->Add(g_lifeVhitsTR);

        mg->SetName("lifeVhits");
        
        mg->Write();

    }

    if(codeConfig == 4){

        TFile f4((saveLoc + dataset  + "_" + configLabel + "/MPVfit_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

        TFile f44((saveLoc + dataset  + "_" + configLabel + "/multiwire_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
        f44.cd();

        TGraphErrors *g_lifeVwires_angXL[2][angBins/2];
        TGraphErrors *g_lifeVwires_angXR[2][angBins/2];
        TGraphErrors *g_lifeVwires_angTL[2][angBins/2];
        TGraphErrors *g_lifeVwires_angTR[2][angBins/2];

        double** angBinLimits = getAngBinLimits(angBins);

         for(int k = 0; k < angBins/2; k++){
            
            for(int j = 0; j <= 1; j++){

                for(int i = 1; i <= maxWireGroup; i++){          

                    TF1 *expoFit_xDriftL = (TF1*)f4.Get(TString::Format("MPVfit_xDriftL_ang%ito%inWires%d", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1], i));
                    TF1 *expoFit_xDriftR = (TF1*)f4.Get(TString::Format("MPVfit_xDriftR_ang%ito%inWires%d", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1], i));
                    TF1 *expoFit_tDriftL = (TF1*)f4.Get(TString::Format("MPVfit_tDriftL_ang%ito%inWires%d", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1], i));
                    TF1 *expoFit_tDriftR = (TF1*)f4.Get(TString::Format("MPVfit_tDriftR_ang%ito%inWires%d", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1], i));

                    double tau_XL = expoFit_xDriftL->GetParameter(1);
                    double etau_XL = expoFit_xDriftL->GetParError(1);
                    double tau_XR = expoFit_xDriftR->GetParameter(1);
                    double etau_XR = expoFit_xDriftR->GetParError(1);
                    double tau_TL = expoFit_tDriftL->GetParameter(1);
                    double etau_TL = expoFit_tDriftL->GetParError(1);
                    double tau_TR = expoFit_tDriftR->GetParameter(1);
                    double etau_TR = expoFit_tDriftR->GetParError(1);

                    if(tau_XL >= 0. && tau_XL <= 50. && etau_XL <= 50){
                        g_lifeVwires_angXL[j][k]->SetPoint(g_lifeVwires_angXL[j][k]->GetN(), i, tau_XL);
                        g_lifeVwires_angXL[j][k]->SetPointError(g_lifeVwires_angXL[j][k]->GetN()-1, 0., etau_XL);
                    }
                    if(tau_XR >= 0. && tau_XR <= 50. && etau_XR <= 50){
                        g_lifeVwires_angXR[j][k]->SetPoint(g_lifeVwires_angXR[j][k]->GetN(), i, tau_XR);
                        g_lifeVwires_angXR[j][k]->SetPointError(g_lifeVwires_angXR[j][k]->GetN()-1, 0., etau_XR);
                    }
                    if(tau_TL >= 0. && tau_TL <= 50. && etau_TL <= 50){
                        g_lifeVwires_angTL[j][k]->SetPoint(g_lifeVwires_angTL[j][k]->GetN(), i, tau_TL);
                        g_lifeVwires_angTL[j][k]->SetPointError(g_lifeVwires_angTL[j][k]->GetN()-1, 0., etau_TL);
                    }
                    if(tau_TR >= 0. && tau_TR <= 50. && etau_TR <= 50){
                        g_lifeVwires_angTR[j][k]->SetPoint(g_lifeVwires_angTR[j][k]->GetN(), i, tau_TR);
                        g_lifeVwires_angTR[j][k]->SetPointError(g_lifeVwires_angTR[j][k]->GetN()-1, 0., etau_TR);
                    }

                }

                g_lifeVwires_angXL[j][k]->SetTitle(TString::Format("x, left TPC, %i#circ to %i#circ", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1]));
                g_lifeVwires_angXR[j][k]->SetTitle(TString::Format("x, right TPC, %i#circ to %i#circ", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1]));
                g_lifeVwires_angTL[j][k]->SetTitle(TString::Format("t, left TPC, %i#circ to %i#circ", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1]));
                g_lifeVwires_angTR[j][k]->SetTitle(TString::Format("t, right TPC, %i#circ to %i#circ", (int)angBinLimits[j][k], (int)angBinLimits[j][k+1]));

            }

        }

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(";N;#tau (ms)");

        for(int k = 0; k < angBins/2; k++){
            
            for(int j = 0; j <= 1; j++){

                mg->Add(g_lifeVwires_angXL[j][k]);
                mg->Add(g_lifeVwires_angXR[j][k]);
                mg->Add(g_lifeVwires_angTL[j][k]);
                mg->Add(g_lifeVwires_angTR[j][k]);

            }
        
        }

        mg->SetName("lifeVhits");
        
        mg->Write();

    }

    return 0;

}