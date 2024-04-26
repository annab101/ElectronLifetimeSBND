//Author: Anna Beever
//Date:   April 2024

//C++ includes
#include<fstream>
#include<filesystem>
#include<valarray>

//ROOT includes
#include "Utilities/ROOTincludes.h"

//Local includes
#include "Helpers/Constants.h"
#include "Helpers/PlottingHelpers.h"
#include "Helpers/PlottingHelpers.cpp"
#include "Helpers/FittingHelpers.h"
#include "Helpers/FittingHelpers.cpp"
#include "Helpers/StatsHelpers.h"
#include "Helpers/StatsHelpers.cpp"
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
	
	//codeConfig (0 = lifetime, 1 = lifetime with wires grouped, 2 = lifetime with hits grouped, 3 = angular bins, 4 = wires and angular bins)
	int codeConfig = 0;
	int maxWireGroup = 20;
	int NBinsX = 100;
	int NBinsT = 100;
	int NBinsdQdx = 75;
	double minX = -200.;
	double maxX = 200.;
	double mindQdx = 200.;
	double maxdQdx = 1800.;
	double minT = 0.;
	double maxT = 1.3;	
	int angBins = 6;
	int writeProjY = 0;

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("codeConfig", codeConfig);
	p->getValue("maxWireGroup", maxWireGroup);
	p->getValue("NBinsX", NBinsX);
	p->getValue("NBinsT", NBinsT);
	p->getValue("NBinsdQdx", NBinsdQdx);
	p->getValue("minX", minX);
	p->getValue("maxX", maxX);
	p->getValue("mindQdx", mindQdx);
	p->getValue("maxdQdx", maxdQdx);
	p->getValue("minT", minT);
	p->getValue("maxT", maxT);
	p->getValue("angBins", angBins);
	p->getValue("writeProjY", writeProjY);


	if(codeConfig == 0){

		std::string configLabel = "classic";

		TFile f0((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

		TH2D* h_dQdx_xDrift_basic = (TH2D*)f0.Get("h_dQdx_xDrift_basic");
		TH2D* h_dQdx_tDriftL_basic = (TH2D*)f0.Get("h_dQdx_tDriftL_basic");
		TH2D* h_dQdx_tDriftR_basic = (TH2D*)f0.Get("h_dQdx_tDriftR_basic");
		
		TH1D* projY_xDrift_basic = new TH1D();
		TH1D* projY_tDriftL_basic = new TH1D();
		TH1D* projY_tDriftR_basic = new TH1D();

		TF1* LGfit_xDrift_basic = new TF1();
		TF1* LGfit_tDriftL_basic = new TF1();
		TF1* LGfit_tDriftR_basic = new TF1();

		fitLGParameters fp_xDrift;
		fitLGParameters fp_tDriftL;
		fitLGParameters fp_tDriftR;

		TGraphErrors *MPVplot_xDrift_basic = new TGraphErrors();
		TGraphErrors *MPVplot_tDriftL_basic = new TGraphErrors();
		TGraphErrors *MPVplot_tDriftR_basic = new TGraphErrors();

		TFile f00((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
		f00.cd();

		for(int i = 1; i <= NBinsX; i++){

			std::cout << "--------- Fitting bin " + to_string(i) + "----------" << std::endl;

			projY_xDrift_basic = h_dQdx_xDrift_basic->ProjectionY(TString::Format("projY_xDrift_basic_bin%d", i), i, i);
			projY_tDriftL_basic = h_dQdx_tDriftL_basic->ProjectionY(TString::Format("projY_tDriftL_basic_bin%d", i), i, i);
			projY_tDriftR_basic = h_dQdx_tDriftR_basic->ProjectionY(TString::Format("projY_tDriftR_basic_bin%d", i), i, i);

			if((bool)writeProjY){

				projY_xDrift_basic->Write();
				projY_tDriftL_basic->Write();
				projY_tDriftR_basic->Write();

			}

			SetLGParameters(projY_xDrift_basic, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
			SetLGParameters(projY_tDriftL_basic, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.lb, fp_tDriftL.ub);
			SetLGParameters(projY_tDriftR_basic, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.lb, fp_tDriftR.ub);

			LGfit_xDrift_basic = fitter(projY_xDrift_basic, fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
			LGfit_tDriftL_basic = fitter(projY_tDriftL_basic, fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.cov, "LG");
			LGfit_tDriftR_basic = fitter(projY_tDriftR_basic, fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.cov, "LG");

			if(fp_xDrift.fp[1] >= 0.){

				MPVplot_xDrift_basic->SetPoint(MPVplot_xDrift_basic->GetN(), h_dQdx_xDrift_basic->GetXaxis()->GetBinCenter(i), fp_xDrift.fp[1]);
				MPVplot_xDrift_basic->SetPointError(MPVplot_xDrift_basic->GetN() - 1, 0., pointError(projY_xDrift_basic, fp_xDrift));
			}

			if(fp_tDriftL.fp[1] >= 0.){

				MPVplot_tDriftL_basic->SetPoint(MPVplot_tDriftL_basic->GetN(), h_dQdx_tDriftL_basic->GetXaxis()->GetBinCenter(i), fp_tDriftL.fp[1]);
				MPVplot_tDriftL_basic->SetPointError(MPVplot_tDriftL_basic->GetN() - 1, 0., pointError(projY_tDriftL_basic, fp_tDriftL));
			}

			if(fp_tDriftR.fp[1] >= 0.){

				MPVplot_tDriftR_basic->SetPoint(MPVplot_tDriftR_basic->GetN(), h_dQdx_tDriftR_basic->GetXaxis()->GetBinCenter(i), fp_tDriftR.fp[1]);
				MPVplot_tDriftR_basic->SetPointError(MPVplot_tDriftR_basic->GetN() - 1, 0., pointError(projY_tDriftR_basic, fp_tDriftR));
			}

		}

		MPVplot_xDrift_basic->SetName("MPVplot_xDrift_basic");
		MPVplot_tDriftL_basic->SetName("MPVplot_tDriftL_basic");
		MPVplot_tDriftR_basic->SetName("MPVplot_tDriftR_basic");

		MPVplot_xDrift_basic->Write();
		MPVplot_tDriftL_basic->Write();
		MPVplot_tDriftR_basic->Write();

	}

	if(codeConfig == 1){

		std::string configLabel = "multiwire";

		TFile f1((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

		TH2D* h_dQdx_xDrift_wires[maxWireGroup];
		TH2D* h_dQdx_tDriftL_wires[maxWireGroup];
		TH2D* h_dQdx_tDriftR_wires[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){

			h_dQdx_xDrift_wires[i-1] = (TH2D*)f1.Get(TString::Format("h_dQdx_xDrift_%dwires", i));
			h_dQdx_tDriftL_wires[i-1] = (TH2D*)f1.Get(TString::Format("h_dQdx_tDriftL_%dwires", i));
			h_dQdx_tDriftR_wires[i-1] = (TH2D*)f1.Get(TString::Format("h_dQdx_tDriftR_%dwires", i));

		}

		TH1D* projY_xDrift_wires[maxWireGroup];
		TH1D* projY_tDriftL_wires[maxWireGroup];
		TH1D* projY_tDriftR_wires[maxWireGroup];

		TF1* LGfit_xDrift_wires[maxWireGroup];
		TF1* LGfit_tDriftL_wires[maxWireGroup];
		TF1* LGfit_tDriftR_wires[maxWireGroup];

		fitLGParameters fp_xDrift;
		fitLGParameters fp_tDriftL;
		fitLGParameters fp_tDriftR;

		TGraphErrors *MPVplot_xDrift_wires[maxWireGroup];
		TGraphErrors *MPVplot_tDriftL_wires[maxWireGroup];
		TGraphErrors *MPVplot_tDriftR_wires[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){

			MPVplot_xDrift_wires[i-1] = new TGraphErrors(); 
			MPVplot_tDriftL_wires[i-1] = new TGraphErrors();
			MPVplot_tDriftR_wires[i-1] = new TGraphErrors();

		}

		TFile f11((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
		f11.cd();

		for(int i = 1; i <= maxWireGroup; i++){

			for(int j = 1; j <= NBinsX; j++){

				std::cout << "--------- Fitting wire group " + to_string(i) + "; x bin " + to_string(j) + "----------" << std::endl;

				fp_xDrift.reset();
				fp_tDriftL.reset();
				fp_tDriftR.reset();

				projY_xDrift_wires[i-1] = h_dQdx_xDrift_wires[i-1]->ProjectionY(TString::Format("projY_xDrift_%dwires_bin%d", i, j), j, j);
				projY_tDriftL_wires[i-1] = h_dQdx_tDriftL_wires[i-1]->ProjectionY(TString::Format("projY_tDriftL_%dwires_bin%d", i, j), j, j);
				projY_tDriftR_wires[i-1] = h_dQdx_tDriftR_wires[i-1]->ProjectionY(TString::Format("projY_tDriftR_%dwires_bin%d", i, j), j, j);
				
				if((bool)writeProjY){

					projY_xDrift_wires[i-1]->Write();
					projY_tDriftL_wires[i-1]->Write();
					projY_tDriftR_wires[i-1]->Write();

				}

				SetLGParameters(projY_xDrift_wires[i-1], fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
				SetLGParameters(projY_tDriftL_wires[i-1], fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.lb, fp_tDriftL.ub);
				SetLGParameters(projY_tDriftR_wires[i-1], fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.lb, fp_tDriftR.ub);

				LGfit_xDrift_wires[i-1] = fitter(projY_xDrift_wires[i-1], fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
				LGfit_tDriftL_wires[i-1] = fitter(projY_tDriftL_wires[i-1], fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.cov, "LG");
				LGfit_tDriftR_wires[i-1] = fitter(projY_tDriftR_wires[i-1], fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.cov, "LG");

				if(fp_xDrift.fp[1] >= 0.){

					MPVplot_xDrift_wires[i-1]->SetPoint(MPVplot_xDrift_wires[i-1]->GetN(), h_dQdx_xDrift_wires[i-1]->GetXaxis()->GetBinCenter(j), fp_xDrift.fp[1]);
					MPVplot_xDrift_wires[i-1]->SetPointError(MPVplot_xDrift_wires[i-1]->GetN() - 1, 0., pointError(projY_xDrift_wires[i-1], fp_xDrift));
				}

				if(fp_tDriftL.fp[1] >= 0.){

					MPVplot_tDriftL_wires[i-1]->SetPoint(MPVplot_tDriftL_wires[i-1]->GetN(), h_dQdx_tDriftL_wires[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftL.fp[1]);
					MPVplot_tDriftL_wires[i-1]->SetPointError(MPVplot_tDriftL_wires[i-1]->GetN() - 1, 0., pointError(projY_tDriftL_wires[i-1], fp_tDriftL));
				}

				if(fp_tDriftR.fp[1] >= 0.){

					MPVplot_tDriftR_wires[i-1]->SetPoint(MPVplot_tDriftR_wires[i-1]->GetN(), h_dQdx_tDriftR_wires[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftR.fp[1]);
					MPVplot_tDriftR_wires[i-1]->SetPointError(MPVplot_tDriftR_wires[i-1]->GetN() - 1, 0., pointError(projY_tDriftR_wires[i-1], fp_tDriftR));
				}

			}

			MPVplot_xDrift_wires[i-1]->SetName(TString::Format("MPVplot_xDrift_%dwires", i));
			MPVplot_tDriftL_wires[i-1]->SetName(TString::Format("MPVplot_tDriftL_%dwires", i));
			MPVplot_tDriftR_wires[i-1]->SetName(TString::Format("MPVplot_tDriftR_%dwires", i));

			MPVplot_xDrift_wires[i-1]->Write();
			MPVplot_tDriftL_wires[i-1]->Write();
			MPVplot_tDriftR_wires[i-1]->Write();

		}

	}

	if(codeConfig == 2){

		std::string configLabel = "multihit";

		TFile f2((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

		TH2D* h_dQdx_xDrift_hits[maxWireGroup];
		TH2D* h_dQdx_tDriftL_hits[maxWireGroup];
		TH2D* h_dQdx_tDriftR_hits[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){

			h_dQdx_xDrift_hits[i-1] = (TH2D*)f2.Get(TString::Format("h_dQdx_xDrift_%dhits", i));
			h_dQdx_tDriftL_hits[i-1] = (TH2D*)f2.Get(TString::Format("h_dQdx_tDriftL_%dhits", i));
			h_dQdx_tDriftR_hits[i-1] = (TH2D*)f2.Get(TString::Format("h_dQdx_tDriftR_%dhits", i));

		}

		TH1D* projY_xDrift_hits[maxWireGroup];
		TH1D* projY_tDriftL_hits[maxWireGroup];
		TH1D* projY_tDriftR_hits[maxWireGroup];

		TF1* LGfit_xDrift_hits[maxWireGroup];
		TF1* LGfit_tDriftL_hits[maxWireGroup];
		TF1* LGfit_tDriftR_hits[maxWireGroup];

		fitLGParameters fp_xDrift;
		fitLGParameters fp_tDriftL;
		fitLGParameters fp_tDriftR;

		TGraphErrors *MPVplot_xDrift_hits[maxWireGroup];
		TGraphErrors *MPVplot_tDriftL_hits[maxWireGroup];
		TGraphErrors *MPVplot_tDriftR_hits[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){

			MPVplot_xDrift_hits[i-1] = new TGraphErrors(); 
			MPVplot_tDriftL_hits[i-1] = new TGraphErrors();
			MPVplot_tDriftR_hits[i-1] = new TGraphErrors();

		}

		TFile f22((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
		f22.cd();

		for(int i = 1; i <= maxWireGroup; i++){

			for(int j = 1; j <= NBinsX; j++){

				std::cout << "--------- Fitting hit group " + to_string(i) + "; x bin " + to_string(j) + "----------" << std::endl;

				fp_xDrift.reset();
				fp_tDriftL.reset();
				fp_tDriftR.reset();

				projY_xDrift_hits[i-1] = h_dQdx_xDrift_hits[i-1]->ProjectionY(TString::Format("projY_xDrift_%dhits_bin%d", i, j), j, j);
				projY_tDriftL_hits[i-1] = h_dQdx_tDriftL_hits[i-1]->ProjectionY(TString::Format("projY_tDriftL_%dhits_bin%d", i, j), j, j);
				projY_tDriftR_hits[i-1] = h_dQdx_tDriftR_hits[i-1]->ProjectionY(TString::Format("projY_tDriftR_%dhits_bin%d", i, j), j, j);
				
				if((bool)writeProjY){

					projY_xDrift_hits[i-1]->Write();
					projY_tDriftL_hits[i-1]->Write();
					projY_tDriftR_hits[i-1]->Write();

				}

				SetLGParameters(projY_xDrift_hits[i-1], fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
				SetLGParameters(projY_tDriftL_hits[i-1], fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.lb, fp_tDriftL.ub);
				SetLGParameters(projY_tDriftR_hits[i-1], fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.lb, fp_tDriftR.ub);

				LGfit_xDrift_hits[i-1] = fitter(projY_xDrift_hits[i-1], fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
				LGfit_tDriftL_hits[i-1] = fitter(projY_tDriftL_hits[i-1], fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.cov, "LG");
				LGfit_tDriftR_hits[i-1] = fitter(projY_tDriftR_hits[i-1], fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.cov, "LG");

				if(fp_xDrift.fp[1] >= 0.){

					MPVplot_xDrift_hits[i-1]->SetPoint(MPVplot_xDrift_hits[i-1]->GetN(), h_dQdx_xDrift_hits[i-1]->GetXaxis()->GetBinCenter(j), fp_xDrift.fp[1]);
					MPVplot_xDrift_hits[i-1]->SetPointError(MPVplot_xDrift_hits[i-1]->GetN() - 1, 0., pointError(projY_xDrift_hits[i-1], fp_xDrift));
				}

				if(fp_tDriftL.fp[1] >= 0.){

					MPVplot_tDriftL_hits[i-1]->SetPoint(MPVplot_tDriftL_hits[i-1]->GetN(), h_dQdx_tDriftL_hits[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftL.fp[1]);
					MPVplot_tDriftL_hits[i-1]->SetPointError(MPVplot_tDriftL_hits[i-1]->GetN() - 1, 0., pointError(projY_tDriftL_hits[i-1], fp_tDriftL));
				}

				if(fp_tDriftR.fp[1] >= 0.){

					MPVplot_tDriftR_hits[i-1]->SetPoint(MPVplot_tDriftR_hits[i-1]->GetN(), h_dQdx_tDriftR_hits[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftR.fp[1]);
					MPVplot_tDriftR_hits[i-1]->SetPointError(MPVplot_tDriftR_hits[i-1]->GetN() - 1, 0., pointError(projY_tDriftR_hits[i-1], fp_tDriftR));
				}

			}

			MPVplot_xDrift_hits[i-1]->SetName(TString::Format("MPVplot_xDrift_%dhits", i));
			MPVplot_tDriftL_hits[i-1]->SetName(TString::Format("MPVplot_tDriftL_%dhits", i));
			MPVplot_tDriftR_hits[i-1]->SetName(TString::Format("MPVplot_tDriftR_%dhits", i));

			MPVplot_xDrift_hits[i-1]->Write();
			MPVplot_tDriftL_hits[i-1]->Write();
			MPVplot_tDriftR_hits[i-1]->Write();

		}

	}

	if(codeConfig == 3){

		std::string configLabel = "angle";

		TFile f3((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

		TH2D* h_dQdx_xDrift_ang[angBins];
		TH2D* h_dQdx_tDriftL_ang[angBins];
		TH2D* h_dQdx_tDriftR_ang[angBins];

		double** angBinLimits = getAngBinLimits(angBins);

		for(int i = 0; i < angBins/2; i++){

			h_dQdx_xDrift_ang[i] = (TH2D*)f3.Get(TString::Format("h_dQdx_xDrift_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
			h_dQdx_tDriftL_ang[i] = (TH2D*)f3.Get(TString::Format("h_dQdx_tDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
			h_dQdx_tDriftR_ang[i] = (TH2D*)f3.Get(TString::Format("h_dQdx_tDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));

			h_dQdx_xDrift_ang[i + angBins/2] = (TH2D*)f3.Get(TString::Format("h_dQdx_xDrift_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
			h_dQdx_tDriftL_ang[i + angBins/2] = (TH2D*)f3.Get(TString::Format("h_dQdx_tDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
			h_dQdx_tDriftR_ang[i + angBins/2] = (TH2D*)f3.Get(TString::Format("h_dQdx_tDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));

		}

		TH1D* projY_xDrift_ang[angBins];
		TH1D* projY_tDriftL_ang[angBins];
		TH1D* projY_tDriftR_ang[angBins];

		TF1* LGfit_xDrift_ang[angBins];
		TF1* LGfit_tDriftL_ang[angBins];
		TF1* LGfit_tDriftR_ang[angBins];

		fitLGParameters fp0_xDrift, fp1_xDrift;
		fitLGParameters fp0_tDriftL, fp1_tDriftL;
		fitLGParameters fp0_tDriftR, fp1_tDriftR;

		TGraphErrors *MPVplot_xDrift_ang[angBins];
		TGraphErrors *MPVplot_tDriftL_ang[angBins];
		TGraphErrors *MPVplot_tDriftR_ang[angBins];

		for(int i = 0; i < angBins; i++){

			MPVplot_xDrift_ang[i] = new TGraphErrors(); 
			MPVplot_tDriftL_ang[i] = new TGraphErrors();
			MPVplot_tDriftR_ang[i] = new TGraphErrors();

		}

		TFile f33((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
		f33.cd();

		for(int i = 0; i < angBins/2; i++){

			for(int j = 1; j <= NBinsX; j++){

				std::cout << "--------- Fitting angular bin " + to_string(i) + "; x bin " + to_string(j) + "----------" << std::endl;

				fp0_xDrift.reset();
				fp0_tDriftL.reset();
				fp0_tDriftR.reset();

				fp1_xDrift.reset();
				fp1_tDriftL.reset();
				fp1_tDriftR.reset();

				projY_xDrift_ang[i] = h_dQdx_xDrift_ang[i]->ProjectionY(TString::Format("projY_xDrift_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j), j, j);
				projY_tDriftL_ang[i] = h_dQdx_tDriftL_ang[i]->ProjectionY(TString::Format("projY_tDriftL_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j), j, j);
				projY_tDriftR_ang[i] = h_dQdx_tDriftR_ang[i]->ProjectionY(TString::Format("projY_tDriftR_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], j), j, j);
				
				projY_xDrift_ang[i+angBins/2] = h_dQdx_xDrift_ang[i+angBins/2]->ProjectionY(TString::Format("projY_xDrift_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j), j, j);
				projY_tDriftL_ang[i+angBins/2] = h_dQdx_tDriftL_ang[i+angBins/2]->ProjectionY(TString::Format("projY_tDriftL_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j), j, j);
				projY_tDriftR_ang[i+angBins/2] = h_dQdx_tDriftR_ang[i+angBins/2]->ProjectionY(TString::Format("projY_tDriftR_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], j), j, j);

				if((bool)writeProjY){

					projY_xDrift_ang[i]->Write();
					projY_tDriftL_ang[i]->Write();
					projY_tDriftR_ang[i]->Write();

					projY_xDrift_ang[i+angBins/2]->Write();
					projY_tDriftL_ang[i+angBins/2]->Write();
					projY_tDriftR_ang[i+angBins/2]->Write();

				}

				SetLGParameters(projY_xDrift_ang[i], fp0_xDrift.fp, fp0_xDrift.efp, fp0_xDrift.lb, fp0_xDrift.ub);
				SetLGParameters(projY_tDriftL_ang[i], fp0_tDriftL.fp, fp0_tDriftL.efp, fp0_tDriftL.lb, fp0_tDriftL.ub);
				SetLGParameters(projY_tDriftR_ang[i], fp0_tDriftR.fp, fp0_tDriftR.efp, fp0_tDriftR.lb, fp0_tDriftR.ub);

				SetLGParameters(projY_xDrift_ang[i+angBins/2], fp1_xDrift.fp, fp1_xDrift.efp, fp1_xDrift.lb, fp1_xDrift.ub);
				SetLGParameters(projY_tDriftL_ang[i+angBins/2], fp1_tDriftL.fp, fp1_tDriftL.efp, fp1_tDriftL.lb, fp1_tDriftL.ub);
				SetLGParameters(projY_tDriftR_ang[i+angBins/2], fp1_tDriftR.fp, fp1_tDriftR.efp, fp1_tDriftR.lb, fp1_tDriftR.ub);

				LGfit_xDrift_ang[i] = fitter(projY_xDrift_ang[i], fp0_xDrift.lb, fp0_xDrift.ub, fp0_xDrift.fp, fp0_xDrift.efp, fp0_xDrift.cov, "LG");
				LGfit_tDriftL_ang[i] = fitter(projY_tDriftL_ang[i], fp0_tDriftL.lb, fp0_tDriftL.ub, fp0_tDriftL.fp, fp0_tDriftL.efp, fp0_tDriftL.cov, "LG");
				LGfit_tDriftR_ang[i] = fitter(projY_tDriftR_ang[i], fp0_tDriftR.lb, fp0_tDriftR.ub, fp0_tDriftR.fp, fp0_tDriftR.efp, fp0_tDriftR.cov, "LG");

				LGfit_xDrift_ang[i+angBins/2] = fitter(projY_xDrift_ang[i+angBins/2], fp1_xDrift.lb, fp1_xDrift.ub, fp1_xDrift.fp, fp1_xDrift.efp, fp1_xDrift.cov, "LG");
				LGfit_tDriftL_ang[i+angBins/2] = fitter(projY_tDriftL_ang[i+angBins/2], fp1_tDriftL.lb, fp1_tDriftL.ub, fp1_tDriftL.fp, fp1_tDriftL.efp, fp1_tDriftL.cov, "LG");
				LGfit_tDriftR_ang[i+angBins/2] = fitter(projY_tDriftR_ang[i+angBins/2], fp1_tDriftR.lb, fp1_tDriftR.ub, fp1_tDriftR.fp, fp1_tDriftR.efp, fp1_tDriftR.cov, "LG");

				if(fp0_xDrift.fp[1] >= 0.){

					MPVplot_xDrift_ang[i]->SetPoint(MPVplot_xDrift_ang[i]->GetN(), h_dQdx_xDrift_ang[i]->GetXaxis()->GetBinCenter(j), fp0_xDrift.fp[1]);
					MPVplot_xDrift_ang[i]->SetPointError(MPVplot_xDrift_ang[i]->GetN() - 1, 0., pointError(projY_xDrift_ang[i], fp0_xDrift));
				}

				if(fp0_tDriftL.fp[1] >= 0.){

					MPVplot_tDriftL_ang[i]->SetPoint(MPVplot_tDriftL_ang[i]->GetN(), h_dQdx_tDriftL_ang[i]->GetXaxis()->GetBinCenter(j), fp0_tDriftL.fp[1]);
					MPVplot_tDriftL_ang[i]->SetPointError(MPVplot_tDriftL_ang[i]->GetN() - 1, 0., pointError(projY_tDriftL_ang[i], fp0_tDriftL));
				}

				if(fp0_tDriftR.fp[1] >= 0.){

					MPVplot_tDriftR_ang[i]->SetPoint(MPVplot_tDriftR_ang[i]->GetN(), h_dQdx_tDriftR_ang[i]->GetXaxis()->GetBinCenter(j), fp0_tDriftR.fp[1]);
					MPVplot_tDriftR_ang[i]->SetPointError(MPVplot_tDriftR_ang[i]->GetN() - 1, 0., pointError(projY_tDriftR_ang[i], fp0_tDriftR));
				}

				if(fp1_xDrift.fp[1] >= 0.){

					MPVplot_xDrift_ang[i+angBins/2]->SetPoint(MPVplot_xDrift_ang[i+angBins/2]->GetN(), h_dQdx_xDrift_ang[i+angBins/2]->GetXaxis()->GetBinCenter(j), fp1_xDrift.fp[1]);
					MPVplot_xDrift_ang[i+angBins/2]->SetPointError(MPVplot_xDrift_ang[i+angBins/2]->GetN() - 1, 0., pointError(projY_xDrift_ang[i+angBins/2], fp1_xDrift));
				}

				if(fp1_tDriftL.fp[1] >= 0.){

					MPVplot_tDriftL_ang[i+angBins/2]->SetPoint(MPVplot_tDriftL_ang[i+angBins/2]->GetN(), h_dQdx_tDriftL_ang[i+angBins/2]->GetXaxis()->GetBinCenter(j), fp1_tDriftL.fp[1]);
					MPVplot_tDriftL_ang[i+angBins/2]->SetPointError(MPVplot_tDriftL_ang[i+angBins/2]->GetN() - 1, 0., pointError(projY_tDriftL_ang[i+angBins/2], fp1_tDriftL));
				}

				if(fp1_tDriftR.fp[1] >= 0.){

					MPVplot_tDriftR_ang[i+angBins/2]->SetPoint(MPVplot_tDriftR_ang[i+angBins/2]->GetN(), h_dQdx_tDriftR_ang[i+angBins/2]->GetXaxis()->GetBinCenter(j), fp1_tDriftR.fp[1]);
					MPVplot_tDriftR_ang[i+angBins/2]->SetPointError(MPVplot_tDriftR_ang[i+angBins/2]->GetN() - 1, 0., pointError(projY_tDriftR_ang[i+angBins/2], fp1_tDriftR));
				}

			}

			MPVplot_xDrift_ang[i]->SetName(TString::Format("MPVplot_xDrift_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
			MPVplot_tDriftL_ang[i]->SetName(TString::Format("MPVplot_tDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));
			MPVplot_tDriftR_ang[i]->SetName(TString::Format("MPVplot_tDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]));

			MPVplot_xDrift_ang[i+angBins/2]->SetName(TString::Format("MPVplot_xDrift_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
			MPVplot_tDriftL_ang[i+angBins/2]->SetName(TString::Format("MPVplot_tDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));
			MPVplot_tDriftR_ang[i+angBins/2]->SetName(TString::Format("MPVplot_tDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]));

			MPVplot_xDrift_ang[i]->Write();
			MPVplot_tDriftL_ang[i]->Write();
			MPVplot_tDriftR_ang[i]->Write();

			MPVplot_xDrift_ang[i+angBins/2]->Write();
			MPVplot_tDriftL_ang[i+angBins/2]->Write();
			MPVplot_tDriftR_ang[i+angBins/2]->Write();

		}

	}

	if(codeConfig == 4){

		std::string configLabel = "wireAndAngle";

		TFile f4((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str());

		TH2D* h_dQdx_xDrift_angWire[angBins][maxWireGroup];
		TH2D* h_dQdx_tDriftL_angWire[angBins][maxWireGroup];
		TH2D* h_dQdx_tDriftR_angWire[angBins][maxWireGroup];

		double** angBinLimits = getAngBinLimits(angBins);

		for(int i = 0; i < angBins/2; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				h_dQdx_xDrift_angWire[i][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_xDrift_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));
				h_dQdx_tDriftL_angWire[i][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_tDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));
				h_dQdx_tDriftR_angWire[i][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_tDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));

				h_dQdx_xDrift_angWire[i + angBins/2][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_xDrift_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));
				h_dQdx_tDriftL_angWire[i + angBins/2][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_tDriftL_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));
				h_dQdx_tDriftR_angWire[i + angBins/2][j-1] = (TH2D*)f4.Get(TString::Format("h_dQdx_tDriftR_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));

			}

		}

		TH1D* projY_xDrift_angWire[angBins][maxWireGroup];
		TH1D* projY_tDriftL_angWire[angBins][maxWireGroup];
		TH1D* projY_tDriftR_angWire[angBins][maxWireGroup];

		TF1* LGfit_xDrift_angWire[angBins][maxWireGroup];
		TF1* LGfit_tDriftL_angWire[angBins][maxWireGroup];
		TF1* LGfit_tDriftR_angWire[angBins][maxWireGroup];

		fitLGParameters fp0_xDrift, fp1_xDrift;
		fitLGParameters fp0_tDriftL, fp1_tDriftL;
		fitLGParameters fp0_tDriftR, fp1_tDriftR;

		TGraphErrors *MPVplot_xDrift_angWire[angBins][maxWireGroup];
		TGraphErrors *MPVplot_tDriftL_angWire[angBins][maxWireGroup];
		TGraphErrors *MPVplot_tDriftR_angWire[angBins][maxWireGroup];

		for(int i = 0; i < angBins; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				MPVplot_xDrift_angWire[i][j-1] = new TGraphErrors(); 
				MPVplot_tDriftL_angWire[i][j-1] = new TGraphErrors();
				MPVplot_tDriftR_angWire[i][j-1] = new TGraphErrors();

			}

		}

		TFile f44((saveLoc + dataset  + "_" + configLabel + "/MPV_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");
		f44.cd();

		for(int i = 0; i < angBins/2; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				for(int k = 1; k <= NBinsX; k++){

					std::cout << "--------- Fitting angular bin " + to_string(i) + "; wire group " + to_string(j) + "; x bin " + to_string(k) + "----------" << std::endl;

					fp0_xDrift.reset();
					fp0_tDriftL.reset();
					fp0_tDriftR.reset();

					fp1_xDrift.reset();
					fp1_tDriftL.reset();
					fp1_tDriftR.reset();

					projY_xDrift_angWire[i][j-1] = h_dQdx_xDrift_angWire[i][j-1]->ProjectionY(TString::Format("projY_xDrift_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], k), k, k);
					projY_tDriftL_angWire[i][j-1] = h_dQdx_tDriftL_angWire[i][j-1]->ProjectionY(TString::Format("projY_tDriftL_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], k), k, k);
					projY_tDriftR_angWire[i][j-1] = h_dQdx_tDriftR_angWire[i][j-1]->ProjectionY(TString::Format("projY_tDriftR_ang%ito%i_bin%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1], k), k, k);
					
					projY_xDrift_angWire[i+angBins/2][j-1] = h_dQdx_xDrift_angWire[i+angBins/2][j-1]->ProjectionY(TString::Format("projY_xDrift_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], k), k, k);
					projY_tDriftL_angWire[i+angBins/2][j-1] = h_dQdx_tDriftL_angWire[i+angBins/2][j-1]->ProjectionY(TString::Format("projY_tDriftL_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], k), k, k);
					projY_tDriftR_angWire[i+angBins/2][j-1] = h_dQdx_tDriftR_angWire[i+angBins/2][j-1]->ProjectionY(TString::Format("projY_tDriftR_ang%ito%i_bin%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1], k), k, k);

					if((bool)writeProjY){

						projY_xDrift_angWire[i][j-1]->Write();
						projY_tDriftL_angWire[i][j-1]->Write();
						projY_tDriftR_angWire[i][j-1]->Write();

						projY_xDrift_angWire[i+angBins/2][j-1]->Write();
						projY_tDriftL_angWire[i+angBins/2][j-1]->Write();
						projY_tDriftR_angWire[i+angBins/2][j-1]->Write();

					}

					SetLGParameters(projY_xDrift_angWire[i][j-1], fp0_xDrift.fp, fp0_xDrift.efp, fp0_xDrift.lb, fp0_xDrift.ub);
					SetLGParameters(projY_tDriftL_angWire[i][j-1], fp0_tDriftL.fp, fp0_tDriftL.efp, fp0_tDriftL.lb, fp0_tDriftL.ub);
					SetLGParameters(projY_tDriftR_angWire[i][j-1], fp0_tDriftR.fp, fp0_tDriftR.efp, fp0_tDriftR.lb, fp0_tDriftR.ub);

					SetLGParameters(projY_xDrift_angWire[i+angBins/2][j-1], fp1_xDrift.fp, fp1_xDrift.efp, fp1_xDrift.lb, fp1_xDrift.ub);
					SetLGParameters(projY_tDriftL_angWire[i+angBins/2][j-1], fp1_tDriftL.fp, fp1_tDriftL.efp, fp1_tDriftL.lb, fp1_tDriftL.ub);
					SetLGParameters(projY_tDriftR_angWire[i+angBins/2][j-1], fp1_tDriftR.fp, fp1_tDriftR.efp, fp1_tDriftR.lb, fp1_tDriftR.ub);

					LGfit_xDrift_angWire[i][j-1] = fitter(projY_xDrift_angWire[i][j-1], fp0_xDrift.lb, fp0_xDrift.ub, fp0_xDrift.fp, fp0_xDrift.efp, fp0_xDrift.cov, "LG");
					LGfit_tDriftL_angWire[i][j-1] = fitter(projY_tDriftL_angWire[i][j-1], fp0_tDriftL.lb, fp0_tDriftL.ub, fp0_tDriftL.fp, fp0_tDriftL.efp, fp0_tDriftL.cov, "LG");
					LGfit_tDriftR_angWire[i][j-1] = fitter(projY_tDriftR_angWire[i][j-1], fp0_tDriftR.lb, fp0_tDriftR.ub, fp0_tDriftR.fp, fp0_tDriftR.efp, fp0_tDriftR.cov, "LG");

					LGfit_xDrift_angWire[i+angBins/2][j-1] = fitter(projY_xDrift_angWire[i+angBins/2][j-1], fp1_xDrift.lb, fp1_xDrift.ub, fp1_xDrift.fp, fp1_xDrift.efp, fp1_xDrift.cov, "LG");
					LGfit_tDriftL_angWire[i+angBins/2][j-1] = fitter(projY_tDriftL_angWire[i+angBins/2][j-1], fp1_tDriftL.lb, fp1_tDriftL.ub, fp1_tDriftL.fp, fp1_tDriftL.efp, fp1_tDriftL.cov, "LG");
					LGfit_tDriftR_angWire[i+angBins/2][j-1] = fitter(projY_tDriftR_angWire[i+angBins/2][j-1], fp1_tDriftR.lb, fp1_tDriftR.ub, fp1_tDriftR.fp, fp1_tDriftR.efp, fp1_tDriftR.cov, "LG");

					if(fp0_xDrift.fp[1] >= 0.){

						MPVplot_xDrift_angWire[i][j-1]->SetPoint(MPVplot_xDrift_angWire[i][j-1]->GetN(), h_dQdx_xDrift_angWire[i][j-1]->GetXaxis()->GetBinCenter(k), fp0_xDrift.fp[1]);
						MPVplot_xDrift_angWire[i][j-1]->SetPointError(MPVplot_xDrift_angWire[i][j-1]->GetN() - 1, 0., pointError(projY_xDrift_angWire[i][j-1], fp0_xDrift));
					}

					if(fp0_tDriftL.fp[1] >= 0.){

						MPVplot_tDriftL_angWire[i][j-1]->SetPoint(MPVplot_tDriftL_angWire[i][j-1]->GetN(), h_dQdx_tDriftL_angWire[i][j-1]->GetXaxis()->GetBinCenter(k), fp0_tDriftL.fp[1]);
						MPVplot_tDriftL_angWire[i][j-1]->SetPointError(MPVplot_tDriftL_angWire[i][j-1]->GetN() - 1, 0., pointError(projY_tDriftL_angWire[i][j-1], fp0_tDriftL));
					}

					if(fp0_tDriftR.fp[1] >= 0.){

						MPVplot_tDriftR_angWire[i][j-1]->SetPoint(MPVplot_tDriftR_angWire[i][j-1]->GetN(), h_dQdx_tDriftR_angWire[i][j-1]->GetXaxis()->GetBinCenter(k), fp0_tDriftR.fp[1]);
						MPVplot_tDriftR_angWire[i][j-1]->SetPointError(MPVplot_tDriftR_angWire[i][j-1]->GetN() - 1, 0., pointError(projY_tDriftR_angWire[i][j-1], fp0_tDriftR));
					}

					if(fp1_xDrift.fp[1] >= 0.){

						MPVplot_xDrift_angWire[i+angBins/2][j-1]->SetPoint(MPVplot_xDrift_angWire[i+angBins/2][j-1]->GetN(), h_dQdx_xDrift_angWire[i+angBins/2][j-1]->GetXaxis()->GetBinCenter(k), fp1_xDrift.fp[1]);
						MPVplot_xDrift_angWire[i+angBins/2][j-1]->SetPointError(MPVplot_xDrift_angWire[i+angBins/2][j-1]->GetN() - 1, 0., pointError(projY_xDrift_angWire[i+angBins/2][j-1], fp1_xDrift));
					}

					if(fp1_tDriftL.fp[1] >= 0.){

						MPVplot_tDriftL_angWire[i+angBins/2][j-1]->SetPoint(MPVplot_tDriftL_angWire[i+angBins/2][j-1]->GetN(), h_dQdx_tDriftL_angWire[i+angBins/2][j-1]->GetXaxis()->GetBinCenter(k), fp1_tDriftL.fp[1]);
						MPVplot_tDriftL_angWire[i+angBins/2][j-1]->SetPointError(MPVplot_tDriftL_angWire[i+angBins/2][j-1]->GetN() - 1, 0., pointError(projY_tDriftL_angWire[i+angBins/2][j-1], fp1_tDriftL));
					}

					if(fp1_tDriftR.fp[1] >= 0.){

						MPVplot_tDriftR_angWire[i+angBins/2][j-1]->SetPoint(MPVplot_tDriftR_angWire[i+angBins/2][j-1]->GetN(), h_dQdx_tDriftR_angWire[i+angBins/2][j-1]->GetXaxis()->GetBinCenter(k), fp1_tDriftR.fp[1]);
						MPVplot_tDriftR_angWire[i+angBins/2][j-1]->SetPointError(MPVplot_tDriftR_angWire[i+angBins/2][j-1]->GetN() - 1, 0., pointError(projY_tDriftR_angWire[i+angBins/2][j-1], fp1_tDriftR));
					}

				}

				MPVplot_xDrift_angWire[i][j-1]->SetName(TString::Format("MPVplot_xDrift_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));
				MPVplot_tDriftL_angWire[i][j-1]->SetName(TString::Format("MPVplot_tDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));
				MPVplot_tDriftR_angWire[i][j-1]->SetName(TString::Format("MPVplot_tDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j));

				MPVplot_xDrift_angWire[i+angBins/2][j-1]->SetName(TString::Format("MPVplot_xDrift_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));
				MPVplot_tDriftL_angWire[i+angBins/2][j-1]->SetName(TString::Format("MPVplot_tDriftL_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));
				MPVplot_tDriftR_angWire[i+angBins/2][j-1]->SetName(TString::Format("MPVplot_tDriftR_ang%ito%inWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j));

				MPVplot_xDrift_angWire[i][j-1]->Write();
				MPVplot_tDriftL_angWire[i][j-1]->Write();
				MPVplot_tDriftR_angWire[i][j-1]->Write();

				MPVplot_xDrift_angWire[i+angBins/2][j-1]->Write();
				MPVplot_tDriftL_angWire[i+angBins/2][j-1]->Write();
				MPVplot_tDriftR_angWire[i+angBins/2][j-1]->Write();

			}

		}

	}

	return 0;
}
