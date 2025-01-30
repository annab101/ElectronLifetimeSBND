//Author: Anna Beever
//Date:   April 2024

//C++ includes
#include<fstream>
#include<filesystem>
#include<valarray>

//ROOT includes
#include "../Utilities/ROOTincludes.h"

//Local includes
#include "../Helpers/Constants.h"
#include "../Helpers/PlottingHelpers.h"
#include "../Helpers/PlottingHelpers.cpp"
#include "../Helpers/FittingHelpers.h"
#include "../Helpers/FittingHelpers.cpp"
#include "../Utilities/ConfigReader.h"
#include "../Utilities/ConfigReader.cpp"

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
	
	int nGroupedWires = 0;
	int NBinsX = 100;
	int NBinsT = 100;
	int NBinsdQdx = 75;
	double minX = -200.;
	double maxX = 200.;
	double mindQdx = 200.;
	double maxdQdx = 1800.;
	double minT = 0.;
	double maxT = 1.3;	
	int writeProjY = 0;
	std::string configLabel = "noConfigLabel";

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("nGroupedWires", nGroupedWires);
	p->getValue("NBinsX", NBinsX);
	p->getValue("NBinsT", NBinsT);
	p->getValue("NBinsdQdx", NBinsdQdx);
	p->getValue("minX", minX);
	p->getValue("maxX", maxX);
	p->getValue("mindQdx", mindQdx);
	p->getValue("maxdQdx", maxdQdx);
	p->getValue("minT", minT);
	p->getValue("maxT", maxT);
	p->getValue("writeProjY", writeProjY);

	TFile f((saveLoc + dataset + "/dQdx_hist_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str());
	TFile ff((saveLoc + dataset  + "/MPV_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str(), "new");


	TH2D* h_dQdx_xDrift = (TH2D*)f.Get(TString::Format("h_dQdx_xDrift_%dwires", nGroupedWires));
	TH2D* h_dQdx_tDriftE = (TH2D*)f.Get(TString::Format("h_dQdx_tDriftE_%dwires", nGroupedWires));
	TH2D* h_dQdx_tDriftW = (TH2D*)f.Get(TString::Format("h_dQdx_tDriftW_%dwires", nGroupedWires));
	
	TH1D* projY_xDrift = new TH1D();
	TH1D* projY_tDriftE = new TH1D();
	TH1D* projY_tDriftW = new TH1D();

	TF1* LGfit_xDrift = new TF1();
	TF1* LGfit_tDriftE = new TF1();
	TF1* LGfit_tDriftW = new TF1();

	fitLGParameters fp_xDrift;
	fitLGParameters fp_tDriftL;
	fitLGParameters fp_tDriftR;

	TGraphAsymmErrors *MPVplot_xDrift = new TGraphAsymmErrors();
	TGraphAsymmErrors *MPVplot_tDriftE = new TGraphAsymmErrors();
	TGraphAsymmErrors *MPVplot_tDriftW = new TGraphAsymmErrors();

	ff.cd();

	for(int i = 1; i <= NBinsX; i++){

		std::cout << "--------- Fitting bin " + to_string(i) + "----------" << std::endl;

		projY_xDrift = h_dQdx_xDrift->ProjectionY(TString::Format("projY_xDrift_%dwires_bin%d", nGroupedWires, i), i, i);
		projY_tDriftE = h_dQdx_tDriftE->ProjectionY(TString::Format("projY_tDriftE_%dwires_bin%d", nGroupedWires, i), i, i);
		projY_tDriftW = h_dQdx_tDriftW->ProjectionY(TString::Format("projY_tDriftW_%dwires_bin%d", nGroupedWires, i), i, i);

		if((bool)writeProjY){

			projY_xDrift->Write();
			projY_tDriftE->Write();
			projY_tDriftW->Write();

		}

		SetLGParameters(projY_xDrift, fp_xDrift.fp, fp_xDrift.ehfp, fp_xDrift.elfp, fp_xDrift.lb, fp_xDrift.ub);
		SetLGParameters(projY_tDriftE, fp_tDriftL.fp, fp_tDriftL.ehfp, fp_tDriftL.elfp, fp_tDriftL.lb, fp_tDriftL.ub);
		SetLGParameters(projY_tDriftW, fp_tDriftR.fp, fp_tDriftR.ehfp, fp_tDriftR.elfp, fp_tDriftR.lb, fp_tDriftR.ub);

		LGfit_xDrift = fitter(projY_xDrift, fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.ehfp, fp_xDrift.elfp, fp_xDrift.cov, "LG");
		LGfit_tDriftE = fitter(projY_tDriftE, fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.ehfp, fp_tDriftL.elfp, fp_tDriftL.cov, "LG");
		LGfit_tDriftW = fitter(projY_tDriftW, fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.ehfp, fp_tDriftR.elfp, fp_tDriftR.cov, "LG");

		if(fp_xDrift.fp[1] >= 0.){

			MPVplot_xDrift->SetPoint(MPVplot_xDrift->GetN(), h_dQdx_xDrift->GetXaxis()->GetBinCenter(i), fp_xDrift.fp[1]);
			MPVplot_xDrift->SetPointError(MPVplot_xDrift->GetN() - 1, 0., 0., -1*fp_xDrift.elfp[1], fp_xDrift.ehfp[1]);
		}

		if(fp_tDriftL.fp[1] >= 0.){

			MPVplot_tDriftE->SetPoint(MPVplot_tDriftE->GetN(), h_dQdx_tDriftE->GetXaxis()->GetBinCenter(i), fp_tDriftL.fp[1]);
			MPVplot_tDriftE->SetPointError(MPVplot_tDriftE->GetN() - 1, 0., 0., -1*fp_tDriftL.elfp[1], fp_tDriftL.ehfp[1]);
		}

		if(fp_tDriftR.fp[1] >= 0.){

			MPVplot_tDriftW->SetPoint(MPVplot_tDriftW->GetN(), h_dQdx_tDriftW->GetXaxis()->GetBinCenter(i), fp_tDriftR.fp[1]);
			MPVplot_tDriftW->SetPointError(MPVplot_tDriftW->GetN() - 1, 0., 0., -1*fp_tDriftR.elfp[1], fp_tDriftR.ehfp[1]);
		}

		if((bool)writeProjY){

			LGfit_xDrift->SetName(TString::Format("LGfit_xDrift_%dwires_bin%d",nGroupedWires,i));
			LGfit_tDriftE->SetName(TString::Format("LGfit_tDriftE_%dwires_bin%d",nGroupedWires,i));
			LGfit_tDriftW->SetName(TString::Format("LGfit_tDriftW_%dwires_bin%d",nGroupedWires,i));

			LGfit_xDrift->Write();
			LGfit_tDriftE->Write();
			LGfit_tDriftW->Write();

		}

	}

	MPVplot_xDrift->SetName(TString::Format("MPVplot_xDrift_%dwires",nGroupedWires));
	MPVplot_tDriftE->SetName(TString::Format("MPVplot_tDriftE_%dwires",nGroupedWires));
	MPVplot_tDriftW->SetName(TString::Format("MPVplot_tDriftW_%dwires",nGroupedWires));

	MPVplot_xDrift->Write();
	MPVplot_tDriftE->Write();
	MPVplot_tDriftW->Write();

	

	return 0;
}
