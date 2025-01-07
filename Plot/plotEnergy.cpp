//Author: Anna Beever
//Date:   February 2024
//Plots energy dist

//C++ includes
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include <chrono>
#include <ctime>    

//ROOT includes
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TAxis.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TError.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TMultiGraph.h"

//Function definitions
void BinLogX(TH1 *hist);

int main(int argc, char**argv) {
	
	//filename is input data file, fileSaveID is (partial) name of output files
	std::string filename = "no file";
	std::string fileSaveID = "noID";

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--filename")){
			filename = argv[i+1];
		}
		if(!strcmp(argv[i], "--fileSaveID")){
			fileSaveID = argv[i+1];
		}
	}
	
	//check input files have been provided
	if(!strcmp("no file", filename.c_str()) || !strcmp("noID", fileSaveID.c_str())){
		std::cout << "!!!!ERROR: MISSING INPUT INFO (file name and ID for saving)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		exit(0);
	}
	
	std::cout << "----Retrieving file----------------------------------------------------" << std::endl;

	//define names of plots
	std::string mydata_cosmics = "/exp/sbnd/data/users/abeever/cosmics_analysis/";

	std::string energyPlotName = "energyTruth_hist_" + fileSaveID + ".png";	
	std::string energySaveLoc = mydata_cosmics + energyPlotName;


	//Read in data files line by line and chain together
	TChain *chain = new TChain;

	std::fstream file_list;
	file_list.open(filename.c_str(),std::ios::in);
	
	if (file_list.is_open()){
		std::string fileUrl;
		while(getline(file_list, fileUrl)){
			//std::cout << fileUrl << std::endl;
			chain->Add((fileUrl + "/caloskim/TrackCaloSkim").c_str());
			
		}
		file_list.close();
	}

	//read in variables of interest
	TTreeReader treereader(chain);
	TTreeReaderValue<Float_t> read_xi(treereader, "trk.start.x");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<Float_t> read_E(treereader, "trk.truth.p.startE");

    std::cout << "----Plotting energy----------------------------------------------------" << std::endl;
    
    //variables of interest
    std::vector<double> E;	
    std::vector<double> weights; //same weights as for angular plots - need a vector of 1s

    //read in data from ROOT files
    treereader.Restart();
    while(treereader.Next()){	
        
        if(*read_xi < -199. || *read_xi > 199. || *read_xf < -199. || *read_xf > 199.){//AC crosser filter
            E.push_back(*read_E);   
            weights.push_back(1.); 
        }				                 
        
    }

    //Plot energy histogram
    TH1D *h5 = new TH1D("h5","Energy (Truth)", 100, 0.1, 10000);
    BinLogX(h5); //logarithmic bins
    TCanvas *c5 = new TCanvas();
    h5->FillN(weights.size(),E.data(), weights.data());
    c5->cd();
    c5->SetLogy();
    c5->SetLogx();
    h5->GetXaxis()->SetTitle("E (GeV)");
    double XbinWidth = h5->GetXaxis()->GetBinWidth(1);
    h5->GetYaxis()->SetTitle("Entries");
    //h5->SetStats(0);
    h5->Draw();
    c5->SaveAs(energySaveLoc.c_str());

    return 0;
}

//Function for making logarithmic bins for a logarithmic axis
void BinLogX(TH1 *hist){

	//Based on code by Roland Kuhn, Technical University of Munich, 2006
   TAxis *axis = hist->GetXaxis();
   int bins = axis->GetNbins();

   double_t xmin = log10(axis->GetXmin());
   double_t xmax = log10(axis->GetXmax());
   double_t width = (xmax - xmin) / bins;
   double_t *new_bins = new double_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = pow(10, xmin + i * width);
   }
   axis->Set(bins, new_bins);
   delete new_bins;
}