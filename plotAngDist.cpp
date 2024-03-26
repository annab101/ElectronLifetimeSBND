//Author: Anna Beever
//Date:   February 2024
//Does angular plots

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
void drawAngleLimits(std::vector<double> bin_limits, TH1 *hist);
void histLabels1D(TH1D *h, std::string units, int precision, bool brackets);

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

	std::string azimuthPlotName = "azimuth_hist_" + fileSaveID + ".png";
	std::string zenithPlotName = "zenith_hist_" + fileSaveID + ".png";
	std::string cosZenithPlotName = "cosZenith_hist_" + fileSaveID + ".png";
	
	std::string azimuthSaveLoc = mydata_cosmics + azimuthPlotName;
	std::string zenithSaveLoc = mydata_cosmics + zenithPlotName;
	std::string cosZenithSaveLoc = mydata_cosmics + cosZenithPlotName;

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
	TTreeReaderValue<Float_t> read_yi(treereader, "trk.start.y");
	TTreeReaderValue<Float_t> read_zi(treereader, "trk.start.z");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<Float_t> read_yf(treereader, "trk.end.y");
	TTreeReaderValue<Float_t> read_zf(treereader, "trk.end.z");
	
    //variables of interest
    std::vector<double> xi;
    std::vector<double> yi;
    std::vector<double> zi;
    std::vector<double> xf;
    std::vector<double> yf;
    std::vector<double> zf;
    std::vector<double> weights; //couldn't find a way of plotting TH1D from a vector of values without a vector of weights,
									 //so just make a weights vector full of 1s

		//fill variables from ROOT data files
		treereader.Restart();
		while(treereader.Next()){
					
            if(*read_xi < -199. || *read_xi > 199. || *read_xf < -199. || *read_xf > 199.){ //AC crosser filter
                weights.push_back(1.);
                xi.push_back(*read_xi);
                yi.push_back(*read_yi);
                zi.push_back(*read_zi);
                xf.push_back(*read_xf);
                yf.push_back(*read_yf);
                zf.push_back(*read_zf);
            }
			
        }	

		std::cout << "AC track count 2: " << xi.size() << std::endl;
		
		//Initialise vectors for angle calculations
		std::vector<double> tanAzimuth;
		std::vector<double> tanZenith;
		std::vector<double> cosZenith;
		std::vector<double> azimuth;
		std::vector<double> zenith;

		std::cout << "----Calculating angular data-------------------------------------------" << std::endl;
		//Calculate angle data
		for(int i = 0; i < weights.size(); i++){
			tanZenith.push_back((sqrt(pow((xf[i]-xi[i]),2)+pow((zf[i]-zi[i]),2)))/(yf[i]-yi[i]));
			cosZenith.push_back(sqrt(1/(1+pow(tanZenith[i],2))));
			azimuth.push_back(atan2((zf[i]-zi[i]),(xf[i]-xi[i]))/M_PI*180);
			zenith.push_back(acos(cosZenith[i])/M_PI*180);
		}

		std::cout << "----Plotting azimuthal angle-------------------------------------------" << std::endl;
		//Azimuthal angle plot

		TH1D *h2 = new TH1D("h2","$\\phi$", 100, -180, 180);
		TCanvas *c2 = new TCanvas();
		h2->FillN(weights.size(),azimuth.data(), weights.data());
		c2->cd();
		//h2->SetStats(0); //turn stats box off
		h2->GetXaxis()->SetTitle("$\\phi$ $(^\\circ)$");
		histLabels1D(h2, "^\\circ", 3, 1);
		h2->Draw();


		//plot angular limits on the same plot
		std::vector<double> azimuth_limits = {-68.2,68.2,-111.8,111.8};		
		drawAngleLimits(azimuth_limits, h2);
		c2->SaveAs(azimuthSaveLoc.c_str());

		std::cout << "----Plotting zenith angle----------------------------------------------" << std::endl;
		//zenith angle plot

		TH1D *h3 = new TH1D("h3","$\\theta$", 100, 0, 90);
		TCanvas *c3 = new TCanvas();
		h3->FillN(weights.size(),zenith.data(), weights.data());
		c3->cd();
		h3->GetXaxis()->SetTitle("$\\theta$ $(^\\circ)$");
		histLabels1D(h3, "^\\circ", 3, 1);
		h3->Draw();

		//plot angular limits on the same plot
		std::vector<double> zenith_limits = {26.6};		
		drawAngleLimits(zenith_limits, h3);
		c3->SaveAs(zenithSaveLoc.c_str());

		std::cout << "----Plotting cosine of zenith angle------------------------------------" << std::endl;
		//cosine of zenith angle plot

		TH1D *h4 = new TH1D("h4","$\\cos \\theta$", 100, 0, 1);
		TCanvas *c4 = new TCanvas();
		h4->FillN(weights.size(),cosZenith.data(), weights.data());
		c4->cd();
		c4->SetLogy();
		h4->GetXaxis()->SetTitle("$\\cos \\theta$");
		histLabels1D(h4, "", 3, 0);
		h4->Draw();

		//plot angular limits on same plot
		std::vector<double> cosZenith_limits = {0.89};		
		drawAngleLimits(cosZenith_limits, h4);
		c4->SaveAs(cosZenithSaveLoc.c_str());
                
    return 0; 
}

//Function for drawing vertical red dashed lines on a histogram
void drawAngleLimits(std::vector<double> bin_limits, TH1 *hist){
	double max = hist->GetMaximum();
	TLine *line = new TLine;
	for(int i=0; i<bin_limits.size();i++){
		line->SetLineColor(2);
		line->SetLineStyle(9);
		line->SetLineWidth(1);
		line->DrawLine(bin_limits[i],0,bin_limits[i],max);
	}
}

void histLabels1D(TH1D *h, std::string units, int precision, bool brackets){
	std::stringstream ss;
	ss << std::setprecision(precision) << h->GetXaxis()->GetBinWidth(1);
	std::string yLabel;
	if(brackets){yLabel = "Entries / (" + ss.str() + " " + units + ")";}
	else{yLabel = "Entries / " + ss.str() + " " + units;}
	h->GetYaxis()->SetTitle(yLabel.c_str());
}
	