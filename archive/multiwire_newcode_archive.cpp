//Author: Anna Beever
//Date:   December 2023
//Takes fileName.list as input, collects data with TChain and then
//does lifetime calc or dQ/dx, energy and angular plots as selected.

//C++ includes
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<chrono>
#include<ctime>   
#include<valarray> 

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


//Classes

class fitLGParameters {
	
	public: 
		

		double fp[4] = {0.}; //fit parameters
		double efp[4] = {0.}; //fit parameter errors
		double cov[4] = {0.}; //covariance matrix
		double lb = -1.; //fit lower bound
		double ub = -1.;

		void reset() {
			std::fill(std::begin(fp), std::end(fp), 0.);
			std::fill(std::begin(efp), std::end(efp), 0.);
			std::fill(std::begin(cov), std::end(cov), 0.);
			lb = -1.;
			ub = -1.;
		}

};

class fitExpoParameters {
	
	public: 
		

		double fp[2] = {0.}; //fit parameters
		double efp[2] = {0.}; //fit parameter errors
		double cov[2] = {0.}; //covariance matrix
		double lb = -1.; //fit lower bound
		double ub = -1.;

		void reset() {
			std::fill(std::begin(fp), std::end(fp), 0.);
			std::fill(std::begin(efp), std::end(efp), 0.);
			std::fill(std::begin(cov), std::end(cov), 0.);
			lb = -1.;
			ub = -1.;
		}
		
};

//Function definitions
int bool_input(std::string bool_name, int arg1, char**arg2);
double_t langaufun(Double_t *x, Double_t *par);
double_t expofunX(Double_t *x, Double_t *par);
double_t expofunT(Double_t *x, Double_t *par);
void SetLGParameters(TH1D *h, Double_t *fp, Double_t *efp, Double_t &lb, Double_t &ub);
void SetExpoParameters(Double_t *fp);
TF1 *fitter(TH1D *h, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, std::string funcName);
TF1 *fitter(TGraphErrors *g, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, std::string funcName);
double pointError(TH1D *proj_y, fitLGParameters fitParams);
void setMarker(TH1D *h, Color_t color, Style_t style, Size_t size);
void setMarker(TGraphErrors *g, Color_t color, Style_t style, Size_t size);
void saveFig(TCanvas *c, std::string plotName);
void openAndClear(TCanvas *c);
void addPointAndError(TGraphErrors *g, double x, double y, double errx, double erry);
void addPointAndError(TGraphErrors *g, int x, double y, double errx, double erry);
TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TF1* f2);
TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1);
TPaveText *statsBox(std::vector<double> pos, int track_count);


int main(int argc, char**argv) {

	auto start = std::chrono::system_clock::now();
	
	//Set options at command line. Note MPV not possible without dQdx
	//bool plot_dQdx = bool_input("--dQdx", argc, argv);
	//bool plot_projY = bool_input("--projY", argc, argv);
	
	//filename is input data file, fileSaveID is (partial) name of output files
	std::string filename = "no file";
	std::string fileSaveID = "noID";
	std::string fileSaveLoc = "noLoc";

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--filename")){
			filename = argv[i+1];
		}
		if(!strcmp(argv[i], "--fileSaveID")){
			fileSaveID = argv[i+1];
		}
		if(!strcmp(argv[i], "--fileSaveLoc")){
			fileSaveLoc = argv[i+1];
		}
	}
	
	//check input files have been provided
	if(!strcmp("no file", filename.c_str()) || !strcmp("noID", fileSaveID.c_str()) || !strcmp("noLoc", fileSaveLoc.c_str())){
		std::cout << "!!!!ERROR: MISSING INPUT INFO (file name and ID for saving)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		exit(0);
	}
	
	std::cout << "----Retrieving file----------------------------------------------------" << std::endl;

	//define names of plots
	std::string mydata_cosmics = "/exp/sbnd/data/users/abeever/cosmics_analysis/";

	auto preamble_time = std::chrono::system_clock::now();
 	std::chrono::duration<double> elapsed_seconds = preamble_time-start;
 
    std::cout << "++++++++++ Preamble time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

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
	TTreeReaderArray<Float_t> read_x(treereader, "trk.hits2.h.sp.x");
	TTreeReaderArray<Float_t> read_y(treereader, "trk.hits2.h.sp.y");
	TTreeReaderArray<Float_t> read_z(treereader, "trk.hits2.h.sp.z");
	TTreeReaderArray<Float_t> read_dqdx(treereader, "trk.hits2.dqdx");
	TTreeReaderArray<Float_t> read_T(treereader, "trk.hits2.h.time");
	TTreeReaderValue<Float_t> read_t0(treereader, "trk.t0");
	TTreeReaderValue<Float_t> read_xi(treereader, "trk.start.x");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<int> read_selected(treereader, "trk.selected");
	TTreeReaderArray<uint16_t> read_wire(treereader, "trk.hits2.h.wire");

	std::chrono::duration<double> iteratorFullTime = start-start;
	std::chrono::duration<double> filterandfillFullTime = start-start;
	std::chrono::duration<double> loopDuration = start-start;
	std::chrono::duration<double> loopOver_time = start-start;
	auto filterAndFill_time = std::chrono::system_clock::now();
	int whileCount = 0;

	auto TChain_time = std::chrono::system_clock::now();
 	elapsed_seconds = TChain_time-preamble_time; 
    std::cout << "++++++++++ TChain read in time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	//Lifetime for each value of N
	std::vector<std::vector<double>> tau_xDrift_LR {{},{}};
	std::vector<std::vector<double>> errTau_xDrift_LR {{},{}};
	std::vector<std::vector<double>> tau_tDrift_LR {{},{}};
	std::vector<std::vector<double>> errTau_tDrift_LR {{},{}};
	std::vector<int> N_wires;
	TGraphErrors *g_lifeVwiresXL = new TGraphErrors;
	TGraphErrors *g_lifeVwiresXR = new TGraphErrors;
	TGraphErrors *g_lifeVwiresTL = new TGraphErrors;
	TGraphErrors *g_lifeVwiresTR = new TGraphErrors;

	int maxWireGroup = 20;
	int N_bins = 100;
	int track_count = 0;

	TH2D** h_dQdx_xDrift = new TH2D*[maxWireGroup];
	TH2D** h_dQdx_tDriftL = new TH2D*[maxWireGroup];
	TH2D** h_dQdx_tDriftR = new TH2D*[maxWireGroup];

	for(int i = 1; i <= maxWireGroup; i++){
		h_dQdx_xDrift[i-1] = new TH2D(TString::Format("h_dQdx_xDrift_%dwires", i),"dQ/dx vs x", N_bins, -200, 200, 75, 200, 1800);
		h_dQdx_tDriftL[i-1] = new TH2D(TString::Format("h_dQdx_tDriftL_%dwires", i),"dQ/dx vs t Left TPC", N_bins, 0, 1.3, 75, 200, 1800);
		h_dQdx_tDriftR[i-1] = new TH2D(TString::Format("h_dQdx_tDriftR_%dwires", i),"dQ/dx vs t Right TPC", N_bins, 0, 1.3, 75, 200, 1800);
	}

	auto DeclHist_time = std::chrono::system_clock::now();
 	elapsed_seconds = DeclHist_time - TChain_time; 
    std::cout << "++++++++++ Declaring hists time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	treereader.Restart();
	while(treereader.Next()){

		if(*read_selected != 1){
			continue;
		}

		std::cout << "track_count: " << track_count << std::endl;
		track_count ++;
		
		auto iteratorStart_time = std::chrono::system_clock::now();
		if(whileCount > 0){loopOver_time += iteratorStart_time - filterAndFill_time;}

		
		auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
		int start_index = std::distance(read_x.begin(), it_start_x);

		if(start_index == read_x.GetSize()){continue;}
		auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
		int end_index = std::distance(read_x.begin(), it_end_x);
		if(end_index < 0){continue;}

		
		auto iteratorEnd_time = std::chrono::system_clock::now();

		iteratorFullTime += iteratorEnd_time - iteratorStart_time;

		int minWire = read_wire[start_index];
		int maxWire = read_wire[end_index];

		if(minWire > maxWire){
			std::swap(minWire, maxWire);
		}

		std::valarray<double> dQdx_sum(0.,maxWireGroup);
		std::valarray<double> x_sum(0.,maxWireGroup);
		std::valarray<double> t_sum(0.,maxWireGroup);
		std::valarray<int> count(0,maxWireGroup);

		
		for(int i = start_index; i <= end_index; i++){

			//Initial thoughts for removing SCE (in truth)
			/*geo::Point_t spacePoint = {read_x[i], read_y[i], read_z[i]};
			geo::Vector_t SCEoffset = SCE->GetCalPosOffsets(spacePoint);
			read_x[i] = read_x[i] - SCEoffset.X();*/

			read_x[i] = read_x[i];

			dQdx_sum += read_dqdx[i];
			x_sum += read_x[i];
			t_sum += read_T[i];
			count += 1;

			if( i < end_index && (read_x[i] * (-read_x[i+1])) < 0.){ //changed for TPC flip experiment

				for (int N = 1; N <= maxWireGroup; N++){
				
					int j = N-1;

					h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
									
					if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
						h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
					}

					if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
						h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
					}
				}				

				dQdx_sum *= 0.;
				x_sum *= 0.;
				t_sum *= 0.;
				count *= 0;

			}
			else if(i == end_index){

				for (int N = 1; N <= maxWireGroup; N++){
				
					int j = N-1;

					h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
						
					if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
						h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
					}

					if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
						h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
					}
				}	

				dQdx_sum *= 0.;
				x_sum *= 0.;
				t_sum *= 0.;
				count *= 0;

			}
			/*else{

				for (int N = 1; N <= maxWireGroup; N++){
					
					if((read_wire[i] - minWire)%N == 0){

						int j = N-1;

						h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
						
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						dQdx_sum [j] = 0.;
						x_sum [j] = 0.;
						t_sum [j] = 0.;
						count [j] = 0;

					}
				}

			}
			else{

				for (int N = 1; N <= maxWireGroup; N++){

					int if_count = 0;
					
					if((read_wire[i] - minWire)/N != (read_wire[i+1] - minWire)/N){

						int m;

						if(end_index - i >= 5){
							m = 5;
						}
						else{
							m = end_index - i;
						}

						for(int k = i + 1; k <= i + m; k++){

							if((read_wire[k] - minWire)/N == (read_wire[i] - minWire)/N){

								dQdx_sum += read_dqdx[k];
								x_sum += read_x[k];
								t_sum += read_T[k];
								count += 1;
								if_count += 1;
								std::swap(read_wire[k-1],read_wire[k]);
								std::swap(read_dqdx[k-1],read_dqdx[k]);
								std::swap(read_x[k-1],read_x[k]);
								std::swap(read_T[k-1],read_T[k]);

							}

						}

						i = i + if_count;

						std::cout << "if count: " << if_count << std::endl;
						std::cout << "i: " << i << std::endl;

						int j = N-1;

						h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
						
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						dQdx_sum [j] = 0.;
						x_sum [j] = 0.;
						t_sum [j] = 0.;
						count [j] = 0;

					}
				}

			}*/
			else{

				for (int N = 1; N <= maxWireGroup; N++){
					
					if((read_wire[i] - minWire)/N != (read_wire[i+1] - minWire)/N){

						int j = N-1;

						h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
						
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						dQdx_sum [j] = 0.;
						x_sum [j] = 0.;
						t_sum [j] = 0.;
						count [j] = 0;

					}
				}

			}

		
		}

		
		filterAndFill_time = std::chrono::system_clock::now();
		filterandfillFullTime += filterAndFill_time - iteratorEnd_time;
		loopDuration += filterAndFill_time - iteratorStart_time;

		whileCount ++;

	}

	
	std::cout << "++++++++++ Iterator time: " << iteratorFullTime.count() << std::endl;
	std::cout << "++++++++++ Filter & fill time: " << filterandfillFullTime.count() << std::endl;
	std::cout << "++++++++++ Loop duration time: " << loopDuration.count() << std::endl;
	std::cout << "++++++++++ Looping time: " << loopOver_time.count() << std::endl;
	auto fullLoop_time = std::chrono::system_clock::now();
	elapsed_seconds = fullLoop_time-DeclHist_time; 
	std::cout << "++++++++++ Full loop time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	//Finding MPVS

	TH1D** projY_xDrift = new TH1D*[maxWireGroup];
	TH1D** projY_tDriftL = new TH1D*[maxWireGroup];
	TH1D** projY_tDriftR = new TH1D*[maxWireGroup];

	TF1** LGfit_xDrift = new TF1*[maxWireGroup];
	TF1** LGfit_tDriftL = new TF1*[maxWireGroup];
	TF1** LGfit_tDriftR = new TF1*[maxWireGroup];

	fitLGParameters fp_xDrift;
	fitLGParameters fp_tDriftL;
	fitLGParameters fp_tDriftR;

	fitExpoParameters expoParams_xDriftL;
	fitExpoParameters expoParams_xDriftR;
	fitExpoParameters expoParams_tDriftL;
	fitExpoParameters expoParams_tDriftR;

	TGraphErrors** MPVplot_xDrift = new TGraphErrors*[maxWireGroup];
	TGraphErrors** MPVplot_tDriftL = new TGraphErrors*[maxWireGroup];
	TGraphErrors** MPVplot_tDriftR = new TGraphErrors*[maxWireGroup];

	TF1 *expoFit_xDriftL = new TF1();
	TF1 *expoFit_xDriftR = new TF1();
	TF1 *expoFit_tDriftL = new TF1();
	TF1 *expoFit_tDriftR = new TF1();

	TCanvas *c_plain = new TCanvas();

	for(int i = 1; i <= maxWireGroup; i++){
		MPVplot_xDrift[i-1] = new TGraphErrors(); 
		MPVplot_tDriftL[i-1] = new TGraphErrors();
		MPVplot_tDriftR[i-1] = new TGraphErrors();
	}

	auto DeclProj_time = std::chrono::system_clock::now();
	elapsed_seconds = DeclProj_time - fullLoop_time; 
	std::cout << "++++++++++ Decl proj time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	for(int i = 1; i <= maxWireGroup; i++){

		expoParams_xDriftL.reset();
		expoParams_xDriftR.reset();
		expoParams_tDriftL.reset();
		expoParams_tDriftR.reset();

		for(int j = 1; j <= N_bins; j++){

			fp_xDrift.reset();
			fp_tDriftL.reset();
			fp_tDriftR.reset();

			projY_xDrift[i-1] = h_dQdx_xDrift[i-1]->ProjectionY(TString::Format("projY_xDrift_%dwires", i), j, j);
			projY_tDriftL[i-1] = h_dQdx_tDriftL[i-1]->ProjectionY(TString::Format("projY_tDriftL_%dwires", i), j, j);
			projY_tDriftR[i-1] = h_dQdx_tDriftR[i-1]->ProjectionY(TString::Format("projY_tDriftR_%dwires", i), j, j);

			SetLGParameters(projY_xDrift[i-1], fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
			SetLGParameters(projY_tDriftL[i-1], fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.lb, fp_tDriftL.ub);
			SetLGParameters(projY_tDriftR[i-1], fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.lb, fp_tDriftR.ub);

			LGfit_xDrift[i-1] = fitter(projY_xDrift[i-1], fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
			LGfit_tDriftL[i-1] = fitter(projY_tDriftL[i-1], fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.cov, "LG");
			LGfit_tDriftR[i-1] = fitter(projY_tDriftR[i-1], fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.cov, "LG");


			if(fp_xDrift.fp[1] >= 0.){

				MPVplot_xDrift[i-1]->SetPoint(MPVplot_xDrift[i-1]->GetN(), h_dQdx_xDrift[i-1]->GetXaxis()->GetBinCenter(j), fp_xDrift.fp[1]);
				MPVplot_xDrift[i-1]->SetPointError(MPVplot_xDrift[i-1]->GetN() - 1, 0., pointError(projY_xDrift[i-1], fp_xDrift));
			}

			if(fp_tDriftL.fp[1] >= 0.){

				MPVplot_tDriftL[i-1]->SetPoint(MPVplot_tDriftL[i-1]->GetN(), h_dQdx_tDriftL[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftL.fp[1]);
				MPVplot_tDriftL[i-1]->SetPointError(MPVplot_tDriftL[i-1]->GetN() - 1, 0., pointError(projY_tDriftL[i-1], fp_tDriftL));
			}

			if(fp_tDriftR.fp[1] >= 0.){

				MPVplot_tDriftR[i-1]->SetPoint(MPVplot_tDriftR[i-1]->GetN(), h_dQdx_tDriftR[i-1]->GetXaxis()->GetBinCenter(j), fp_tDriftR.fp[1]);
				MPVplot_tDriftR[i-1]->SetPointError(MPVplot_tDriftR[i-1]->GetN() - 1, 0., pointError(projY_tDriftR[i-1], fp_tDriftR));
			}

		}

		auto fitMPV_time = std::chrono::system_clock::now();
		elapsed_seconds = fitMPV_time - DeclProj_time; 
		std::cout << "++++++++++ Fit MPV time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

		SetExpoParameters(expoParams_xDriftL.fp);
		SetExpoParameters(expoParams_xDriftR.fp);
		SetExpoParameters(expoParams_tDriftL.fp);
		SetExpoParameters(expoParams_tDriftR.fp);

		expoParams_xDriftL.lb = -160.;
		expoParams_xDriftL.ub = -40.;
		expoParams_xDriftR.lb = 40.;
		expoParams_xDriftR.ub = 160.;
		expoParams_tDriftL.lb = 0.2;
		expoParams_tDriftL.ub = 1.0;
		expoParams_tDriftR.lb = 0.2;
		expoParams_tDriftR.ub = 1.0;

		expoFit_xDriftL = fitter(MPVplot_xDrift[i-1], expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.efp, expoParams_xDriftL.cov, "expoX");
		expoFit_xDriftR = fitter(MPVplot_xDrift[i-1], expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.efp, expoParams_xDriftR.cov, "expoX");
		expoFit_tDriftL = fitter(MPVplot_tDriftL[i-1], expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.efp, expoParams_tDriftL.cov, "expoT");
		expoFit_tDriftR = fitter(MPVplot_tDriftR[i-1], expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.efp, expoParams_tDriftR.cov, "expoT");


		
		/*std::string saveLocName = mydata_cosmics + "NewCodeTest/MPVplots/MPV_xDrift_" + std::to_string(i) + "wires";
		
		openAndClear(c_plain);
		setMarker(MPVplot_xDrift[i-1], kAzure - 3, 21, 0.5);
		MPVplot_xDrift[i-1]->SetTitle("dQ/dx MPV vs x");
		MPVplot_xDrift[i-1]->GetXaxis()->SetTitle("x (cm)");
		MPVplot_xDrift[i-1]->GetYaxis()->SetTitle("dQ/dx MPV (ADC/cm)");
		MPVplot_xDrift[i-1]->Draw("AP");
		expoFit_xDriftL->SetLineColor(kRed);
		expoFit_xDriftL->Draw("SAME");
		expoFit_xDriftR->SetLineColor(kViolet+4);
		expoFit_xDriftR->Draw("SAME");
		TPaveStats *stats = statsBox({.31,.65,.69,.88}, track_count, expoFit_xDriftL, expoFit_xDriftR);
		saveFig(c_plain, mydata_cosmics + fileSaveLoc + "MPVplots/MPV_xDrift_" + fileSaveID);*/

		auto fitExpo_time = std::chrono::system_clock::now();
		elapsed_seconds = fitExpo_time - fitMPV_time; 
		std::cout << "++++++++++ Fit Expo time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

		std::cout << "====================================================================" << std::endl;
		std::cout << "===============================RESULTS==============================" << std::endl;
		std::cout << "Lifetime, x, L: " << expoFit_xDriftL->GetParameter(1) << " \u00b1 " << expoFit_xDriftL->GetParError(1) << std::endl;
		std::cout << "Lifetime, x, R: " << expoFit_xDriftR->GetParameter(1) << " \u00b1 " << expoFit_xDriftR->GetParError(1) << std::endl;
		std::cout << "Lifetime, t, L: " << expoFit_tDriftL->GetParameter(1) << " \u00b1 " << expoFit_tDriftL->GetParError(1) << std::endl;
		std::cout << "Lifetime, t, R: " << expoFit_tDriftR->GetParameter(1) << " \u00b1 " << expoFit_tDriftR->GetParError(1) << std::endl;
		
		N_wires.push_back(i);
		tau_xDrift_LR[0].push_back(expoFit_xDriftL->GetParameter(1));
		errTau_xDrift_LR[0].push_back(expoFit_xDriftL->GetParError(1));
		tau_xDrift_LR[1].push_back(expoFit_xDriftR->GetParameter(1));
		errTau_xDrift_LR[1].push_back(expoFit_xDriftR->GetParError(1));

		tau_tDrift_LR[0].push_back(expoFit_tDriftL->GetParameter(1));
		errTau_tDrift_LR[0].push_back(expoFit_tDriftL->GetParError(1));
		tau_tDrift_LR[1].push_back(expoFit_tDriftR->GetParameter(1));
		errTau_tDrift_LR[1].push_back(expoFit_tDriftR->GetParError(1));

		if(expoParams_xDriftL.fp[1] >= 0. && expoParams_xDriftL.fp[1] <= 50. && expoParams_xDriftL.efp[1] <= 50){
			g_lifeVwiresXL->SetPoint(g_lifeVwiresXL->GetN(), i, expoFit_xDriftL->GetParameter(1));
			g_lifeVwiresXL->SetPointError(g_lifeVwiresXL->GetN()-1, 0., expoFit_xDriftL->GetParError(1));
		}
		if(expoParams_xDriftR.fp[1] >= 0. && expoParams_xDriftR.fp[1] <= 50. && expoParams_xDriftR.efp[1] <= 50){
			g_lifeVwiresXR->SetPoint(g_lifeVwiresXR->GetN(), i, expoFit_xDriftR->GetParameter(1));
			g_lifeVwiresXR->SetPointError(g_lifeVwiresXR->GetN()-1, 0., expoFit_xDriftR->GetParError(1));
		}
		if(expoParams_tDriftL.fp[1] >= 0. && expoParams_tDriftL.fp[1] <= 50. && expoParams_tDriftL.efp[1] <= 50){
			g_lifeVwiresTL->SetPoint(g_lifeVwiresTL->GetN(), i, expoFit_tDriftL->GetParameter(1));
			g_lifeVwiresTL->SetPointError(g_lifeVwiresTL->GetN()-1, 0., expoFit_tDriftL->GetParError(1));
		}
		if(expoParams_tDriftR.fp[1] >= 0. && expoParams_tDriftR.fp[1] <= 50. && expoParams_tDriftR.efp[1] <= 50){
			g_lifeVwiresTR->SetPoint(g_lifeVwiresTR->GetN(), i, expoFit_tDriftR->GetParameter(1));
			g_lifeVwiresTR->SetPointError(g_lifeVwiresTR->GetN()-1, 0., expoFit_tDriftR->GetParError(1));
		}

		auto makePlots_time = std::chrono::system_clock::now();
		elapsed_seconds = makePlots_time - fitExpo_time; 
		std::cout << "++++++++++ Make plots time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	}

	openAndClear(c_plain);
	TMultiGraph *mg = new TMultiGraph();
	setMarker(g_lifeVwiresXL, kRed, 21, 1.0);
	setMarker(g_lifeVwiresXR, kViolet+4, 21, 1.0);
	setMarker(g_lifeVwiresTL, kRed, 47, 1.0);
	setMarker(g_lifeVwiresTR, kViolet+4, 47, 1.0);

	g_lifeVwiresXL->SetTitle("x, left TPC");
	g_lifeVwiresXR->SetTitle("x, right TPC");
	g_lifeVwiresTL->SetTitle("t, left TPC");
	g_lifeVwiresTR->SetTitle("t, right TPC");

	mg->SetTitle("Lifetime vs number of combined wires;N; #tau (ms)");
	mg->Add(g_lifeVwiresXL);
	mg->Add(g_lifeVwiresXR);
	mg->Add(g_lifeVwiresTL);
	mg->Add(g_lifeVwiresTR);
	mg->Draw("AP");
	c_plain->BuildLegend(0.4,0.7,0.6,0.88);
	saveFig(c_plain, mydata_cosmics + fileSaveLoc + "MultiWirePlots/LifetimeVwireNum" + fileSaveID);
	
	return 0;
}

//The Landau-Gaussian function
double_t langaufun(Double_t *x, Double_t *par) {

	//Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
	//Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
	//Markus Friedl (Markus.Friedl@cern.ch) "langaus.C"

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

double_t expofunX(Double_t *x, Double_t *par){

	double vDrift = 156.267;
	Double_t xx = x[0];
	Double_t f = par[0]*exp(-((200-TMath::Abs(xx))/(par[1]*vDrift)));
   	return f;
}

double_t expofunT(Double_t *x, Double_t *par){

	Double_t xx = x[0];
	Double_t f = par[0]*exp(-(xx/par[1]));
   	return f;
}

//Function to set LG fit initial values
void SetLGParameters(TH1D *h, Double_t *fp, Double_t *efp, Double_t &lb, Double_t &ub){
	//Adapted from Lan Nguyen's code which was adapted from Dom Barker's code etc.

    // Set fit parameters and errors  
    double nentries = h->GetEntries();
    double max = h->GetBinCenter(h->GetMaximumBin()); //x position of max dQdx value 
    double area = h->GetEntries()*h->GetBinWidth(1)*0.8 ; //normalisation constant //starting value for the fit - fit doesn't include whole hist
    double rms = h->GetRMS();
	
	//Starting guesses for parameters
    fp[0] = max / 100 * 4; //scale parameter of Landau density, starting guess about 5% of MPV x value
    fp[1] = max; //MPV: max dQdx value 		(position)
    fp[2] = area; // total area or integral //NORM (bin width) getNentries*Getbinwidth
    fp[3] = rms / 10;  //GSigma typically smaller than distribution rms
    
	//Starting guesses for errors
    efp[0] = max / 100 * 0.5; //0.5% of max dQdx value  	
    efp[1] = max / 100 * 8;	// 8% of max dQdx value
    efp[2] = area * 0.1;	
    efp[3] = rms * 0.01;
  
	//Lower and upper bounds for fit: sufficient to find the peak of the distribution precisely,
	//whilst staying in high stats bins to reduce statistical fluctuations
    double dQdxpeak = h->GetMaximumBin(); 
    lb = h->GetBinCenter(dQdxpeak-8);	
	ub = h->GetBinCenter(dQdxpeak+12);

}

void SetExpoParameters(Double_t *fp){
	fp[0] = 1000.;
	fp[1] = 10.;

}

TF1 *fitter(TH1D *h, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, std::string funcName){

	Double_t(*func)(Double_t *,Double_t *);
	int func_index;
	int nParams;

	if(!strcmp(funcName.c_str(), "expoX")){
		func = expofunX;
		func_index = 1;
		nParams = 2;
	}
	else if(!strcmp(funcName.c_str(), "expoT")){
		func = expofunT;
		func_index = 2;
		nParams = 2;
	}
	else if(!strcmp(funcName.c_str(), "LG")){
		func = langaufun;
		func_index = 3;
		nParams = 4;
	}
	else{
		std::cout << "UNKNOWN FUNCTION" << std::endl;
	}

	Char_t FitFuncName[100]; 
  	sprintf(FitFuncName,"Fitfcn_%s",h->GetName());

	TF1 *fitfunc = new TF1(FitFuncName,func,lbound,ubound, nParams);
	
	if(func_index == 1 || func_index == 2){
		fitfunc->SetParameters(fitparams[0], fitparams[1]);
		fitfunc->SetParError(0,fiterrors[0]);
		fitfunc->SetParError(1,fiterrors[1]);
		fitfunc->SetParLimits(1,0.,25.);
		fitfunc->SetParNames("Norm","Lifetime");
	}
	else if(func_index == 3){
		fitfunc->SetParameters(fitparams[0], fitparams[1], fitparams[2], fitparams[3]);
		fitfunc->SetParError(0,fiterrors[0]);
		if (fiterrors[0]==0) fitfunc->FixParameter(0,fitparams[0]); //if scale parameter error is 0 scale parameter is fixed
		fitfunc->SetParError(1,fiterrors[1]);
		fitfunc->SetParError(2,fiterrors[2]);
		fitfunc->SetParError(3,fiterrors[3]);
		fitfunc->SetParLimits(0,20,60);
		fitfunc->SetParLimits(3,10,200);
		fitfunc->SetParNames("Width","MPV","TotalArea","GSigma"); 
	}

	h->Fit(FitFuncName,"LREQ");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode
  	//Other fitting options https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html (7.1.1)
  
	TString fitOutcome = gMinuit->fCstatu.Data();

	for(int i = 0; i < nParams; i++){
		fitparams[i] = h->GetFunction(FitFuncName)->GetParameter(i);	
		fiterrors[i] = h->GetFunction(FitFuncName)->GetParError(i);
	}

	if (!fitOutcome.BeginsWith("SUCC")) { 
		for(int i = 0; i < nParams; i++){
			fitparams[i] = -1000.0;	
			fiterrors[i] = -1000.0;
			covmat[i] = -1000.0;
		} 
	}

	return(fitfunc);
	
}

TF1 *fitter(TGraphErrors *g, Double_t lbound, Double_t ubound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, std::string funcName){

	Double_t(*func)(Double_t *,Double_t *);
	int func_index;
	int nParams;

	if(!strcmp(funcName.c_str(), "expoX")){
		func = expofunX;
		func_index = 1;
		nParams = 2;
	}
	else if(!strcmp(funcName.c_str(), "expoT")){
		func = expofunT;
		func_index = 2;
		nParams = 2;
	}
	else if(!strcmp(funcName.c_str(), "LG")){
		func = langaufun;
		func_index = 3;
		nParams = 4;
	}
	else{
		std::cout << "UNKNOWN FUNCTION" << std::endl;
	}

	Char_t FitFuncName[100]; 
  	sprintf(FitFuncName,"Fitfcn_%s",g->GetName());

	TF1 *fitfunc = new TF1(FitFuncName,func,lbound,ubound, nParams);
	
	if(func_index == 1 || func_index == 2){
		fitfunc->SetParameters(fitparams[0], fitparams[1]);
		fitfunc->SetParError(0,fiterrors[0]);
		fitfunc->SetParError(1,fiterrors[1]);
		fitfunc->SetParNames("Norm","Lifetime");
	}
	else if(func_index == 3){
		fitfunc->SetParameters(fitparams[0], fitparams[1], fitparams[2], fitparams[3]);
		fitfunc->SetParError(0,fiterrors[0]);
		if (fiterrors[0]==0) fitfunc->FixParameter(0,fitparams[0]); //if scale parameter error is 0 scale parameter is fixed
		fitfunc->SetParError(1,fiterrors[1]);
		fitfunc->SetParError(2,fiterrors[2]);
		fitfunc->SetParError(3,fiterrors[3]);
		fitfunc->SetParLimits(0,20,60);
		fitfunc->SetParLimits(3,10,200);
		fitfunc->SetParNames("Width","MPV","TotalArea","GSigma"); 
	}

	g->Fit(FitFuncName,"LREQ");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode
  	//Other fitting options https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html (7.1.1)
  
	TString fitOutcome = gMinuit->fCstatu.Data();

	for(int i = 0; i < nParams; i++){
		fitparams[i] = g->GetFunction(FitFuncName)->GetParameter(i);	
		fiterrors[i] = g->GetFunction(FitFuncName)->GetParError(i);
	}

	if (!fitOutcome.BeginsWith("SUCC")) { 
		for(int i = 0; i < nParams; i++){
			fitparams[i] = -1000.0;	
			fiterrors[i] = -1000.0;
			covmat[i] = -1000.0;
		} 
	}

	return(fitfunc);
	
}

double pointError(TH1D *proj_y, fitLGParameters fitParams){

	int Nslice = proj_y->GetEntries();
	return sqrt(pow(fitParams.efp[1],2)+pow((fitParams.fp[1]/sqrt(Nslice)),2));

}

void setMarker(TH1D *h, Color_t color, Style_t style, Size_t size){
	h->SetMarkerColor(color);
	h->SetMarkerStyle(style);
	h->SetMarkerSize(size);
}

void setMarker(TGraphErrors *g, Color_t color, Style_t style, Size_t size){
	g->SetMarkerColor(color);
	g->SetMarkerStyle(style);
	g->SetMarkerSize(size);
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

void addPointAndError(TGraphErrors *g, double x, double y, double errx, double erry){
	g->SetPoint(g->GetN(), x, y);
	g->SetPointError(g->GetN()-1, errx, erry);
}

void addPointAndError(TGraphErrors *g, int x, double y, double errx, double erry){
	g->SetPoint(g->GetN(), x, y);
	g->SetPointError(g->GetN()-1, errx, erry);
}

TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TF1* f2){
	
	TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

	pt->SetBorderSize(1);
	pt->SetFillColor(0);
	pt->AddText(Form("AC cosmics: %d",track_count));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
	pt->AddText(Form("#chi^{2} / DoF: %g / %g",f1->GetChisquare(),f1->GetNDF()));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
	pt->AddText(Form("Lifetime: %g#pm%g",f1->GetParameter(1),f1->GetParError(1)));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(f1->GetLineColor());
	pt->AddText(Form("#chi^{2} / DoF: %g / %g",f2->GetChisquare(),f2->GetNDF()));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
	pt->AddText(Form("Lifetime: %g#pm%g",f2->GetParameter(1),f2->GetParError(1)));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(f2->GetLineColor());
	
	return pt;
}

TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1){
	TPaveText *pt = new TPaveText(pos[0], pos[1], pos[2], pos[3], "blNDC");

	pt->SetBorderSize(1);
	pt->SetFillColor(0);
	pt->AddText(Form("AC cosmics: %d",track_count));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
	pt->AddText(Form("#chi^{2} / DoF: %g / %g",f1->GetChisquare(),f1->GetNDF()));
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

	/*for(int j=1; j<=20; j++){
		//Histograms
		auto loopStart_time = std::chrono::system_clock::now();

		TH2D *h_dQdx_xDrift = new TH2D("h_dQdx_xDrift","dQ/dx vs x", 100, -200, 200, 75, 200, 1800);
		TGraph *g_XvT = new TGraph;
		TH2D *h_dQdx_tDriftL = new TH2D("h_dQdx_tDriftL","dQ/dx vs t Left TPC", 100, 0, 1.3, 75, 200, 1800);
		TH2D *h_dQdx_tDriftR = new TH2D("h_dQdx_tDriftR","dQ/dx vs t Right TPC", 100, 0, 1.3, 75, 200, 1800);
				
		int track_count = 0;
		int overall_count = 0;
		int N = j; //no. of wires to group together

		//vector for grouping wires together
		std::vector<double> dQdx;
		std::vector<double> wire;
		std::vector<double> x;
		std::vector<double> t;
		std::vector<double> dQdx_wireSum;
		std::vector<int> count;
		std::vector<double> x_wireSum;
		std::vector<double> t_wireSum;

		auto loopPreAmble_time = std::chrono::system_clock::now();
		elapsed_seconds = loopPreAmble_time - loopStart_time; 
    	std::cout << "++++++++++ Loop preamble time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

		treereader.Restart();
		while(treereader.Next()){
			
			//if(*read_xi < -200. || *read_xi > 200. || *read_xf < -200. || *read_xf > 200.){
			//if(true){
			if(*read_selected == 1){
				auto startFillingVectors_time = std::chrono::system_clock::now();
				//Remove all nan entries
				for (int i = 0; i < read_dqdx.GetSize(); i++){
					if(std::isnan(read_x[i]) == 0){ //get rid of nan entries
						dQdx.push_back(read_dqdx[i]);
						x.push_back(read_x[i]);
						t.push_back(read_T[i]);
						wire.push_back(read_wire[i]);
					}
				}

				auto removeNaNentries_time = std::chrono::system_clock::now();
				elapsed_seconds = removeNaNentries_time - startFillingVectors_time; 
    			std::cout << "++++++++++ Remove NaN entries time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

				//now have dQdx, x, t and wire for all elements that didn't have x as nan. 


				auto maxWire = std::max_element(wire.begin(), wire.end());
				auto minWire = std::min_element(wire.begin(), wire.end());
				int quotient = (*maxWire - *minWire) / N;
				dQdx_wireSum.resize(quotient + 1, 0.);
				x_wireSum.resize(quotient + 1, 0.);
				t_wireSum.resize(quotient + 1, 0.);
				count.resize(quotient + 1, 0);

				for (int i = 0; i < dQdx.size(); i++){
					dQdx_wireSum[ (int)((wire[i] - *minWire)/N) ] += dQdx[i];
					count[ (int)((wire[i] - *minWire)/N) ] += 1;
					x_wireSum[ (int)((wire[i] - *minWire)/N) ] += x[i];
					t_wireSum[ (int)((wire[i] - *minWire)/N) ] += t[i];
				}

				auto fillQuotientVector_time = std::chrono::system_clock::now();
				elapsed_seconds = fillQuotientVector_time - removeNaNentries_time; 
    			std::cout << "++++++++++ Fill quotient vector time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;
				

				for(int i=0; i < dQdx_wireSum.size(); i++){
					if(count[i] > 0){
						h_dQdx_xDrift->Fill(x_wireSum[i]/count[i], dQdx_wireSum[i]/count[i]);
						g_XvT->SetPoint(g_XvT->GetN(), x_wireSum[i]/count[i], t_wireSum[i]/(count[i]*2000) - 0.2 - *read_t0/1000000);
						
						//Left TPC
						if(-200. < x_wireSum[i]/count[i] && x_wireSum[i]/count[i] < 0.){
							h_dQdx_tDriftL->Fill(t_wireSum[i]/(count[i]*2000) - 0.2 - *read_t0/1000000, dQdx_wireSum[i]/count[i]);
						}

						//Right TPC
						if(0. < x_wireSum[i]/count[i] && x_wireSum[i]/count[i] < 200.){
							h_dQdx_tDriftR->Fill(t_wireSum[i]/(count[i]*2000) - 0.2 - *read_t0/1000000, dQdx_wireSum[i]/count[i]);
						}
					}
				}

				auto fillHist_time = std::chrono::system_clock::now();
				elapsed_seconds = fillHist_time - fillQuotientVector_time; 
    			std::cout << "++++++++++ Fill hist time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

				//clear all vectors ready for next loop
				dQdx_wireSum.clear();
				wire.clear();
				x.clear();
				dQdx.clear();
				count.clear();
				x_wireSum.clear();
				t.clear();
				t_wireSum.clear();
				track_count += 1;

				auto clearAndReset_time = std::chrono::system_clock::now();
				elapsed_seconds = clearAndReset_time - fillHist_time; 
				std::cout << "++++++++++ Clear & reset time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;
			}

			
			
			overall_count += 1;
		}

		auto TTreeReader_time = std::chrono::system_clock::now();
		elapsed_seconds = TTreeReader_time-loopPreAmble_time; 
    	std::cout << "++++++++++ TTreeReader overall time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;
		
		std::cout << "AC track count is " << track_count << std::endl;
		std::cout << "Overall track count is " << overall_count << std::endl;

		if(plot_dQdx){

			std::cout << "----Plotting dQ/dx-----------------------------------------------------" << std::endl;

			TCanvas *c_dQdx_xDrift = new TCanvas();
			h_dQdx_xDrift->SetStats(0);
			h_dQdx_xDrift->GetXaxis()->SetTitle("x [cm]");
			h_dQdx_xDrift->GetYaxis()->SetTitle("dQ/dx [ADc/cm]");

			std::stringstream ssX;
			std::stringstream ssY;
			ssX << std::setprecision(3) << h_dQdx_xDrift->GetXaxis()->GetBinWidth(1);
			ssY << std::setprecision(3) << h_dQdx_xDrift->GetYaxis()->GetBinWidth(1);

			std::string colorLabel = "Entries / (" + ssX.str() + " " + "cm" + " #times " + ssY.str() + " " + "ADC/cm" + ")";
			TLatex *colorTitle = new TLatex(0.95,0.5,colorLabel.c_str());
			colorTitle->SetNDC(1);
			colorTitle->SetTextFont(42);
			colorTitle->SetTextAngle(90);
			colorTitle->SetTextSize(0.035);
			colorTitle->SetTextAlign(22);
			colorTitle->Draw();
			
			c_dQdx_xDrift->SetRightMargin(0.18);
			c_dQdx_xDrift->cd();
			h_dQdx_xDrift->Draw("COLZ");
			TPaveText *stats_dQdx = new TPaveText(.55,.80,.80,.88,"blNDC");

			stats_dQdx->SetBorderSize(1);
			stats_dQdx->SetFillColor(0);
			stats_dQdx->AddText(Form("AC cosmics: %d",track_count));
			((TText*)stats_dQdx->GetListOfLines()->Last())->SetTextColor(1);
			stats_dQdx->Draw();

			c_dQdx_xDrift->SaveAs((mydata_cosmics + std::to_string(N) + "_dQdx_multiwireTEST.png").c_str());

			TCanvas *c_XvT = new TCanvas();		
			g_XvT->SetMarkerColor(kViolet + 3);
			g_XvT->SetMarkerStyle(kPlus);
			g_XvT->SetMarkerSize(0.1);
			g_XvT->GetXaxis()->SetTitle("x (cm)");
			g_XvT->GetYaxis()->SetTitle("t (ms)");
			g_XvT->SetTitle("drift time vs drift distance");
			g_XvT->Draw("AP");
			TPaveText *stats_XvT = new TPaveText(.63,.80,.88,.88,"blNDC");

			stats_XvT->SetBorderSize(1);
			stats_XvT->SetFillColor(0);
			stats_XvT->AddText(Form("AC cosmics: %d",track_count));
			((TText*)stats_XvT->GetListOfLines()->Last())->SetTextColor(1);
			stats_XvT->Draw();

			c_XvT->SaveAs((mydata_cosmics + std::to_string(N) + "_XvTmultiwireTEST.png").c_str());

			TCanvas *c_dQdx_tDrift = new TCanvas();
			c_dQdx_tDrift->SetWindowSize(1500,500);
			TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
			TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
			c_dQdx_tDrift->cd();
			pad1->SetRightMargin(0.18);
			pad1->SetLeftMargin(0.15);
			pad1->Draw();
			pad1->cd();
			h_dQdx_tDriftL->SetStats(0);
			h_dQdx_tDriftL->GetXaxis()->SetTitle("t [ms]");
			h_dQdx_tDriftL->GetYaxis()->SetTitle("dQ/dx [ADc/cm]");
			h_dQdx_tDriftL->Draw("COLZ");

			std::stringstream ssX1;
			std::stringstream ssY1;
			ssX1 << std::setprecision(3) << h_dQdx_tDriftL->GetXaxis()->GetBinWidth(1);
			ssY1 << std::setprecision(3) << h_dQdx_tDriftL->GetYaxis()->GetBinWidth(1);

			std::string colorLabel1 = "Entries / (" + ssX.str() + " " + "ms" + " #times " + ssY.str() + " " + "ADC/cm" + ")";
			TLatex *colorTitle1 = new TLatex(0.95,0.5,colorLabel1.c_str());
			colorTitle1->SetNDC(1);
			colorTitle1->SetTextFont(42);
			colorTitle1->SetTextAngle(90);
			colorTitle1->SetTextSize(0.035);
			colorTitle1->SetTextAlign(22);
			colorTitle1->Draw();
			
			TPaveText *stats_dQdx_L = new TPaveText(.55,.80,.80,.88,"blNDC");

			stats_dQdx_L->SetBorderSize(1);
			stats_dQdx_L->SetFillColor(0);
			stats_dQdx_L->AddText(Form("AC cosmics: %d",track_count));
			((TText*)stats_dQdx_L->GetListOfLines()->Last())->SetTextColor(1);
			stats_dQdx_L->Draw();

			c_dQdx_tDrift->cd();
			pad2->SetRightMargin(0.18);
			pad2->SetLeftMargin(0.15);
			pad2->Draw();
			pad2->cd();
			h_dQdx_tDriftR->SetStats(0);
			h_dQdx_tDriftR->GetXaxis()->SetTitle("t [ms]");
			h_dQdx_tDriftR->GetYaxis()->SetTitle("dQ/dx [ADc/cm]");
			h_dQdx_tDriftR->Draw("COLZ");

			std::stringstream ssX2;
			std::stringstream ssY2;
			ssX2 << std::setprecision(3) << h_dQdx_tDriftR->GetXaxis()->GetBinWidth(1);
			ssY2 << std::setprecision(3) << h_dQdx_tDriftR->GetYaxis()->GetBinWidth(1);

			std::string colorLabel2 = "Entries / (" + ssX.str() + " " + "ms" + " #times " + ssY.str() + " " + "ADC/cm" + ")";
			TLatex *colorTitle2 = new TLatex(0.95,0.5,colorLabel2.c_str());
			colorTitle2->SetNDC(1);
			colorTitle2->SetTextFont(42);
			colorTitle2->SetTextAngle(90);
			colorTitle2->SetTextSize(0.035);
			colorTitle2->SetTextAlign(22);
			colorTitle2->Draw();

			TPaveText *stats_dQdx_R = new TPaveText(.55,.80,.80,.88,"blNDC");

			stats_dQdx_R->SetBorderSize(1);
			stats_dQdx_R->SetFillColor(0);
			stats_dQdx_R->AddText(Form("AC cosmics: %d",track_count));
			((TText*)stats_dQdx_R->GetListOfLines()->Last())->SetTextColor(1);
			stats_dQdx_R->Draw();

			c_dQdx_tDrift->SaveAs((mydata_cosmics + std::to_string(N) + "_dQdxVt_multiwireTEST.png").c_str());

		}

		auto plotdQdx_time = std::chrono::system_clock::now();
		elapsed_seconds = plotdQdx_time-TTreeReader_time; 
    	std::cout << "++++++++++ Plot dQdx time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;


		std::cout << "----Finding MPV--------------------------------------------------------" << std::endl;
		
		TGraphErrors *MPVplot = new TGraphErrors(); 

		//Fit parameters
		double fp[4]; //fit parameters
		double efp[4]; //fit parameter errors
		double cov[4]; //covariance matrix
		double lb; //fit lower bound
		double ub; //fit upper bound

		//loop over each bin. For overflow bins change to i=0 to num_xbins + 2.
		for(int i=1; i<h_dQdx_xDrift->GetNbinsX()+1; i++){
			
			//define filenames for projection plots
			std::string bin_number = "Bin" + std::to_string(i);		
			std::string saveLocTemp = mydata_cosmics + "ProjYplots/" + "xDrift_" + dQdxProjYPlotName + bin_number +".png";
			
			//perform the Landau-Gaussian fit
			TH1D *proj_y = h_dQdx_xDrift->ProjectionY("proj_y", i, i);
			SetLGParameters(proj_y, fp, efp, lb, ub);
			TF1 *fitres = LGfitter(proj_y, lb, ub, fp, efp, cov, 0);

			//fill vectors with fit values
			int Nslice = proj_y->GetEntries(); //for stat error
		
			if(fp[1] >= 0){ //Get rid of any points where the fit failed, in which case fp[1] is set to -1000
				MPVplot->SetPoint(MPVplot->GetN(), h_dQdx_xDrift->GetXaxis()->GetBinCenter(i), fp[1]);
				MPVplot->SetPointError(MPVplot->GetN()-1, 0., sqrt(pow(efp[1],2)+pow((fp[1]/sqrt(Nslice)),2)));
			}

			if(plot_projY){

				std::cout << "----Plotting projY-----------------------------------------------------" << std::endl;

				std::string plotTitleTemp;
				plotTitleTemp = "dQ/dx for " + std::to_string(h_dQdx_xDrift->GetXaxis()->GetBinLowEdge(i))
				+ " cm < x < " + std::to_string(h_dQdx_xDrift->GetXaxis()->GetBinLowEdge(i) + h_dQdx_xDrift->GetXaxis()->GetBinWidth(i)) + " cm";

				TCanvas *c9 = new TCanvas();		
				proj_y->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
				std::stringstream ss;
				ss << std::setprecision(3) << proj_y->GetXaxis()->GetBinWidth(1);
				std::string yLabel;
				yLabel = "Entries / (" + ss.str() + " ADC/cm)";
				proj_y->GetYaxis()->SetTitle(yLabel.c_str());
			
				proj_y->SetTitle(plotTitleTemp.c_str()); 
				proj_y->Draw();
				fitres->Draw("same");
				gStyle->SetOptFit(1111); //stats box
				c9->SaveAs(saveLocTemp.c_str());
			}
			
		}

		//Plot and fit MPV
		TCanvas *c_MPV = new TCanvas();
		c_MPV->cd();
        MPVplot->SetMarkerColor(kAzure - 3);
        MPVplot->SetMarkerStyle(21);
        MPVplot->SetMarkerSize(0.5);
        MPVplot->GetXaxis()->SetTitle("x (cm)");
        MPVplot->SetTitle("dQ/dx MPV vs x");
        MPVplot->GetYaxis()->SetTitle("dQ/dx MPV (ADC/cm)");
        MPVplot->Draw("AP");		

		//Fit exponential in both TPCs
		
		double vDrift = 156.267;
		gROOT->GetFunction("expo");
		std::string Lformula = "[0]*exp(-((200+x)/([1]*" + std::to_string(vDrift) + ")))";
		std::string Rformula = "[0]*exp(-((200-x)/([1]*" + std::to_string(vDrift) + ")))";
		TF1 *LeftTPCFit = new TF1("LeftTPCFit", Lformula.c_str(), -160., -40.);
		TF1 *RightTPCFit = new TF1("RightTPCFit", Rformula.c_str(), 40., 160.);
		LeftTPCFit->SetParNames("Norm","Lifetime");
		RightTPCFit->SetParNames("Norm","Lifetime");

		LeftTPCFit->SetLineColor(kRed);
		RightTPCFit->SetLineColor(kGreen + 1);

		//Set initial parameters based on roughly the expected values (lifetime simulated as 10ms)
		LeftTPCFit->SetParameter(0,1000);
		LeftTPCFit->SetParameter(1,10);
		RightTPCFit->SetParameter(0,1000);
		RightTPCFit->SetParameter(1,10);
		
		//Fit function
		std::cout << "----Fitting exponentials-----------------------------------------------" << std::endl;
		MPVplot->Fit(LeftTPCFit, "R");
		MPVplot->Fit(RightTPCFit, "R+");

		//Fit results
		double Lpar0 = LeftTPCFit->GetParameter(0);
		double Lpar1 = LeftTPCFit->GetParameter(1);
		double Lerr0 = LeftTPCFit->GetParError(0);
		double Lerr1 = LeftTPCFit->GetParError(1);
		double Rpar0 = RightTPCFit->GetParameter(0);
		double Rpar1 = RightTPCFit->GetParameter(1);
		double Rerr0 = RightTPCFit->GetParError(0);
		double Rerr1 = RightTPCFit->GetParError(1);

		//Making stats box
		MPVplot->SetStats(0);
		
		TPaveText *pt = new TPaveText(.31,.65,.69,.88,"blNDC");
		pt->SetBorderSize(1);
		pt->SetFillColor(0);
		pt->AddText(Form("AC cosmics: %d",track_count));
		((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
		pt->AddText(Form("#chi^{2} / DoF: %g / %g",LeftTPCFit->GetChisquare(),LeftTPCFit->GetNDF()));
		((TText*)pt->GetListOfLines()->Last())->SetTextColor(LeftTPCFit->GetLineColor());
		pt->AddText(Form("Lifetime: %g#pm%g",Lpar1,Lerr1));
		((TText*)pt->GetListOfLines()->Last())->SetTextColor(LeftTPCFit->GetLineColor());
		pt->AddText(Form("#chi^{2} / DoF: %g / %g",RightTPCFit->GetChisquare(),RightTPCFit->GetNDF()));
		((TText*)pt->GetListOfLines()->Last())->SetTextColor(RightTPCFit->GetLineColor());
		pt->AddText(Form("Lifetime: %g#pm%g",Rpar1,Rerr1));
		((TText*)pt->GetListOfLines()->Last())->SetTextColor(RightTPCFit->GetLineColor());

		//Draw stats box and save plot
		pt->Draw();
		
		//Print fit parameters
		std::cout << "Fit L TPC : par0 = " << Lpar0 << " ; par1 = " << Lpar1 << " ; err0 = " << Lerr0 << " ; err1 = " << Lerr1 << std::endl;
		std::cout << "Fit R TPC : par0 = " << Rpar0 << " ; par1 = " << Rpar1 << " ; err0 = " << Rerr0 << " ; err1 = " << Rerr1 << std::endl;

		std::cout << "Lifetime Left TPC : " << Lpar1<< " \u00b1 " << Lerr1 << "Norm Left TPC : " << Lpar0 << " \u00b1 " << Lerr0 << std::endl;
		std::cout << "Lifetime Right TPC : " << Rpar1 <<  " \u00b1 " << Rerr1 << "Norm Right TPC : " << Rpar0 << " \u00b1 " << Rerr0 << std::endl;

		std::vector<double> fitParams_xDrift = {Lpar1, Lerr1, Rpar1, Rerr1};
	
		c_MPV->SaveAs((mydata_cosmics + std::to_string(N) + "_MPV_multiwireTEST.png").c_str());

		//Histogram
		TGraphErrors *MPVplot_L = new TGraphErrors(); 

		//Fit parameters
		double fp[4]; //fit parameters
		double efp[4]; //fit parameter errors
		double cov[4]; //covariance matrix
		double lb; //fit lower bound
		double ub; //fit upper bound

		//loop over each bin. For overflow bins change to i=0 to num_xbins + 2.
		for(int i=1; i<h_dQdx_tDriftL->GetNbinsX()+1; i++){
			
			//define filenames for projection plots
			std::string bin_number = "Bin" + std::to_string(i);		
			std::string saveLocTemp = mydata_cosmics + "ProjYplots/" + "tDriftL" + dQdxProjYPlotName + bin_number +".png";
			
			//perform the Landau-Gaussian fit
			TH1D *proj_y = h_dQdx_tDriftL->ProjectionY("proj_y", i, i);
			SetLGParameters(proj_y, fp, efp, lb, ub);
			TF1 *fitres = LGfitter(proj_y, lb, ub, fp, efp, cov, 0);

			//fill vectors with fit values
			int Nslice = proj_y->GetEntries(); //for stat error
		
			if(fp[1] >= 0){ //Get rid of any points where the fit failed, in which case fp[1] is set to -1000
				MPVplot_L->SetPoint(MPVplot_L->GetN(), h_dQdx_tDriftL->GetXaxis()->GetBinCenter(i), fp[1]);
				MPVplot_L->SetPointError(MPVplot_L->GetN()-1, 0., sqrt(pow(efp[1],2)+pow((fp[1]/sqrt(Nslice)),2)));
			}

			if(plot_projY){

				std::cout << "----Plotting projY-----------------------------------------------------" << std::endl;

				std::string plotTitleTemp;
				
				plotTitleTemp = "dQ/dx for " + std::to_string(h_dQdx_tDriftL->GetXaxis()->GetBinLowEdge(i))
				+ " ms < t < " + std::to_string(h_dQdx_tDriftL->GetXaxis()->GetBinLowEdge(i) + h_dQdx_tDriftL->GetXaxis()->GetBinWidth(i)) + " ms";

				TCanvas *c9 = new TCanvas();		
				proj_y->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
				std::stringstream ss;
				ss << std::setprecision(3) << proj_y->GetXaxis()->GetBinWidth(1);
				std::string yLabel;
				yLabel = "Entries / (" + ss.str() + " ADC/cm)";
				proj_y->GetYaxis()->SetTitle(yLabel.c_str());
			
				proj_y->SetTitle(plotTitleTemp.c_str()); 
				proj_y->Draw();
				fitres->Draw("same");
				gStyle->SetOptFit(1111); //stats box
				c9->SaveAs(saveLocTemp.c_str());
			}
			
		}

		TGraphErrors *MPVplot_R = new TGraphErrors(); 

		//Fit parameters
		double fp[4]; //fit parameters
		double efp[4]; //fit parameter errors
		double cov[4]; //covariance matrix
		double lb; //fit lower bound
		double ub; //fit upper bound

		//loop over each bin. For overflow bins change to i=0 to num_xbins + 2.
		for(int i=1; i<h_dQdx_tDriftR->GetNbinsX()+1; i++){
			
			//define filenames for projection plots
			std::string bin_number = "Bin" + std::to_string(i);		
			std::string saveLocTemp = mydata_cosmics + "ProjYplots/" + "tDriftR" + dQdxProjYPlotName + bin_number +".png";
			
			//perform the Landau-Gaussian fit
			TH1D *proj_y = h_dQdx_tDriftR->ProjectionY("proj_y", i, i);
			SetLGParameters(proj_y, fp, efp, lb, ub);
			TF1 *fitres = LGfitter(proj_y, lb, ub, fp, efp, cov, 0);

			//fill vectors with fit values
			int Nslice = proj_y->GetEntries(); //for stat error
		
			if(fp[1] >= 0){ //Get rid of any points where the fit failed, in which case fp[1] is set to -1000
				MPVplot_R->SetPoint(MPVplot_R->GetN(), h_dQdx_tDriftR->GetXaxis()->GetBinCenter(i), fp[1]);
				MPVplot_R->SetPointError(MPVplot_R->GetN()-1, 0., sqrt(pow(efp[1],2)+pow((fp[1]/sqrt(Nslice)),2)));
			}

			if(plot_projY){

				std::cout << "----Plotting projY-----------------------------------------------------" << std::endl;

				std::string plotTitleTemp;
		
				plotTitleTemp = "dQ/dx for " + std::to_string(h_dQdx_tDriftR->GetXaxis()->GetBinLowEdge(i))
				+ " ms < t < " + std::to_string(h_dQdx_tDriftR->GetXaxis()->GetBinLowEdge(i) + h_dQdx_tDriftR->GetXaxis()->GetBinWidth(i)) + " ms";
				

				TCanvas *c9 = new TCanvas();		
				proj_y->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
				std::stringstream ss;
				ss << std::setprecision(3) << proj_y->GetXaxis()->GetBinWidth(1);
				std::string yLabel;
				yLabel = "Entries / (" + ss.str() + " ADC/cm)";
				proj_y->GetYaxis()->SetTitle(yLabel.c_str());
			
				proj_y->SetTitle(plotTitleTemp.c_str()); 
				proj_y->Draw();
				fitres->Draw("same");
				gStyle->SetOptFit(1111); //stats box
				c9->SaveAs(saveLocTemp.c_str());
			}
			
		}


		TCanvas *c_tMPV = new TCanvas();
		c_tMPV->SetWindowSize(1500,500);
		TPad *pad1_tMPV = new TPad("pad1_tMPV", "", 0, 0, 0.5, 1.0);
		TPad *pad2_tMPV = new TPad("pad2_tMPV", "", 0.5, 0, 1.0, 1.0);
		c_tMPV->cd();
		pad1_tMPV->SetRightMargin(0.18);
		pad1_tMPV->SetLeftMargin(0.15);
		pad1_tMPV->Draw();
		pad1_tMPV->cd();
        MPVplot_L->SetMarkerColor(kAzure - 3);
        MPVplot_L->SetMarkerStyle(21);
        MPVplot_L->SetMarkerSize(0.5);
        MPVplot_L->GetXaxis()->SetTitle("t (ms)");
        MPVplot_L->SetTitle("dQ/dx MPV vs t Left TPC");
        MPVplot_L->GetYaxis()->SetTitle("dQ/dx MPV (ADC/cm)");
        MPVplot_L->Draw("AP");	

		//Function:
		std::string t_formula = "[0]*exp(-(x/[1]))";
		TF1 *TPCFit = new TF1("TPCFit", t_formula.c_str(), 0.2, 1.0);
		TPaveText *pt_L = new TPaveText(.50,.75,.80,.88,"blNDC");

		gROOT->GetFunction("expo");
		TPCFit->SetParNames("Norm","Lifetime");
		TPCFit->SetLineColor(kRed);

		//Set initial parameters based on roughly the expected values (lifetime simulated as 10ms)
		TPCFit->SetParameter(0,1000);
		TPCFit->SetParameter(1,10);

		//Fit function
		std::cout << "----Fitting exponential in time---------------------------------------" << std::endl;
		MPVplot_L->Fit(TPCFit, "R");

		//Making stats box
		MPVplot_L->SetStats(0);

		std::cout << TPCFit->GetChisquare() << std::endl;
		std::cout << TPCFit->GetNDF() << std::endl;
		//Stats box position and text
		pt_L->SetBorderSize(1);
		pt_L->SetFillColor(0);
		pt_L->AddText(Form("AC cosmics: %d",track_count));
		((TText*)pt_L->GetListOfLines()->Last())->SetTextColor(1);
		pt_L->AddText(Form("#chi^{2} / DoF: %g / %d",TPCFit->GetChisquare(),TPCFit->GetNDF()));
		((TText*)pt_L->GetListOfLines()->Last())->SetTextColor(TPCFit->GetLineColor());
		pt_L->AddText(Form("Lifetime: %g#pm%g",TPCFit->GetParameter(1),TPCFit->GetParError(1)));
		((TText*)pt_L->GetListOfLines()->Last())->SetTextColor(TPCFit->GetLineColor());
		
		//Draw stats box and save plot
		pt_L->Draw();
		
		//Print fit parameters
		std::cout << "Fit L TPC : par0 = " << TPCFit->GetParameter(0) << " ; par1 = " << TPCFit->GetParameter(1) 
		<< " ; err0 = " << TPCFit->GetParError(0) << " ; err1 = " << TPCFit->GetParError(1) << std::endl;

		std::cout << "Lifetime Left TPC : " << TPCFit->GetParameter(1) << " \u00b1 " << TPCFit->GetParError(1) << std::endl;

		std::vector<double> fitParams_tDriftL = {TPCFit->GetParameter(1), TPCFit->GetParError(1)};


		c_tMPV->cd();
		pad2_tMPV->SetRightMargin(0.18);
		pad2_tMPV->SetLeftMargin(0.15);
		pad2_tMPV->Draw();
		pad2_tMPV->cd();
        MPVplot_R->SetMarkerColor(kAzure - 3);
        MPVplot_R->SetMarkerStyle(21);
        MPVplot_R->SetMarkerSize(0.5);
        MPVplot_R->GetXaxis()->SetTitle("t (ms)");
        MPVplot_R->SetTitle("dQ/dx MPV vs t Right TPC");
        MPVplot_R->GetYaxis()->SetTitle("dQ/dx MPV (ADC/cm)");
        MPVplot_R->Draw("AP");	

		TPaveText *pt_R = new TPaveText(.50,.75,.80,.88,"blNDC");
		gROOT->GetFunction("expo");
		TPCFit->SetParNames("Norm","Lifetime");
		TPCFit->SetLineColor(kGreen + 1);

		//Set initial parameters based on roughly the expected values (lifetime simulated as 10ms)
		TPCFit->SetParameter(0,1000);
		TPCFit->SetParameter(1,10);

		//Fit function
		std::cout << "----Fitting exponential in time---------------------------------------" << std::endl;
		MPVplot_R->Fit(TPCFit, "R");

		//Making stats box
		MPVplot_R->SetStats(0);

		std::cout << TPCFit->GetChisquare() << std::endl;
		std::cout << TPCFit->GetNDF() << std::endl;
		//Stats box position and text
		pt_R->SetBorderSize(1);
		pt_R->SetFillColor(0);
		pt_R->AddText(Form("AC cosmics: %d",track_count));
		((TText*)pt_R->GetListOfLines()->Last())->SetTextColor(1);
		pt_R->AddText(Form("#chi^{2} / DoF: %g / %d",TPCFit->GetChisquare(),TPCFit->GetNDF()));
		((TText*)pt_R->GetListOfLines()->Last())->SetTextColor(TPCFit->GetLineColor());
		pt_R->AddText(Form("Lifetime: %g#pm%g",TPCFit->GetParameter(1),TPCFit->GetParError(1)));
		((TText*)pt_R->GetListOfLines()->Last())->SetTextColor(TPCFit->GetLineColor());
		
		//Draw stats box and save plot
		pt_R->Draw();
		
		//Print fit parameters
		std::cout << "Fit L TPC : par0 = " << TPCFit->GetParameter(0) << " ; par1 = " << TPCFit->GetParameter(1) 
		<< " ; err0 = " << TPCFit->GetParError(0) << " ; err1 = " << TPCFit->GetParError(1) << std::endl;

		std::cout << "Lifetime Left TPC : " << TPCFit->GetParameter(1) << " \u00b1 " << TPCFit->GetParError(1) << std::endl;

		std::vector<double> fitParams_tDriftR = {TPCFit->GetParameter(1), TPCFit->GetParError(1)};

		c_tMPV->SaveAs((mydata_cosmics + std::to_string(N) + "_MPVt_multiwireTEST.png").c_str());

		N_wires.push_back(N);
		tau_xDrift_LR[0].push_back(fitParams_xDrift[0]);
		errTau_xDrift_LR[0].push_back(fitParams_xDrift[1]);
		tau_xDrift_LR[1].push_back(fitParams_xDrift[2]);
		errTau_xDrift_LR[1].push_back(fitParams_xDrift[3]);

		tau_tDrift_LR[0].push_back(fitParams_tDriftL[0]);
		errTau_tDrift_LR[0].push_back(fitParams_tDriftL[1]);
		tau_tDrift_LR[1].push_back(fitParams_tDriftR[0]);
		errTau_tDrift_LR[1].push_back(fitParams_tDriftR[1]);

		g_lifeVwiresXL->SetPoint(g_lifeVwiresXL->GetN(), N, fitParams_xDrift[0]);
		g_lifeVwiresXL->SetPointError(g_lifeVwiresXL->GetN()-1, 0., fitParams_xDrift[1]);
		g_lifeVwiresXR->SetPoint(g_lifeVwiresXR->GetN(), N, fitParams_xDrift[2]);
		g_lifeVwiresXR->SetPointError(g_lifeVwiresXR->GetN()-1, 0., fitParams_xDrift[3]);

		g_lifeVwiresTL->SetPoint(g_lifeVwiresTL->GetN(), N, fitParams_tDriftL[0]);
		g_lifeVwiresTL->SetPointError(g_lifeVwiresTL->GetN()-1, 0., fitParams_tDriftL[1]);
		g_lifeVwiresTR->SetPoint(g_lifeVwiresTR->GetN(), N, fitParams_tDriftR[0]);
		g_lifeVwiresTR->SetPointError(g_lifeVwiresTR->GetN()-1, 0., fitParams_tDriftR[1]);
		

		auto plotMPV_time = std::chrono::system_clock::now();
		elapsed_seconds = plotMPV_time - plotdQdx_time; 
    	std::cout << "++++++++++ Plot MPV time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	}

	auto postForLoop_time = std::chrono::system_clock::now();

	elapsed_seconds = postForLoop_time-preForLoop_time; 
    std::cout << "++++++++++ For Loop time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	TCanvas *C_mg = new TCanvas();
	C_mg->cd();
	TMultiGraph *mg = new TMultiGraph();
	g_lifeVwiresXL->SetMarkerColor(kRed);
	g_lifeVwiresTL->SetMarkerColor(kRed);
	g_lifeVwiresXR->SetMarkerColor(kGreen+1);
	g_lifeVwiresTR->SetMarkerColor(kGreen+1);
	g_lifeVwiresXL->SetMarkerStyle(21);
	g_lifeVwiresXR->SetMarkerStyle(21);
	g_lifeVwiresTL->SetMarkerStyle(47);
	g_lifeVwiresTR->SetMarkerStyle(47);
	g_lifeVwiresXL->SetTitle("x, left TPC");
	g_lifeVwiresXR->SetTitle("x, right TPC");
	g_lifeVwiresTL->SetTitle("t, left TPC");
	g_lifeVwiresTR->SetTitle("t, right TPC");
	mg->SetTitle("Lifetime vs number of combined wires;N; #tau (ms)");
	mg->Add(g_lifeVwiresXL);
	mg->Add(g_lifeVwiresXR);
	mg->Add(g_lifeVwiresTL);
	mg->Add(g_lifeVwiresTR);
	mg->Draw("AP");
	C_mg->BuildLegend(0.4,0.7,0.6,0.88);
	C_mg->SaveAs((mydata_cosmics + "LifetimeVwireNum2022AfullSelectedTEST.png").c_str());

	auto finalResults_time = std::chrono::system_clock::now();

	elapsed_seconds = finalResults_time-postForLoop_time; 
    std::cout << "++++++++++ Final plots time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;

	std::cout << "----------Lifetime Results-------------" << std::endl;
	std::cout << "----------Drift Distance-------------" << std::endl;
	std::cout << " N wires | Lifetime L | Error L | Lifetime R | Error R " << std::endl;

	for (int i = 0; i < 20; i++){
		
		std::cout << N_wires[i] << " | " << tau_xDrift_LR[0][i] << " | " << errTau_xDrift_LR[0][i] << " | " << tau_xDrift_LR[1][i] << " | " << errTau_xDrift_LR[1][i] << std::endl;
	}

	std::cout << "----------Drift Time-----------------" << std::endl;
	std::cout << " Lifetime L | Error L | Lifetime R | Error R " << std::endl;

	for (int i = 0; i < 20; i++){
		
		std::cout << N_wires[i] << " | " << tau_tDrift_LR[0][i] << " | " << errTau_tDrift_LR[0][i] << " | " << tau_tDrift_LR[1][i] << " | " << errTau_tDrift_LR[1][i] << std::endl;
	}

return 0;
	
}

//Function for taking a boolean input
int bool_input(std::string bool_name, int arg1, char**arg2){
	int input = 0;
	for(int i=0; i<arg1; ++i){
		if(!strcmp(arg2[i], bool_name.c_str())){
			input = std::stoi(arg2[i+1]);
		}
	}
	return input;
}

//The Landau-Gaussian function
double_t langaufun(Double_t *x, Double_t *par) {

	//Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
	//Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
	//Markus Friedl (Markus.Friedl@cern.ch) "langaus.C"

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

//Function to set LG fit initial values
void SetLGParameters(TH1D *h, Double_t *fp, Double_t *efp, Double_t &lb, Double_t &ub){
	//Adapted from Lan Nguyen's code which was adapted from Dom Barker's code etc.

    // Set fit parameters and errors  
    double nentries = h->GetEntries();
    double max = h->GetBinCenter(h->GetMaximumBin()); //x position of max dQdx value 
    double area = h->GetEntries()*h->GetBinWidth(1)*0.8 ; //normalisation constant //starting value for the fit - fit doesn't include whole hist
    double rms = h->GetRMS();
	
	//Starting guesses for parameters
    fp[0] = max / 100 * 4; //scale parameter of Landau density, starting guess about 5% of MPV x value
    fp[1] = max; //MPV: max dQdx value 		(position)
    fp[2] = area; // total area or integral //NORM (bin width) getNentries*Getbinwidth
    fp[3] = rms / 10;  //GSigma typically smaller than distribution rms
    
	//Starting guesses for errors
    efp[0] = max / 100 * 0.5; //0.5% of max dQdx value  	
    efp[1] = max / 100 * 8;	// 8% of max dQdx value
    efp[2] = area * 0.1;	
    efp[3] = rms * 0.01;
  
	//Lower and upper bounds for fit: sufficient to find the peak of the distribution precisely,
	//whilst staying in high stats bins to reduce statistical fluctuations
    double dQdxpeak = h->GetMaximumBin(); 
    lb = h->GetBinCenter(dQdxpeak-8);	
	ub = h->GetBinCenter(dQdxpeak+12);

}

//Function to fit the LG function to the data
TF1 *LGfitter(TH1 *h,  Double_t lbound, Double_t rbound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, Int_t ief){
	//Adapted by code from H.Pernegger (Heinz.Pernegger@cern.ch) and
	//Markus Friedl (Markus.Friedl@cern.ch): "langaus.C"
  
  // Fit the histogram h to a Landau function and pass back the parameters
  // and their errors.  
  
  //  Fit histogram to Landau/Gaussian conv
  Char_t FunName[100]; 
  sprintf(FunName,"Fitfcn_%s",h->GetName()); //(some kind of cycling through function names?)

  // The fitting function
  TF1 *ffit = new TF1(FunName,langaufun,lbound,rbound,4); 
  ffit->SetParameters(fitparams[0],fitparams[1],fitparams[2],fitparams[3]);
  ffit->SetParError(0,fiterrors[0]);
  if (fiterrors[0]==0) ffit->FixParameter(0,fitparams[0]); //if scale parameter is 0 it is fixed at 0
  ffit->SetParError(1,fiterrors[1]);
  ffit->SetParError(2,fiterrors[2]);
  ffit->SetParError(3,fiterrors[3]);

  //set LG scale parameter and GSigma limits to avoid too many or too few oscillations
  ffit->SetParLimits(0,20,60);
  ffit->SetParLimits(3,10,200);

  ffit->SetParNames("Width","MPV","TotalArea","GSigma");  
  

  TFitResultPtr r; 
  r = h->Fit(FunName,"LREQ");  //L = log likelihood method, E = error estimations using the Minos techniques, R = specied range, Q = quiet mode
  //Other fitting options https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html (7.1.1)
  
  //Fit parameter values if fit unsuccessful
  fiterrors[0] = -1000.0;
  fiterrors[1] = -1000.0;
  fiterrors[2] = -1000.0;
  fiterrors[3] = -1000.0;
  fitparams[0] = -1000.0;
  fitparams[1] = -1000.0;
  fitparams[2] = -1000.0;
  fitparams[3] = -1000.0;
  covmat[0] = -1000.0;
  covmat[1] = -1000.0;
  covmat[2] = -1000.0;  
  covmat[3] = -1000.0;  
  
  // Check fit status - options 'SUCCESSFUL', 'FAILED' and many many more.
  TString test =  gMinuit->fCstatu.Data();
  //std::cout << "LG fit " << test << std::endl;

  //Fill fit parameters if fit was successful
  if (test.BeginsWith("SUCC")) { 

  // Get Fit Parameters, their errors and cov matrix
  fitparams[0] = h->GetFunction(FunName)->GetParameter(0);
  fitparams[1] = h->GetFunction(FunName)->GetParameter(1);
  fitparams[2] = h->GetFunction(FunName)->GetParameter(2);
  fitparams[3] = h->GetFunction(FunName)->GetParameter(3);
  fiterrors[0] = h->GetFunction(FunName)->GetParError(0);
  fiterrors[1] = h->GetFunction(FunName)->GetParError(1);
  fiterrors[2] = h->GetFunction(FunName)->GetParError(2);
  fiterrors[3] = h->GetFunction(FunName)->GetParError(3);
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  TMatrixD matrix(4,4, fitter->GetCovarianceMatrix());
  covmat[0] = fitter->GetCovarianceMatrixElement(0, 1);
  covmat[1] = fitter->GetCovarianceMatrixElement(0, 2);
  covmat[2] = fitter->GetCovarianceMatrixElement(1, 2);
  covmat[3] = fitter->GetCovarianceMatrixElement(0, 3);
  //missing covariance terms here !
  

  }

  return(ffit);

}*/