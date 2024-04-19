//Author: Anna Beever
//Date:   December 2023
//Takes fileName.list as input, collects data with TChain and then
//does lifetime calc or dQ/dx, energy and angular plots as selected.

//C++ includes
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
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

//Local includes
#include "Helpers/PlottingHelpers.h"
#include "Helpers/PlottingHelpers.cpp"
#include "Helpers/FittingHelpers.h"
#include "Helpers/FittingHelpers.cpp"
#include "Helpers/StatsHelpers.h"
#include "Helpers/StatsHelpers.cpp"

using namespace calib;


int main(int argc, char**argv) {
	
	
	//read in filename (data), fileSaveID (name of files saved), filesaveloc (where files saved), 
	//codeConfig (0 = lifetime, 1 = lifetime with wires grouped, 2 = lifetime with hits grouped)
	std::string filename = "no file";
	std::string fileSaveID = "noID";
	std::string fileSaveLoc = "noLoc";
	int codeConfig = 0;
	bool plotProjY = false;

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
		if(!strcmp(argv[i], "--config")){
			codeConfig = std::stoi(argv[i+1]);
		}
	}
	
	//check input files have been provided
	if(!strcmp("no file", filename.c_str()) || !strcmp("noID", fileSaveID.c_str()) || !strcmp("noLoc", fileSaveLoc.c_str())){
		std::cout << "!!!!ERROR: MISSING INPUT INFO (file name and ID for saving)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		exit(0);
	}

	std::cout << "----Running configuration " << codeConfig << " ---------------------------------" << std::endl;
	
	std::cout << "----Retrieving file----------------------------------------------------" << std::endl;

	//define names of plots
	std::string mydata_cosmics = "/exp/sbnd/data/users/abeever/cosmics_analysis/";

	//Read in data files line by line and chain together
	TChain *chain = new TChain;

	std::fstream file_list;
	file_list.open(filename.c_str(),std::ios::in);
	
	if (file_list.is_open()){
		std::string fileUrl;
		while(getline(file_list, fileUrl)){

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


	//HISTOGRAMS AND CONTAINERS FOR RESULTS

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

	TH2D* h_dQdx_xDrift_basic = new TH2D("h_dQdx_xDrift_basic", "dQ/dx vs x", N_bins, -200, 200, 75, 200, 1800);
	TH2D* h_dQdx_tDriftL_basic = new TH2D("h_dQdx_tDriftL_basic", "dQ/dx vs t Left TPC", N_bins, 0, 1.3, 75, 200, 1800);
	TH2D* h_dQdx_tDriftR_basic = new TH2D("h_dQdx_tDriftR_basic", "dQ/dx vs t Right TPC", N_bins, 0, 1.3, 75, 200, 1800);

	TH2D* h_dQdx_xDrift_basic_corrected = new TH2D("h_dQdx_xDrift_basic_corrected", "dQ/dx vs x", N_bins, -200, 200, 75, 200, 1800);

	TH1D** projY_xDrift = new TH1D*[maxWireGroup];
	TH1D** projY_tDriftL = new TH1D*[maxWireGroup];
	TH1D** projY_tDriftR = new TH1D*[maxWireGroup];

	TF1** LGfit_xDrift = new TF1*[maxWireGroup];
	TF1** LGfit_tDriftL = new TF1*[maxWireGroup];
	TF1** LGfit_tDriftR = new TF1*[maxWireGroup];

	TH1D* projY_xDrift_basic = new TH1D();
	TH1D* projY_tDriftL_basic = new TH1D();
	TH1D* projY_tDriftR_basic = new TH1D();

	TF1* LGfit_xDrift_basic = new TF1();
	TF1* LGfit_tDriftL_basic = new TF1();
	TF1* LGfit_tDriftR_basic = new TF1();

	fitLGParameters fp_xDrift;
	fitLGParameters fp_tDriftL;
	fitLGParameters fp_tDriftR;

	fitExpoParameters expoParams_xDriftL;
	fitExpoParameters expoParams_xDriftR;
	fitExpoParameters expoParams_tDriftL;
	fitExpoParameters expoParams_tDriftR;

	fitExpoConstParameters expoConstParams;

	TGraphErrors** MPVplot_xDrift = new TGraphErrors*[maxWireGroup];
	TGraphErrors** MPVplot_tDriftL = new TGraphErrors*[maxWireGroup];
	TGraphErrors** MPVplot_tDriftR = new TGraphErrors*[maxWireGroup];

	for(int i = 1; i <= maxWireGroup; i++){
		MPVplot_xDrift[i-1] = new TGraphErrors(); 
		MPVplot_tDriftL[i-1] = new TGraphErrors();
		MPVplot_tDriftR[i-1] = new TGraphErrors();
	}

	TGraphErrors *MPVplot_xDrift_basic = new TGraphErrors();
	TGraphErrors *MPVplot_tDriftL_basic = new TGraphErrors();
	TGraphErrors *MPVplot_tDriftR_basic = new TGraphErrors();

	TF1 *expoFit_xDriftL = new TF1();
	TF1 *expoFit_xDriftR = new TF1();
	TF1 *expoFit_tDriftL = new TF1();
	TF1 *expoFit_tDriftR = new TF1();

	TF1 *expoConstFit = new TF1();

	TCanvas *c_plain = new TCanvas();
	c_plain->SetLeftMargin(0.15);
	c_plain->SetRightMargin(0.18);
	c_plain->SetBottomMargin(0.12);

	TCanvas *c_split = new TCanvas();
	c_split->SetWindowSize(2000,500);
	TPad *pad1 = new TPad("pad1", "", 0, 0, 0.5, 1.0);
	TPad *pad2 = new TPad("pad2", "", 0.5, 0, 1.0, 1.0);
	c_split->cd();
	pad1->SetRightMargin(0.20);
	pad1->SetLeftMargin(0.18);
	pad2->SetRightMargin(0.20);
	pad2->SetLeftMargin(0.18);


	TCanvas *c_split2 = new TCanvas();
	c_split2->SetWindowSize(2000,500);
	TPad *pad12 = new TPad("pad12", "", 0, 0, 0.5, 1.0);
	TPad *pad22 = new TPad("pad22", "", 0.5, 0, 1.0, 1.0);
	c_split2->cd();
	pad12->SetRightMargin(0.1);
	pad12->SetLeftMargin(0.18);
	pad22->SetRightMargin(0.1);
	pad22->SetLeftMargin(0.18);

	//Sheffield colour scheme
	Int_t deepViolet = TColor::GetColor("#440099");
	Int_t coral = TColor::GetColor("#E7004C");
	Int_t powderBlue = TColor::GetColor("#9ADBE8");
	Int_t flamingo = TColor::GetColor("#FF6371");
	Int_t mintGreen = TColor::GetColor("#00CE7C");


	treereader.Restart();
	while(treereader.Next()){

		if(*read_selected != 1){
			continue;
		}

		if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
		track_count ++;	

		
		auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
		int start_index = std::distance(read_x.begin(), it_start_x);

		if(start_index == read_x.GetSize()){continue;}
		auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
		int end_index = std::distance(read_x.begin(), it_end_x);
		if(end_index < 0){continue;}

		int minWire = read_wire[start_index];
		int maxWire = read_wire[end_index];

		if(minWire > maxWire){
			std::swap(minWire, maxWire);
		}

		std::valarray<double> dQdx_sum(0.,maxWireGroup);
		std::valarray<double> x_sum(0.,maxWireGroup);
		std::valarray<double> t_sum(0.,maxWireGroup);
		std::valarray<int> count(0,maxWireGroup);

		
		if(codeConfig == 0){
			for(int i = start_index; i <= end_index; i++){

				h_dQdx_xDrift_basic->Fill(read_x[i], read_dqdx[i]);
				h_dQdx_xDrift_basic_corrected->Fill(read_x[i], read_dqdx[i]*exp((200-TMath::Abs(read_x[i]))/(12.8*156.267)));
				
										
				if(-200. < read_x[i] && read_x[i] < 0.){
					h_dQdx_tDriftL_basic->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

				if(0. < read_x[i] && read_x[i] < 200.){
					h_dQdx_tDriftR_basic->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

			}
		}
		if(codeConfig == 1){
			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

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
		}
		if(codeConfig == 2){
			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

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
				else{

					for (int N = 1; N <= maxWireGroup; N++){
						
						int j = N-1;

						if(count[j] == N){

							h_dQdx_xDrift[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
							if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
								h_dQdx_tDriftL[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
							}

							if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
								h_dQdx_tDriftR[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
							}

							dQdx_sum[j] = 0.;
							x_sum[j] = 0.;
							t_sum[j] = 0.;
							count[j] = 0;

						}
					}

				}

			
			}
		}

	}

	//Basic analysis
	if(codeConfig == 0){

		openAndClear(c_plain);
		h_dQdx_xDrift_basic_corrected->SetStats(0);
		h_dQdx_xDrift_basic_corrected->SetTitle(";x (cm);dQ/dx (ADC/cm);Entries");
		h_dQdx_xDrift_basic_corrected->GetXaxis()->SetLabelOffset(0.01);
		h_dQdx_xDrift_basic_corrected->GetXaxis()->SetNdivisions(505);
		h_dQdx_xDrift_basic_corrected->GetYaxis()->SetLabelOffset(0.01);
		h_dQdx_xDrift_basic_corrected->GetYaxis()->SetNdivisions(505);
		h_dQdx_xDrift_basic_corrected->GetYaxis()->SetTitleOffset(1.8);
		setFontSize<TH2>(h_dQdx_xDrift_basic_corrected, 133, 25);
		setFontSizeZ(h_dQdx_xDrift_basic_corrected, 133, 25);
		h_dQdx_xDrift_basic_corrected->GetZaxis()->CenterTitle(true);
		h_dQdx_xDrift_basic_corrected->GetZaxis()->SetTitleOffset(2.0);
		h_dQdx_xDrift_basic_corrected->Draw("COLZ");
		saveFig(c_plain, mydata_cosmics + fileSaveLoc + "dQdx_xDrift_corrected_" + fileSaveID);

		openAndClear(c_plain);
		h_dQdx_xDrift_basic->SetStats(0);
		h_dQdx_xDrift_basic->SetTitle(";x (cm);dQ/dx (ADC/cm);Entries");
		h_dQdx_xDrift_basic->GetXaxis()->SetLabelOffset(0.01);
		h_dQdx_xDrift_basic->GetXaxis()->SetNdivisions(505);
		h_dQdx_xDrift_basic->GetYaxis()->SetLabelOffset(0.01);
		h_dQdx_xDrift_basic->GetYaxis()->SetNdivisions(505);
		h_dQdx_xDrift_basic->GetYaxis()->SetTitleOffset(1.8);
		setFontSize<TH2>(h_dQdx_xDrift_basic, 133, 25);
		setFontSizeZ(h_dQdx_xDrift_basic, 133, 25);
		h_dQdx_xDrift_basic->GetZaxis()->CenterTitle(true);
		h_dQdx_xDrift_basic->GetZaxis()->SetTitleOffset(2.0);
		h_dQdx_xDrift_basic->Draw("COLZ");
		saveFig(c_plain, mydata_cosmics + fileSaveLoc + "dQdx_xDrift_" + fileSaveID);

		h_dQdx_tDriftL_basic->SetStats(0);
		h_dQdx_tDriftL_basic->SetTitle(";t (ms);dQ/dx (ADC/cm);Entries");
		h_dQdx_tDriftL_basic->GetXaxis()->SetLabelOffset(0.01);
		h_dQdx_tDriftL_basic->GetXaxis()->SetNdivisions(505);
		h_dQdx_tDriftL_basic->GetYaxis()->SetLabelOffset(0.01);
		h_dQdx_tDriftL_basic->GetYaxis()->SetNdivisions(505);
		h_dQdx_tDriftL_basic->GetYaxis()->SetTitleOffset(1.8);
		setFontSize<TH2>(h_dQdx_tDriftL_basic, 133, 25);
		setFontSizeZ(h_dQdx_tDriftL_basic, 133, 25);
		h_dQdx_tDriftL_basic->GetZaxis()->CenterTitle(true);
		h_dQdx_tDriftL_basic->GetZaxis()->SetTitleOffset(2.0);

		h_dQdx_tDriftR_basic->SetStats(0);
		h_dQdx_tDriftR_basic->SetTitle(";t (ms);dQ/dx (ADC/cm);Entries");
		h_dQdx_tDriftR_basic->GetXaxis()->SetLabelOffset(0.01);
		h_dQdx_tDriftR_basic->GetXaxis()->SetNdivisions(505);
		h_dQdx_tDriftR_basic->GetYaxis()->SetLabelOffset(0.01);
		h_dQdx_tDriftR_basic->GetYaxis()->SetNdivisions(505);
		h_dQdx_tDriftR_basic->GetYaxis()->SetTitleOffset(1.8);
		setFontSize<TH2>(h_dQdx_tDriftR_basic, 133, 25);
		setFontSizeZ(h_dQdx_tDriftR_basic, 133, 25);
		h_dQdx_tDriftR_basic->GetZaxis()->CenterTitle(true);
		h_dQdx_tDriftR_basic->GetZaxis()->SetTitleOffset(2.0);

		c_split->cd();		
		pad1->Draw();		
		pad1->cd();	
		h_dQdx_tDriftL_basic->Draw("COLZ");
		c_split->cd();
		pad2->Draw();
		pad2->cd();
		h_dQdx_tDriftR_basic->Draw("COLZ");
		saveFig(c_split, mydata_cosmics + fileSaveLoc + "dQdx_tDrift_" + fileSaveID);


		

		expoParams_xDriftL.reset();
		expoParams_xDriftR.reset();
		expoParams_tDriftL.reset();
		expoParams_tDriftR.reset();

		for(int j = 1; j <= N_bins; j++){

			fp_xDrift.reset();
			fp_tDriftL.reset();
			fp_tDriftR.reset();

			projY_xDrift_basic = h_dQdx_xDrift_basic->ProjectionY("projY_xDrift_basic", j, j);
			projY_tDriftL_basic = h_dQdx_tDriftL_basic->ProjectionY("projY_tDriftL_basic", j, j);
			projY_tDriftR_basic = h_dQdx_tDriftR_basic->ProjectionY("projY_tDriftR_basic", j, j);

			SetLGParameters(projY_xDrift_basic, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.lb, fp_xDrift.ub);
			SetLGParameters(projY_tDriftL_basic, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.lb, fp_tDriftL.ub);
			SetLGParameters(projY_tDriftR_basic, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.lb, fp_tDriftR.ub);

			LGfit_xDrift_basic = fitter(projY_xDrift_basic, fp_xDrift.lb, fp_xDrift.ub, fp_xDrift.fp, fp_xDrift.efp, fp_xDrift.cov, "LG");
			LGfit_tDriftL_basic = fitter(projY_tDriftL_basic, fp_tDriftL.lb, fp_tDriftL.ub, fp_tDriftL.fp, fp_tDriftL.efp, fp_tDriftL.cov, "LG");
			LGfit_tDriftR_basic = fitter(projY_tDriftR_basic, fp_tDriftR.lb, fp_tDriftR.ub, fp_tDriftR.fp, fp_tDriftR.efp, fp_tDriftR.cov, "LG");


			if(fp_xDrift.fp[1] >= 0.){

				MPVplot_xDrift_basic->SetPoint(MPVplot_xDrift_basic->GetN(), h_dQdx_xDrift_basic->GetXaxis()->GetBinCenter(j), fp_xDrift.fp[1]);
				MPVplot_xDrift_basic->SetPointError(MPVplot_xDrift_basic->GetN() - 1, 0., pointError(projY_xDrift_basic, fp_xDrift));
			}

			if(fp_tDriftL.fp[1] >= 0.){

				MPVplot_tDriftL_basic->SetPoint(MPVplot_tDriftL_basic->GetN(), h_dQdx_tDriftL_basic->GetXaxis()->GetBinCenter(j), fp_tDriftL.fp[1]);
				MPVplot_tDriftL_basic->SetPointError(MPVplot_tDriftL_basic->GetN() - 1, 0., pointError(projY_tDriftL_basic, fp_tDriftL));
			}

			if(fp_tDriftR.fp[1] >= 0.){

				MPVplot_tDriftR_basic->SetPoint(MPVplot_tDriftR_basic->GetN(), h_dQdx_tDriftR_basic->GetXaxis()->GetBinCenter(j), fp_tDriftR.fp[1]);
				MPVplot_tDriftR_basic->SetPointError(MPVplot_tDriftR_basic->GetN() - 1, 0., pointError(projY_tDriftR_basic, fp_tDriftR));
			}

			if(j==35 && plotProjY == true){
				std::string projYtitle;
				std::stringstream s1;
				s1 << std::setprecision(4) << h_dQdx_xDrift_basic->GetXaxis()->GetBinLowEdge(j);
				std::string str1 = s1.str();
				std::stringstream s2;
				s2 << std::setprecision(4) << (h_dQdx_xDrift_basic->GetXaxis()->GetBinLowEdge(j) + h_dQdx_xDrift_basic->GetXaxis()->GetBinWidth(j));
				std::string str2 = s2.str();
				projYtitle = "dQ/dx for " + str1 + " cm < x < " + str2 + " cm;dQ/dx (ADC/cm);Entries";
				projY_xDrift_basic->SetTitle(projYtitle.c_str());
				projY_xDrift_basic->SetTitleFont(133);
				projY_xDrift_basic->SetTitleSize(28);

				projY_xDrift_basic->GetXaxis()->SetLabelOffset(0.01);
				projY_xDrift_basic->GetXaxis()->SetNdivisions(505);
				setFontSize<TH1>(projY_xDrift_basic, 133, 25);
				projY_xDrift_basic->GetYaxis()->SetLabelOffset(0.01);

				projY_xDrift_basic->SetFillColor(flamingo);
				projY_xDrift_basic->SetLineColor(flamingo);
				projY_xDrift_basic->SetStats(0);
				LGfit_xDrift_basic->SetLineColor(deepViolet);
				LGfit_xDrift_basic->SetLineWidth(2.3);

				openAndClear(c_plain); 
				projY_xDrift_basic->Draw();
				LGfit_xDrift_basic->Draw("same");
				saveFig(c_plain, mydata_cosmics + fileSaveLoc + "projY_xDrift_" + fileSaveID);
			}

		}

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

		expoFit_xDriftL = fitter(MPVplot_xDrift_basic, expoParams_xDriftL.lb, expoParams_xDriftL.ub, expoParams_xDriftL.fp, expoParams_xDriftL.efp, expoParams_xDriftL.cov, "expoX");
		expoFit_xDriftR = fitter(MPVplot_xDrift_basic, expoParams_xDriftR.lb, expoParams_xDriftR.ub, expoParams_xDriftR.fp, expoParams_xDriftR.efp, expoParams_xDriftR.cov, "expoX");
		expoFit_tDriftL = fitter(MPVplot_tDriftL_basic, expoParams_tDriftL.lb, expoParams_tDriftL.ub, expoParams_tDriftL.fp, expoParams_tDriftL.efp, expoParams_tDriftL.cov, "expoT");
		expoFit_tDriftR = fitter(MPVplot_tDriftR_basic, expoParams_tDriftR.lb, expoParams_tDriftR.ub, expoParams_tDriftR.fp, expoParams_tDriftR.efp, expoParams_tDriftR.cov, "expoT");


		std::cout << "====================================================================" << std::endl;
		std::cout << "===============================RESULTS==============================" << std::endl;
		std::cout << "Lifetime, x, L: " << expoFit_xDriftL->GetParameter(1) << " \u00b1 " << expoFit_xDriftL->GetParError(1) << std::endl;
		std::cout << "Lifetime, x, R: " << expoFit_xDriftR->GetParameter(1) << " \u00b1 " << expoFit_xDriftR->GetParError(1) << std::endl;
		std::cout << "Lifetime, t, L: " << expoFit_tDriftL->GetParameter(1) << " \u00b1 " << expoFit_tDriftL->GetParError(1) << std::endl;
		std::cout << "Lifetime, t, R: " << expoFit_tDriftR->GetParameter(1) << " \u00b1 " << expoFit_tDriftR->GetParError(1) << std::endl;

		//PLOTS
		
		openAndClear(c_plain);
		setMarker<TGraphErrors>(MPVplot_xDrift_basic, kAzure - 3, 21, 0.6);
		MPVplot_xDrift_basic->SetTitle(";x (cm);dQ/dx MPV (ADC/cm)");
		MPVplot_xDrift_basic->GetXaxis()->SetNdivisions(505);
		setFontSize<TGraphErrors>(MPVplot_xDrift_basic, 133, 25);
		MPVplot_xDrift_basic->Draw("AP");
		expoFit_xDriftL->SetLineColor(coral);
		expoFit_xDriftL->SetLineWidth(2.3);
		expoFit_xDriftL->Draw("SAME");
		expoFit_xDriftR->SetLineColor(deepViolet);
		expoFit_xDriftR->SetLineWidth(2.3);
		expoFit_xDriftR->Draw("SAME");
		TPaveText *stats = statsBox({.34,.62,.63,.88}, track_count, expoFit_xDriftL, expoFit_xDriftR);
		stats->Draw();
		saveFig(c_plain, mydata_cosmics + fileSaveLoc + "MPV_xDrift_" + fileSaveID);

		setMarker<TGraphErrors>(MPVplot_tDriftL_basic, kAzure - 3, 21, 0.6);
		MPVplot_tDriftL_basic->SetTitle("dQ/dx MPV vs t Left TPC;t (ms); dQ/dx MPV (ADC/cm)");
		expoFit_tDriftL->SetLineColor(coral);
		setFontSize<TGraphErrors>(MPVplot_tDriftL_basic, 133, 25);
		TPaveText *statsL = statsBox({.45,.65,.80,.88}, track_count, expoFit_tDriftL);
		setMarker<TGraphErrors>(MPVplot_tDriftR_basic, kAzure - 3, 21, 0.5);
		MPVplot_tDriftR_basic->SetTitle("dQ/dx MPV vs t Right TPC;t (ms); dQ/dx MPV (ADC/cm)");
		expoFit_tDriftR->SetLineColor(deepViolet);
		setFontSize<TGraphErrors>(MPVplot_tDriftR_basic, 133, 25);
		TPaveText *statsR = statsBox({.45,.65,.80,.88}, track_count, expoFit_tDriftR);
		
		c_split2->cd();		
		pad12->Draw();		
		pad12->cd();	
		MPVplot_tDriftL_basic->Draw("AP");	
		expoFit_tDriftL->Draw("SAME");
		statsL->Draw();
		c_split2->cd();
		pad22->Draw();
		pad22->cd();
		MPVplot_tDriftR_basic->Draw("AP");
		expoFit_tDriftR->Draw("SAME");
		statsR->Draw();
		saveFig(c_split2, mydata_cosmics + fileSaveLoc + "MPV_tDrift_" + fileSaveID);
			
	}

	//Multiwire analysis
	if(codeConfig == 1 || codeConfig == 2){
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
			setMarker<TGraphErrors>(MPVplot_xDrift[i-1], kAzure - 3, 21, 0.5);
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

		}

		openAndClear(c_plain);
		TMultiGraph *mg = new TMultiGraph();
		setMarker<TGraphErrors>(g_lifeVwiresXL, coral, 21, 1.0);
		setMarker<TGraphErrors>(g_lifeVwiresXR, deepViolet, 21, 1.0);
		setMarker<TGraphErrors>(g_lifeVwiresTL, coral, 47, 1.0);
		setMarker<TGraphErrors>(g_lifeVwiresTR, deepViolet, 47, 1.0);

		//g_lifeVwiresXL->SetTitle("TPC East");
		//g_lifeVwiresXR->SetTitle("TPC West");
		g_lifeVwiresXL->SetTitle("x, left TPC");
		g_lifeVwiresXR->SetTitle("x, right TPC");
		g_lifeVwiresTL->SetTitle("t, left TPC");
		g_lifeVwiresTR->SetTitle("t, right TPC");

		//mg->SetTitle("Lifetime vs number of combined wires;N; #tau (ms)");
		mg->SetTitle(";N;#tau (ms)");
		mg->Add(g_lifeVwiresXL);
		mg->Add(g_lifeVwiresXR);
		//mg->Add(g_lifeVwiresTL);
		//mg->Add(g_lifeVwiresTR);
		setFontSize<TMultiGraph>(mg, 133, 25);
		mg->Draw("AP");
		c_plain->BuildLegend(0.4,0.7,0.6,0.88);
		if(codeConfig == 1){
			saveFig(c_plain, mydata_cosmics + fileSaveLoc + "MultiWirePlots/LifetimeVwireNum" + fileSaveID);
		}
		if(codeConfig == 2){
			saveFig(c_plain, mydata_cosmics + fileSaveLoc + "MultiWirePlots/LifetimeVnumHits" + fileSaveID);
		}

		expoConstParams.reset();
		SetExpoConstParameters(expoConstParams.fp);
		expoConstParams.lb = 0.5;
		expoConstParams.ub = 20.5;
		expoConstFit = fitter(g_lifeVwiresXL,expoConstParams.lb, expoConstParams.ub, expoConstParams.fp, expoConstParams.efp, expoConstParams.cov, "expoConst");
		expoConstFit->SetLineColor(mintGreen);
		expoConstFit->SetLineWidth(2.3);

		openAndClear(c_plain);
		g_lifeVwiresXL->SetTitle(";N;#tau (ms)");
		setFontSize<TGraphErrors>(g_lifeVwiresXL, 133, 25);
		g_lifeVwiresXL->Draw("AP");
		expoConstFit->Draw("SAME");
		saveFig(c_plain, mydata_cosmics + fileSaveLoc + "diffusionExpo_" + fileSaveID);

		std::cout << expoConstFit->GetParameter(1) << " pm " << expoConstFit->GetParError(1) << std::endl;

	}
	
	return 0;
}
