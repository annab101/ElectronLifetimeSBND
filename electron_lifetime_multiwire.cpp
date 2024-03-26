//Author: Anna Beever
//Date:   December 2023
//Takes fileName.list as input, collects data with TChain and then
//does lifetime calc or dQ/dx, energy and angular plots as selected.

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
int bool_input(std::string bool_name, int arg1, char**arg2);
double_t langaufun(Double_t *x, Double_t *par);
void SetLGParameters(TH1D *h, Double_t *fp, Double_t *efp, Double_t &lb, Double_t &ub);
TF1 *LGfitter(TH1 *h,  Double_t lbound, Double_t rbound, Double_t *fitparams, Double_t *fiterrors, Double_t *covmat, Int_t ief);
void histLabels1D(TH1D *h, std::string units, int precision, bool brackets);
void histLabels2D(TH2D *h, std::string Xunits, std::string Yunits, int Xprecision, int Yprecision);
TH2D *hist_2D(std::string name, std::vector<double> x_data, std::vector<double> y_data, std::vector<double> x_binning, std::vector<double> y_binning, std::vector<std::string> labels);
void drawLifetimeStatBox(int track_count, double Lchisq, double Lndf, double Lpar1, double Lerr1, double Rchisq, double Rndf, double Rpar1, double Rerr1, TF1 *LeftTPCFit, TF1 *RightTPCFit);
void drawCosmicStats(TPaveText *stats_h, std::string text, int track_count);
TGraphErrors *findMPV(TGraphErrors *MPVplot, TH2D *h_dQdx, bool plot_projY, std::string mydata_cosmics, std::string dQdxProjYPlotName, bool XorT);
std::vector<double> fitExpoInX(TGraphErrors *MPVplot, int track_count);
std::vector<double> fitExpoInT(TF1 *TPCFit, TGraphErrors *MPVPlot_L, std::string t_formula, TPaveText *pt_L, int track_count, Color_t expoLineColor);
void MPVfancyDraw(TGraphErrors *MPVplot, bool x_var, int whichTPC = 0);


int main(int argc, char**argv) {

	auto start = std::chrono::system_clock::now();
	
	//Set options at command line. Note MPV not possible without dQdx
	bool plot_dQdx = bool_input("--dQdx", argc, argv);
	bool plot_projY = bool_input("--projY", argc, argv);
	
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

	std::string dQdxPlotName = "dQdx_hist_" + fileSaveID + ".png";
	std::string dQdxProjYPlotName = "dQdxProjY_hist_" + fileSaveID; //no .png since files labelled also by bin number
	std::string MPVPlotName = "MPV_" + fileSaveID + ".png";	
	std::string XvTPlotName = "XvT_" + fileSaveID + ".png";
	std::string dQdxVtPlotName = "dQdxVt_hist_" + fileSaveID + ".png";
	std::string dQdxVtProjYPlotName = "dQdxVtProjY_hist_" + fileSaveID;
	std::string MPVtPlotName = "MPVt_" + fileSaveID + ".png";

	//define location to save plots to (including plot name(!))
	std::string dQdxSaveLoc = mydata_cosmics + dQdxPlotName;
	std::string MPVSaveLoc = mydata_cosmics + MPVPlotName;
	std::string XvTSaveLoc = mydata_cosmics + XvTPlotName;
	std::string dQdxVtSaveLoc = mydata_cosmics + dQdxVtPlotName;
	std::string MPVtSaveLoc = mydata_cosmics + MPVtPlotName;

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
	TTreeReaderArray<Float_t> read_dqdx(treereader, "trk.hits2.dqdx");
	TTreeReaderArray<Float_t> read_T(treereader, "trk.hits2.h.time");
	TTreeReaderValue<Float_t> read_t0(treereader, "trk.t0");
	TTreeReaderValue<Float_t> read_xi(treereader, "trk.start.x");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<int> read_selected(treereader, "trk.selected");
	TTreeReaderArray<uint16_t> read_wire(treereader, "trk.hits2.h.wire");

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

	auto preForLoop_time = std::chrono::system_clock::now();
	
	for(int j=1; j<=20; j++){
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
			histLabels2D(h_dQdx_xDrift, "cm", "ADC/cm", 3, 3);
			c_dQdx_xDrift->SetRightMargin(0.18);
			c_dQdx_xDrift->cd();
			h_dQdx_xDrift->Draw("COLZ");
			TPaveText *stats_dQdx = new TPaveText(.55,.80,.80,.88,"blNDC");
			drawCosmicStats(stats_dQdx,"AC cosmics: %d",track_count);
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
			drawCosmicStats(stats_XvT,"AC cosmics: %d",track_count);		
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
			histLabels2D(h_dQdx_tDriftL, "ms", "ADC/cm", 3, 3);
			TPaveText *stats_dQdx_L = new TPaveText(.55,.80,.80,.88,"blNDC");
			drawCosmicStats(stats_dQdx_L,"AC cosmics: %d",track_count);
			c_dQdx_tDrift->cd();
			pad2->SetRightMargin(0.18);
			pad2->SetLeftMargin(0.15);
			pad2->Draw();
			pad2->cd();
			h_dQdx_tDriftR->SetStats(0);
			h_dQdx_tDriftR->GetXaxis()->SetTitle("t [ms]");
			h_dQdx_tDriftR->GetYaxis()->SetTitle("dQ/dx [ADc/cm]");
			h_dQdx_tDriftR->Draw("COLZ");
			histLabels2D(h_dQdx_tDriftR, "ms", "ADC/cm", 3, 3);
			TPaveText *stats_dQdx_R = new TPaveText(.55,.80,.80,.88,"blNDC");
			drawCosmicStats(stats_dQdx_R,"AC cosmics: %d",track_count);
			c_dQdx_tDrift->SaveAs((mydata_cosmics + std::to_string(N) + "_dQdxVt_multiwireTEST.png").c_str());

		}

		auto plotdQdx_time = std::chrono::system_clock::now();
		elapsed_seconds = plotdQdx_time-TTreeReader_time; 
    	std::cout << "++++++++++ Plot dQdx time: " << elapsed_seconds.count() << "s ++++++++++" << std::endl;


		std::cout << "----Finding MPV--------------------------------------------------------" << std::endl;
		
		TGraphErrors *MPVplot = new TGraphErrors(); 
		MPVplot = findMPV(MPVplot, h_dQdx_xDrift, plot_projY, mydata_cosmics, "xDrift_" + dQdxProjYPlotName, 0);

		//Plot and fit MPV
		TCanvas *c_MPV = new TCanvas();
		c_MPV->cd();
		MPVfancyDraw(MPVplot, 0); //0 for distance, 1 for time			
		std::vector<double> fitParams_xDrift = fitExpoInX(MPVplot, track_count);
		c_MPV->SaveAs((mydata_cosmics + std::to_string(N) + "_MPV_multiwireTEST.png").c_str());

		//Histogram
		TGraphErrors *MPVplot_L = new TGraphErrors(); 
		MPVplot_L = findMPV(MPVplot_L, h_dQdx_tDriftL, plot_projY, mydata_cosmics, "tDriftL" + dQdxProjYPlotName, 1);
		TGraphErrors *MPVplot_R = new TGraphErrors(); 
		MPVplot_R = findMPV(MPVplot_R, h_dQdx_tDriftR, plot_projY, mydata_cosmics, "tDriftR" + dQdxProjYPlotName, 1);


		TCanvas *c_tMPV = new TCanvas();
		c_tMPV->SetWindowSize(1500,500);
		TPad *pad1_tMPV = new TPad("pad1_tMPV", "", 0, 0, 0.5, 1.0);
		TPad *pad2_tMPV = new TPad("pad2_tMPV", "", 0.5, 0, 1.0, 1.0);
		c_tMPV->cd();
		pad1_tMPV->SetRightMargin(0.18);
		pad1_tMPV->SetLeftMargin(0.15);
		pad1_tMPV->Draw();
		pad1_tMPV->cd();
		MPVfancyDraw(MPVplot_L, 1, 0);

		//Function:
		std::string t_formula = "[0]*exp(-(x/[1]))";
		TF1 *TPCFit = new TF1("TPCFit", t_formula.c_str(), 0.2, 1.0);
		TPaveText *pt_L = new TPaveText(.50,.75,.80,.88,"blNDC");
		std::vector<double> fitParams_tDriftL = fitExpoInT(TPCFit, MPVplot_L, t_formula, pt_L, track_count, kRed);


		c_tMPV->cd();
		pad2_tMPV->SetRightMargin(0.18);
		pad2_tMPV->SetLeftMargin(0.15);
		pad2_tMPV->Draw();
		pad2_tMPV->cd();
		MPVfancyDraw(MPVplot_R, 1, 1);
		TPaveText *pt_R = new TPaveText(.50,.75,.80,.88,"blNDC");
		std::vector<double> fitParams_tDriftR = fitExpoInT(TPCFit, MPVplot_R, t_formula, pt_R, track_count, kGreen + 1);

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

void histLabels1D(TH1D *h, std::string units, int precision, bool brackets){
	std::stringstream ss;
	ss << std::setprecision(precision) << h->GetXaxis()->GetBinWidth(1);
	std::string yLabel;
	if(brackets){yLabel = "Entries / (" + ss.str() + " " + units + ")";}
	else{yLabel = "Entries / " + ss.str() + " " + units;}
	h->GetYaxis()->SetTitle(yLabel.c_str());
}

void histLabels2D(TH2D *h, std::string Xunits, std::string Yunits, int Xprecision, int Yprecision){
	std::stringstream ssX;
	std::stringstream ssY;
	ssX << std::setprecision(Xprecision) << h->GetXaxis()->GetBinWidth(1);
	ssY << std::setprecision(Yprecision) << h->GetYaxis()->GetBinWidth(1);

	std::string colorLabel = "Entries / (" + ssX.str() + " " + Xunits + " #times " + ssY.str() + " " + Yunits + ")";
	TLatex *colorTitle = new TLatex(0.95,0.5,colorLabel.c_str());
	colorTitle->SetNDC(1);
	colorTitle->SetTextFont(42);
	colorTitle->SetTextAngle(90);
	colorTitle->SetTextSize(0.035);
	colorTitle->SetTextAlign(22);
	colorTitle->Draw();
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

}

TH2D *hist_2D(std::string name, std::vector<double> x_data, std::vector<double> y_data, std::vector<double> x_binning, std::vector<double> y_binning, std::vector<std::string> labels){
	TH2D *h = new TH2D(name.c_str(),(labels[2] + " vs " + labels[0]).c_str(), x_binning[0], x_binning[1], x_binning[2], y_binning[0], y_binning[1], y_binning[2]);
	std::vector<double> weights(x_data.size(), 1.);
	h->FillN(weights.size(), x_data.data(), y_data.data(), weights.data());
	h->SetStats(0);
	h->GetXaxis()->SetTitle((labels[0] + " [" + labels[1] + "]").c_str());
	h->GetYaxis()->SetTitle((labels[2] + " [" + labels[3] + "]").c_str());
	histLabels2D(h, labels[1], labels[3], 3, 3);

	return h;
}

void drawLifetimeStatBox(int track_count, double Lchisq, double Lndf, double Lpar1, double Lerr1, double Rchisq, double Rndf, double Rpar1, double Rerr1, TF1 *LeftTPCFit, TF1 *RightTPCFit){
	TPaveText *pt = new TPaveText(.31,.65,.69,.88,"blNDC");
	pt->SetBorderSize(1);
	pt->SetFillColor(0);
	pt->AddText(Form("AC cosmics: %d",track_count));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(1);
	pt->AddText(Form("#chi^{2} / DoF: %g / %g",Lchisq,Lndf));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(LeftTPCFit->GetLineColor());
	pt->AddText(Form("Lifetime: %g#pm%g",Lpar1,Lerr1));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(LeftTPCFit->GetLineColor());
	pt->AddText(Form("#chi^{2} / DoF: %g / %g",Rchisq,Rndf));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(RightTPCFit->GetLineColor());
	pt->AddText(Form("Lifetime: %g#pm%g",Rpar1,Rerr1));
	((TText*)pt->GetListOfLines()->Last())->SetTextColor(RightTPCFit->GetLineColor());

	//Draw stats box and save plot
	pt->Draw();
}

void drawCosmicStats(TPaveText *stats_h, std::string text, int track_count){

	stats_h->SetBorderSize(1);
	stats_h->SetFillColor(0);
	stats_h->AddText(Form(text.c_str(),track_count));
	((TText*)stats_h->GetListOfLines()->Last())->SetTextColor(1);
	stats_h->Draw();

}

TGraphErrors *findMPV(TGraphErrors *MPVplot, TH2D *h_dQdx, bool plot_projY, std::string mydata_cosmics, std::string dQdxProjYPlotName, bool XorT){

	//Fit parameters
	double fp[4]; //fit parameters
	double efp[4]; //fit parameter errors
	double cov[4]; //covariance matrix
	double lb; //fit lower bound
	double ub; //fit upper bound

	//loop over each bin. For overflow bins change to i=0 to num_xbins + 2.
	for(int i=1; i<h_dQdx->GetNbinsX()+1; i++){
		
		//define filenames for projection plots
		std::string bin_number = "Bin" + std::to_string(i);		
		std::string saveLocTemp = mydata_cosmics + "ProjYplots/" + dQdxProjYPlotName + bin_number +".png";
		
		//perform the Landau-Gaussian fit
		TH1D *proj_y = h_dQdx->ProjectionY("proj_y", i, i);
		SetLGParameters(proj_y, fp, efp, lb, ub);
		TF1 *fitres = LGfitter(proj_y, lb, ub, fp, efp, cov, 0);

		//fill vectors with fit values
		int Nslice = proj_y->GetEntries(); //for stat error
	
		if(fp[1] >= 0){ //Get rid of any points where the fit failed, in which case fp[1] is set to -1000
			MPVplot->SetPoint(MPVplot->GetN(), h_dQdx->GetXaxis()->GetBinCenter(i), fp[1]);
			MPVplot->SetPointError(MPVplot->GetN()-1, 0., sqrt(pow(efp[1],2)+pow((fp[1]/sqrt(Nslice)),2)));
		}

		if(plot_projY){

			std::cout << "----Plotting projY-----------------------------------------------------" << std::endl;

			std::string plotTitleTemp;
			if(XorT == 0){
				plotTitleTemp = "dQ/dx for " + std::to_string(h_dQdx->GetXaxis()->GetBinLowEdge(i))
				+ " cm < x < " + std::to_string(h_dQdx->GetXaxis()->GetBinLowEdge(i) + h_dQdx->GetXaxis()->GetBinWidth(i)) + " cm";
			}
			else if(XorT == 1){
				plotTitleTemp = "dQ/dx for " + std::to_string(h_dQdx->GetXaxis()->GetBinLowEdge(i))
				+ " ms < t < " + std::to_string(h_dQdx->GetXaxis()->GetBinLowEdge(i) + h_dQdx->GetXaxis()->GetBinWidth(i)) + " ms";
			}

			TCanvas *c9 = new TCanvas();		
			proj_y->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
			histLabels1D(proj_y, "ADC/cm", 3, 1);			
			proj_y->SetTitle(plotTitleTemp.c_str()); 
			proj_y->Draw();
			fitres->Draw("same");
			gStyle->SetOptFit(1111); //stats box
			c9->SaveAs(saveLocTemp.c_str());
		}
		
	}

	return MPVplot;

}

std::vector<double> fitExpoInX(TGraphErrors *MPVplot, int track_count){
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
	drawLifetimeStatBox(track_count, LeftTPCFit->GetChisquare(), LeftTPCFit->GetNDF(), 
		Lpar1, Lerr1, RightTPCFit->GetChisquare(), RightTPCFit->GetNDF(), Rpar1, Rerr1, LeftTPCFit, RightTPCFit);
	
	//Print fit parameters
	std::cout << "Fit L TPC : par0 = " << Lpar0 << " ; par1 = " << Lpar1 << " ; err0 = " << Lerr0 << " ; err1 = " << Lerr1 << std::endl;
	std::cout << "Fit R TPC : par0 = " << Rpar0 << " ; par1 = " << Rpar1 << " ; err0 = " << Rerr0 << " ; err1 = " << Rerr1 << std::endl;

	std::cout << "Lifetime Left TPC : " << Lpar1<< " \u00b1 " << Lerr1 << "Norm Left TPC : " << Lpar0 << " \u00b1 " << Lerr0 << std::endl;
	std::cout << "Lifetime Right TPC : " << Rpar1 <<  " \u00b1 " << Rerr1 << "Norm Right TPC : " << Rpar0 << " \u00b1 " << Rerr0 << std::endl;

	std::vector<double> results_x = {Lpar1, Lerr1, Rpar1, Rerr1};

	return results_x;
}

std::vector<double> fitExpoInT(TF1 *TPCFit, TGraphErrors *MPVplot_L, std::string t_formula, TPaveText *pt_L, int track_count, Color_t expoLineColor){

	gROOT->GetFunction("expo");
	TPCFit->SetParNames("Norm","Lifetime");

	TPCFit->SetLineColor(expoLineColor);

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

	std::vector<double> results_T = {TPCFit->GetParameter(1), TPCFit->GetParError(1)};

	return results_T;

}

void MPVfancyDraw(TGraphErrors *MPVplot, bool x_var, int whichTPC){ //0 for distance, 1 for time
	MPVplot->SetMarkerColor(kAzure - 3);
	MPVplot->SetMarkerStyle(21);
	MPVplot->SetMarkerSize(0.5);
	if(x_var == 0){
		MPVplot->GetXaxis()->SetTitle("x (cm)");
		MPVplot->SetTitle("dQ/dx MPV vs x");
	}
	else if(x_var == 1){
		MPVplot->GetXaxis()->SetTitle("t (ms)");
		if(whichTPC==0){
			MPVplot->SetTitle("dQ/dx MPV vs t Left TPC");
		}
		else if (whichTPC==1){
			MPVplot->SetTitle("dQ/dx MPV vs t Right TPC");
		}
	}
	MPVplot->GetYaxis()->SetTitle("dQ/dx MPV (ADC/cm)");
	MPVplot->Draw("AP");		
}