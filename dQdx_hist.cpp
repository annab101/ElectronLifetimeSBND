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
	
	//codeConfig (0 = lifetime, 1 = lifetime with wires grouped, 2 = lifetime with hits grouped, 3 = angular bins, 4 = wires and angular bins, 5 = octants)
	int choppyLim = -1;
	int codeConfig = 0;
	int maxWireGroup = 20;
	int NBinsX = 100;
	int NBinsT = 100;
	int NBinsHalfX = 50;
	int NBinsdQdx = 75;
	double minX = -200.;
	double halfX = 0.;
	double maxX = 200.;
	double mindQdx = 200.;
	double maxdQdx = 1800.;
	double minT = 0.;
	double maxT = 1.3;	
	int angBins = 6;
	std::string configLabel = "noConfigLabel";

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("choppyLim", choppyLim);
	p->getValue("codeConfig", codeConfig);
	p->getValue("maxWireGroup", maxWireGroup);
	p->getValue("NBinsX", NBinsX);
	p->getValue("NBinsT", NBinsT);
	p->getValue("NBinsHalfX", NBinsHalfX);
	p->getValue("NBinsdQdx", NBinsdQdx);
	p->getValue("minX", minX);
	p->getValue("halfX", halfX);
	p->getValue("maxX", maxX);
	p->getValue("mindQdx", mindQdx);
	p->getValue("maxdQdx", maxdQdx);
	p->getValue("minT", minT);
	p->getValue("maxT", maxT);
	p->getValue("angBins", angBins);
	configLabel = configMap[codeConfig];

	//Read in data files line by line and chain together
	TChain *chain = new TChain;

	std::fstream file_list;
	file_list.open(inputData.c_str(),std::ios::in);
	
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
	TTreeReaderValue<Float_t> read_yi(treereader, "trk.start.y");
	TTreeReaderValue<Float_t> read_zi(treereader, "trk.start.z");
	TTreeReaderValue<Float_t> read_xf(treereader, "trk.end.x");
	TTreeReaderValue<Float_t> read_yf(treereader, "trk.end.y");
	TTreeReaderValue<Float_t> read_zf(treereader, "trk.end.z");
	TTreeReaderValue<int> read_event(treereader, "trk.meta.evt");
	TTreeReaderValue<int> read_selected(treereader, "trk.selected");
	TTreeReaderArray<uint16_t> read_wire(treereader, "trk.hits2.h.wire");

	std::filesystem::create_directory((saveLoc + dataset + "_" + configLabel).c_str());

	TFile f((saveLoc + dataset  + "_" + configLabel + "/dQdx_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");

	int track_count = 0;
	TH1D* h_stats = new TH1D("h_stats", "Number of tracks", 9, 0, 9);

	if(codeConfig == 0){

		TH2D* h_dQdx_xDrift_basic = new TH2D("h_dQdx_xDrift_basic", "dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
		TH2D* h_dQdx_tDriftL_basic = new TH2D("h_dQdx_tDriftL_basic", "dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
		TH2D* h_dQdx_tDriftR_basic = new TH2D("h_dQdx_tDriftR_basic", "dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);

			auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
			int start_index = std::distance(read_x.begin(), it_start_x);

			if(start_index == read_x.GetSize()){continue;}
			auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
			int end_index = std::distance(read_x.begin(), it_end_x);
			if(end_index < 0){continue;}

			
			for(int i = start_index; i <= end_index; i++){

				h_dQdx_xDrift_basic->Fill(read_x[i], read_dqdx[i]);
				
										
				if(-200. < read_x[i] && read_x[i] < 0.){
					h_dQdx_tDriftL_basic->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

				if(0. < read_x[i] && read_x[i] < 200.){
					h_dQdx_tDriftR_basic->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

			}

		}

		h_dQdx_xDrift_basic->Write();
		h_dQdx_tDriftL_basic->Write();
		h_dQdx_tDriftR_basic->Write();
		h_stats->Write();

		f.Close();
				
	}

	if(codeConfig == 1){

		TH2D** h_dQdx_xDrift_wires = new TH2D*[maxWireGroup];
		TH2D** h_dQdx_tDriftL_wires = new TH2D*[maxWireGroup];
		TH2D** h_dQdx_tDriftR_wires = new TH2D*[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){

			h_dQdx_xDrift_wires[i-1] = new TH2D(TString::Format("h_dQdx_xDrift_%dwires", i),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftL_wires[i-1] = new TH2D(TString::Format("h_dQdx_tDriftL_%dwires", i),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftR_wires[i-1] = new TH2D(TString::Format("h_dQdx_tDriftR_%dwires", i),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
		
		}

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			//Track counter
			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);

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

			
			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

					for (int N = 1; N <= maxWireGroup; N++){
					
						int j = N-1;

						h_dQdx_xDrift_wires[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
										
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

						h_dQdx_xDrift_wires[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

							h_dQdx_xDrift_wires[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
							if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
								h_dQdx_tDriftL_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
							}

							if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
								h_dQdx_tDriftR_wires[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

		for(int i = 1; i <= maxWireGroup; i++){

			h_dQdx_xDrift_wires[i-1]->Write();
			h_dQdx_tDriftL_wires[i-1]->Write();
			h_dQdx_tDriftR_wires[i-1]->Write();

		}

		h_stats->Write();

		f.Close();

	}

	if(codeConfig == 2){

		TH2D** h_dQdx_xDrift_hits = new TH2D*[maxWireGroup];
		TH2D** h_dQdx_tDriftL_hits = new TH2D*[maxWireGroup];
		TH2D** h_dQdx_tDriftR_hits = new TH2D*[maxWireGroup];

		for(int i = 1; i <= maxWireGroup; i++){
			h_dQdx_xDrift_hits[i-1] = new TH2D(TString::Format("h_dQdx_xDrift_%dhits", i),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftL_hits[i-1] = new TH2D(TString::Format("h_dQdx_tDriftL_%dhits", i),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftR_hits[i-1] = new TH2D(TString::Format("h_dQdx_tDriftR_%dhits", i),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
		}

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			//Track counter
			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);

			auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
			int start_index = std::distance(read_x.begin(), it_start_x);

			if(start_index == read_x.GetSize()){continue;}
			auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
			int end_index = std::distance(read_x.begin(), it_end_x);
			if(end_index < 0){continue;}

			std::valarray<double> dQdx_sum(0.,maxWireGroup);
			std::valarray<double> x_sum(0.,maxWireGroup);
			std::valarray<double> t_sum(0.,maxWireGroup);
			std::valarray<int> count(0,maxWireGroup);

			
			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

					for (int N = 1; N <= maxWireGroup; N++){
					
						int j = N-1;

						h_dQdx_xDrift_hits[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
										
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

						h_dQdx_xDrift_hits[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

							h_dQdx_xDrift_hits[j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
							if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
								h_dQdx_tDriftL_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
							}

							if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
								h_dQdx_tDriftR_hits[j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

		for(int i = 1; i <= maxWireGroup; i++){

			h_dQdx_xDrift_hits[i-1]->Write();
			h_dQdx_tDriftL_hits[i-1]->Write();
			h_dQdx_tDriftR_hits[i-1]->Write();

		}

		h_stats->Write();

		f.Close();

	}

	if(codeConfig == 3){

		double** angBinLimits = getAngBinLimits(angBins);

		TH2D** h_dQdx_xDrift_ang = new TH2D*[angBins];
		TH2D** h_dQdx_tDriftL_ang = new TH2D*[angBins];
		TH2D** h_dQdx_tDriftR_ang = new TH2D*[angBins];

		for(int i = 0; i < angBins/2; i++){

			h_dQdx_xDrift_ang[i] = new TH2D(TString::Format("h_dQdx_xDrift_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftL_ang[i] = new TH2D(TString::Format("h_dQdx_tDriftL_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftR_ang[i] = new TH2D(TString::Format("h_dQdx_tDriftR_ang%ito%i", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1]),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

			h_dQdx_xDrift_ang[i + angBins/2] = new TH2D(TString::Format("h_dQdx_xDrift_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftL_ang[i + angBins/2] = new TH2D(TString::Format("h_dQdx_tDriftL_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftR_ang[i + angBins/2] = new TH2D(TString::Format("h_dQdx_tDriftR_ang%ito%i", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1]),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

		}

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);

			auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
			int start_index = std::distance(read_x.begin(), it_start_x);

			if(start_index == read_x.GetSize()){continue;}
			auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
			int end_index = std::distance(read_x.begin(), it_end_x);
			if(end_index < 0){continue;}

			double track_azimuth = atan2((*read_zf-*read_zi),(*read_xf-*read_xi))/M_PI*180;

			if(track_azimuth < -90.){

				track_azimuth = track_azimuth + 360.;

			}

			int angle_index = angBins + 1;

			for(int i = 0; i < angBins/2; i++){

				if(angBinLimits[0][i] <= track_azimuth && track_azimuth < angBinLimits[0][i+1]){

					angle_index = i;

				}
				else if(angBinLimits[1][i] <= track_azimuth && track_azimuth < angBinLimits[1][i+1]){

					angle_index = i + angBins/2;

				}

			}

			for(int i = start_index; i <= end_index; i++){

				h_dQdx_xDrift_ang[angle_index]->Fill(read_x[i], read_dqdx[i]);
				
				if(-200. < read_x[i] && read_x[i] < 0.){
					h_dQdx_tDriftL_ang[angle_index]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

				if(0. < read_x[i] && read_x[i] < 200.){
					h_dQdx_tDriftR_ang[angle_index]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

			}

		}

		for(int i = 0; i < angBins; i++){

			h_dQdx_xDrift_ang[i]->Write();
			h_dQdx_tDriftL_ang[i]->Write();
			h_dQdx_tDriftR_ang[i]->Write();

		}

		h_stats->Write();

		f.Close();

	}

	if(codeConfig == 4){

		double** angBinLimits = getAngBinLimits(angBins);

		TH2D *h_dQdx_xDrift_angWire[angBins][maxWireGroup];
		TH2D *h_dQdx_tDriftL_angWire[angBins][maxWireGroup];
		TH2D *h_dQdx_tDriftR_angWire[angBins][maxWireGroup];


		for(int i = 0; i < angBins/2; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				h_dQdx_xDrift_angWire[i][j-1] = new TH2D(TString::Format("h_dQdx_xDrift_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
				h_dQdx_tDriftL_angWire[i][j-1] = new TH2D(TString::Format("h_dQdx_tDriftL_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
				h_dQdx_tDriftR_angWire[i][j-1] = new TH2D(TString::Format("h_dQdx_tDriftR_ang%ito%inWires%d", (int)angBinLimits[0][i], (int)angBinLimits[0][i+1],j),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

				h_dQdx_xDrift_angWire[i + angBins/2][j-1] = new TH2D(TString::Format("h_dQdx_xDrift_ang%dto%dnWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j),"dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
				h_dQdx_tDriftL_angWire[i + angBins/2][j-1] = new TH2D(TString::Format("h_dQdx_tDriftL_ang%dto%dnWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
				h_dQdx_tDriftR_angWire[i + angBins/2][j-1] = new TH2D(TString::Format("h_dQdx_tDriftR_ang%dto%dnWires%d", (int)angBinLimits[1][i], (int)angBinLimits[1][i+1],j),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

			}

		}

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);

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

			double track_azimuth = atan2((*read_zf-*read_zi),(*read_xf-*read_xi))/M_PI*180;

			if(track_azimuth < -90.){

				track_azimuth = track_azimuth + 360.;

			}

			int angle_index = angBins + 1;

			for(int i = 0; i < angBins/2; i++){

				if(angBinLimits[0][i] <= track_azimuth && track_azimuth < angBinLimits[0][i+1]){

					angle_index = i;

				}
				else if(angBinLimits[1][i] <= track_azimuth && track_azimuth < angBinLimits[1][i+1]){

					angle_index = i + angBins/2;

				}

			}

			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

					for (int N = 1; N <= maxWireGroup; N++){
					
						int j = N-1;

						h_dQdx_xDrift_angWire[angle_index][j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
										
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

						h_dQdx_xDrift_angWire[angle_index][j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
						if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
							h_dQdx_tDriftL_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
						}

						if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
							h_dQdx_tDriftR_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

							h_dQdx_xDrift_angWire[angle_index][j]->Fill(x_sum[j]/count[j], dQdx_sum[j]/count[j]);
							
							if(-200. < x_sum[j]/count[j] && x_sum[j]/count[j] < 0.){
								h_dQdx_tDriftL_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
							}

							if(0. < x_sum[j]/count[j] && x_sum[j]/count[j] < 200.){
								h_dQdx_tDriftR_angWire[angle_index][j]->Fill(t_sum[j]/(count[j]*2000) - 0.2 - *read_t0/1000000, dQdx_sum[j]/count[j]);
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

		for(int i = 0; i < angBins; i++){

			for(int j = 1; j <= maxWireGroup; j++){

				h_dQdx_xDrift_angWire[i][j-1]->Write();
				h_dQdx_tDriftL_angWire[i][j-1]->Write();
				h_dQdx_tDriftR_angWire[i][j-1]->Write();

			}

		}

		h_stats->Write();

		f.Close();

	}

	if(codeConfig == 5){

		TH2D** h_dQdx_xDriftL_oct = new TH2D*[4];
		TH2D** h_dQdx_xDriftR_oct = new TH2D*[4];
		TH2D** h_dQdx_tDriftL_oct = new TH2D*[4];
		TH2D** h_dQdx_tDriftR_oct = new TH2D*[4];

		for(int i = 1; i <= 4; i++){
			
			h_dQdx_xDriftL_oct[i-1] = new TH2D(TString::Format("h_dQdx_xDriftL_oct%i", 2*i),"dQ/dx vs x", NBinsHalfX, minX, halfX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_xDriftR_oct[i-1] = new TH2D(TString::Format("h_dQdx_xDriftR_oct%i", (2*i)-1),"dQ/dx vs x", NBinsHalfX, halfX, maxX, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftL_oct[i-1] = new TH2D(TString::Format("h_dQdx_tDriftL_oct%i", 2*i),"dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
			h_dQdx_tDriftR_oct[i-1] = new TH2D(TString::Format("h_dQdx_tDriftR_oct%i", (2*i)-1),"dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

		}

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1 || (choppyLim > 0 && *read_event >= choppyLim) ){
				continue;
			}

			if(track_count % 100 == 0){std::cout << "track_count: " << track_count << std::endl;}
			track_count ++;	
			h_stats->Fill(1.);
			
			auto it_start_x = std::find_if(read_x.begin(), read_x.end(), [](float f){return !std::isnan(f);});
			int start_index = std::distance(read_x.begin(), it_start_x);

			if(start_index == read_x.GetSize()){continue;}
			auto it_end_x = std::find_if(it_start_x, read_x.end(), [](float f){return std::isnan(f);}) - 1;
			int end_index = std::distance(read_x.begin(), it_end_x);
			if(end_index < 0){continue;}

			for(int i = start_index; i <= end_index; i++){

				//Octant 1, 0 < x < 200, 0 < y < 200, 250 < z < 500
				if(0. < read_x[i] && read_x[i] < 200. && 0. < read_y[i] && read_y[i] < 200. && 250. < read_z[i] && read_z[i] < 500.){
					
					h_dQdx_xDriftR_oct[0]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftR_oct[0]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 2, -200 < x < 0, 0 < y < 200, 250 < z < 500
				if(-200. < read_x[i] && read_x[i] < 0. && 0. < read_y[i] && read_y[i] < 200. && 250. < read_z[i] && read_z[i] < 500.){
					
					h_dQdx_xDriftL_oct[0]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftL_oct[0]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 3, 0 < x < 200, 0 < y < 200, 0 < z < 250
				if(0. < read_x[i] && read_x[i] < 200. && 0. < read_y[i] && read_y[i] < 200. && 0. < read_z[i] && read_z[i] < 250.){
					
					h_dQdx_xDriftR_oct[1]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftR_oct[1]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 4, -200 < x < 0, 0 < y < 200, 0 < z < 250
				if(-200. < read_x[i] && read_x[i] < 0. && 0. < read_y[i] && read_y[i] < 200. && 0. < read_z[i] && read_z[i] < 250.){
					
					h_dQdx_xDriftL_oct[1]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftL_oct[1]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 5, 0 < x < 200, -200 < y < -0, 250 < z < 500
				if(0. < read_x[i] && read_x[i] < 200. && -200. < read_y[i] && read_y[i] < 0. && 250. < read_z[i] && read_z[i] < 500.){
					
					h_dQdx_xDriftR_oct[2]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftR_oct[2]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 6, -200 < x < 0, -200 < y < 0, 250 < z < 500
				if(-200. < read_x[i] && read_x[i] < 0. && -200. < read_y[i] && read_y[i] < 0. && 250. < read_z[i] && read_z[i] < 500.){
					
					h_dQdx_xDriftL_oct[2]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftL_oct[2]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 7, 0 < x < 200, -200 < y < 0, 0 < z < 250
				if(0. < read_x[i] && read_x[i] < 200. && -200. < read_y[i] && read_y[i] < 0. && 0. < read_z[i] && read_z[i] < 250.){
					
					h_dQdx_xDriftR_oct[3]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftR_oct[3]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

				//Octant 8, -200 < x < 0, -200 < y < 0, 0 < z < 250
				if(-200. < read_x[i] && read_x[i] < 0. && -200. < read_y[i] && read_y[i] < 0. && 0. < read_z[i] && read_z[i] < 250.){
					
					h_dQdx_xDriftL_oct[3]->Fill(read_x[i], read_dqdx[i]);
					h_dQdx_tDriftL_oct[3]->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);

				}

			}

		}

		for(int i = 1; i <= 4; i++){

			h_dQdx_xDriftL_oct[i-1]->Write();
			h_dQdx_xDriftR_oct[i-1]->Write();
			h_dQdx_tDriftL_oct[i-1]->Write();
			h_dQdx_tDriftR_oct[i-1]->Write();

		}

		h_stats->Write();

		f.Close();

	}

	return 0;

}