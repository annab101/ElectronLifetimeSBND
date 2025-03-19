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
	
	int choppyLim = -1;
	int nGroupedWires = 1;
	int NBinsX = 100;
	int NBinsT = 100;
	int NBinsdQdx = 75;
	double minX = -200.;
	double halfX = 0.;
	double maxX = 200.;
	double mindQdx = 200.;
	double maxdQdx = 1800.;
	double minT = 0.;
	double maxT = 1.3;

	p->getValue("inputData", inputData);
	p->getValue("tag", tag);
	p->getValue("dataset", dataset);
	p->getValue("saveLoc", saveLoc);
	p->getValue("choppyLim", choppyLim);
	p->getValue("nGroupedWires", nGroupedWires);
	p->getValue("NBinsX", NBinsX);
	p->getValue("NBinsT", NBinsT);
	p->getValue("NBinsdQdx", NBinsdQdx);
	p->getValue("minX", minX);
	p->getValue("halfX", halfX);
	p->getValue("maxX", maxX);
	p->getValue("mindQdx", mindQdx);
	p->getValue("maxdQdx", maxdQdx);
	p->getValue("minT", minT);
	p->getValue("maxT", maxT);

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

	std::filesystem::create_directory((saveLoc + dataset).c_str());

	TFile f((saveLoc + dataset + "/dQdx_hist_" + std::to_string(nGroupedWires) + "wires_" + dataset + "_" + tag + ".root").c_str(), "new");

	int track_count = 0;
	int notTrack_count = 0;
	TH1D* h_stats = new TH1D("h_stats", "Number of tracks", 9, 0, 9);

	TH2D* h_dQdx_xDrift = new TH2D(TString::Format("h_dQdx_xDrift_%dwires", nGroupedWires), "dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
	TH2D* h_dQdx_tDriftE = new TH2D(TString::Format("h_dQdx_tDriftE_%dwires", nGroupedWires), "dQ/dx vs t East TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
	TH2D* h_dQdx_tDriftW = new TH2D(TString::Format("h_dQdx_tDriftW_%dwires", nGroupedWires), "dQ/dx vs t West TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

	if(nGroupedWires == 1){

		treereader.Restart();
		while(treereader.Next()){

			
			if(!((*read_xi < -195. && *read_xf > -5.) || (*read_xf < -195. && *read_xi > -5.) || (*read_xi < 5. && *read_xf > 195.) || (*read_xf < 5. && *read_xi > 195.)) || (choppyLim > 0 && *read_event >= choppyLim) ){
				notTrack_count++;
				if(notTrack_count % 100 == 0){std::cout << "notTrack_count: " << notTrack_count << std::endl;}
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

				h_dQdx_xDrift->Fill(read_x[i], read_dqdx[i]);
				
										
				if(-200. < read_x[i] && read_x[i] < 0.){
					h_dQdx_tDriftE->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

				if(0. < read_x[i] && read_x[i] < 200.){
					h_dQdx_tDriftW->Fill(read_T[i]/2000 - 0.2 - *read_t0/1000000, read_dqdx[i]);
				}

			}

		}

		h_dQdx_xDrift->Write();
		h_dQdx_tDriftE->Write();
		h_dQdx_tDriftW->Write();
		h_stats->Write();

		f.Close();
				
	}

	else {

		treereader.Restart();
		while(treereader.Next()){

			if(!((*read_xi < -195. && *read_xf > -5.) || (*read_xf < -195. && *read_xi > -5.) || (*read_xi < 5. && *read_xf > 195.) || (*read_xf < 5. && *read_xi > 195.)) || (choppyLim > 0 && *read_event >= choppyLim) ){
				notTrack_count++;
				if(notTrack_count % 100 == 0){std::cout << "notTrack_count: " << notTrack_count << std::endl;}
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

			double dQdx_sum=0;
			double x_sum=0;
			double t_sum=0;
			int count=0;
			
			for(int i = start_index; i <= end_index; i++){

				dQdx_sum += read_dqdx[i];
				x_sum += read_x[i];
				t_sum += read_T[i];
				count += 1;

				//Fill hist if swap TPC
				if( i < end_index && (read_x[i] * read_x[i+1]) < 0.){

					h_dQdx_xDrift->Fill(x_sum/count, dQdx_sum/count);
										
					if(-200. < x_sum/count && x_sum/count < 0.){
						h_dQdx_tDriftE->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
					}

					if(0. < x_sum/count && x_sum/count < 200.){
						h_dQdx_tDriftW->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
					}				

				dQdx_sum = 0.;
				x_sum = 0.;
				t_sum = 0.;
				count = 0;

				}
				else if(i == end_index){

					h_dQdx_xDrift->Fill(x_sum/count, dQdx_sum/count);
						
					if(-200. < x_sum/count && x_sum/count < 0.){
						h_dQdx_tDriftE->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
					}

					if(0. < x_sum/count && x_sum/count < 200.){
						h_dQdx_tDriftW->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
					}

					dQdx_sum = 0.;
					x_sum = 0.;
					t_sum = 0.;
					count = 0;

				}
				else{
						
					if((read_wire[i] - minWire)/nGroupedWires != (read_wire[i+1] - minWire)/nGroupedWires){

						h_dQdx_xDrift->Fill(x_sum/count, dQdx_sum/count);
						
						if(-200. < x_sum/count && x_sum/count < 0.){
							h_dQdx_tDriftE->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
						}

						if(0. < x_sum/count && x_sum/count < 200.){
							h_dQdx_tDriftW->Fill(t_sum/(count*2000) - 0.2 - *read_t0/1000000, dQdx_sum/count);
						}

						dQdx_sum = 0.;
						x_sum = 0.;
						t_sum = 0.;
						count = 0;

					}

				}

			}

		}

		h_dQdx_xDrift->Write();
		h_dQdx_tDriftE->Write();
		h_dQdx_tDriftW->Write();

		h_stats->Write();

		f.Close();

	}

	return 0;

}