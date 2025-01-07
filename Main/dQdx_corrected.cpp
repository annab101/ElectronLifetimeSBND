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
	double tE;
    double tW;

	for(int i=0; i<argc; ++i){
		if(!strcmp(argv[i], "--config")){
			filename = argv[i+1];
		}
		if(!strcmp(argv[i], "--tE")){
			tE = std::stod(argv[i+1]);
		}
        if(!strcmp(argv[i], "--tW")){
			tW = std::stoi(argv[i+1]);
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
	std::string configLabel = "noConfigLabel";

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
	TTreeReaderValue<int> read_selected(treereader, "trk.selected");
	TTreeReaderArray<uint16_t> read_wire(treereader, "trk.hits2.h.wire");

	TFile f((saveLoc + dataset  + "_" + configLabel + "/dQdxCorrected_hist_" + configLabel + "_" + dataset + "_" + tag + ".root").c_str(), "new");

	int track_count = 0;

	if(codeConfig == 0){

		TH2D* h_dQdx_xDrift_basic = new TH2D("h_dQdx_xDrift_basic", "dQ/dx vs x", NBinsX, minX, maxX, NBinsdQdx, mindQdx, maxdQdx);
		TH2D* h_dQdx_tDriftL_basic = new TH2D("h_dQdx_tDriftL_basic", "dQ/dx vs t Left TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);
		TH2D* h_dQdx_tDriftR_basic = new TH2D("h_dQdx_tDriftR_basic", "dQ/dx vs t Right TPC", NBinsT, minT, maxT, NBinsdQdx, mindQdx, maxdQdx);

		TH1D* h_stats = new TH1D("h_stats", "Number of tracks", 1, 0, 2);

		treereader.Restart();
		while(treereader.Next()){

			if(*read_selected != 1){
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
										
				if(-200. < read_x[i] && read_x[i] < 0.){
					double t = read_T[i]/2000 - 0.2 - *read_t0/1000000;
					double dQdx_x = read_dqdx[i]*TMath::Exp((200 - TMath::Abs(read_x[i]))/(vDrift*tE));
					double dQdx_T = read_dqdx[i]*TMath::Exp(t/tE);
					h_dQdx_xDrift_basic->Fill(read_x[i], dQdx_x);
					h_dQdx_tDriftL_basic->Fill(t, dQdx_T);
				}

				if(0. < read_x[i] && read_x[i] < 200.){
					double t = read_T[i]/2000 - 0.2 - *read_t0/1000000;
					double dQdx_x = read_dqdx[i]*TMath::Exp((200 - TMath::Abs(read_x[i]))/(vDrift*tW));
					double dQdx_T = read_dqdx[i]*TMath::Exp(t/tW);
					h_dQdx_xDrift_basic->Fill(read_x[i], dQdx_x);
					h_dQdx_tDriftR_basic->Fill(t, dQdx_T);
				}

			}

		}

		h_dQdx_xDrift_basic->Write();
		h_dQdx_tDriftL_basic->Write();
		h_dQdx_tDriftR_basic->Write();
		h_stats->Write();

		f.Close();
				
	}

	return 0;

}