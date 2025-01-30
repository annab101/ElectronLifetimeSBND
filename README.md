# Lifetime

Code for measuring, analysing and plotting the electron lifetime from calibration ntuples. 

## Structure

* Main - contains the macros to measure the lifetime
* Utilities - config reader script & ROOT includes
* Helpers - useful functions & constants
* Plot - various scripts for plotting various things (some more useful than others)
* Analysis - old scripts & analysis stuff (mostly bad and hacky)

## Setup

To run the macros in Main, set up sbndcode and then compile as 

```bash
g++ -o macro macro.cpp $(root-config --cflags --glibs) -lMinuit
```

Run the macro as

```bash
./macro --config ../Config/YourConfigFileName.txt
```

Macros in Plot are similar but probably have more command line inputs than just the config file.

## Running

Run lifetime analysis as 
1. dQdx_hist.cpp
* Input: calib ntuples
* Output: root file with 4 TH2Ds: dQ/dx vs x, dQ/dx vs t (East), dQ/dx vs t (West), stats (number of anode-cathode-crossing cosmics)
2. MPV.cpp
* Input: root file with dQ/dx histograms
* Output: root file with 3 TGraphAsymmErrors: MPVs vs x, MPVs vs t (East), MPVs vs t (West). Option to also save y projection and LG fit for each bin.
3. fitMPV.cpp (fits an exponential to the MPVs)
* Input: root file with MPVs
* Output: root file with 4 TF1s: expo fit for x (East, West) and t (East, West). Also, text file with lifetime values and errors. **Note: if fit status isn't "SUCCESSFUL" these values are saved as -1000**

## The config file

The config file structure & script is based heavily on similar code from Rhiannon Jones (hero!).

For each lifetime measurement, add a config text file to the 'Config' folder. Inputs needed (if there's a default, it's optional) are:

| Input    | Description |
| -------- | ------- |
| inputData  | location of a text/list file that lists the calib ntuple files    |
| dataset    | a name for the dataset you're analysing - the output folder will have this name    |
| saveLoc    | a directory where you want all the outputs to go    |
| tag | another label (e.g if you do more than one analysis per dataset, they're saved in the same folder but distinguishable by this tag)    |
| nGroupedWires    | default = 1    |
| NBinsX    | default = 100    |
| NBinsT    | default = 100    |
| minX    | default = -200.    |
| maxX    | default = 200.    |
| NBinsdQdx    | default = 75    |
| mindQdx    | default = 200.    |
| maxdQdx    | default = 1800.    |
| minT    | default = 0.   |
| maxT    | default = 1.3    |
| expo_xDriftE_lb    | lower bound for the lifetime fit in distance in the East TPC    |
| expo_xDriftW_lb    | lower bound for the lifetime fit in distance in the West TPC    |
| expo_xDriftE_ub    | upper bound for the lifetime fit in distance in the East TPC    |
| expo_xDriftW_ub    | upper bound for the lifetime fit in distance in the West TPC    |
| expo_tDriftE_lb    | lower bound for the lifetime fit in time in the East TPC    |
| expo_tDriftW_lb    | lower bound for the lifetime fit in time in the West TPC    |
| expo_tDriftE_ub    | upper bound for the lifetime fit in time in the East TPC    |
| expo_tDriftW_ub    | upper bound for the lifetime fit in time in the West TPC    |
| writeProjY    | Whether to save the y projection & LG fit in each bin (default = 0)    |