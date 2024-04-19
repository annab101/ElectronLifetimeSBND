#include "StatsHelpers.h"

namespace calib {

    double pointError(TH1D *proj_y, fitLGParameters fitParams){

        int Nslice = proj_y->GetEntries();
        return sqrt(pow(fitParams.efp[1],2)+pow((fitParams.fp[1]/sqrt(Nslice)),2));

    }

}
