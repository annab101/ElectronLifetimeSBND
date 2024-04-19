/*
* Helper functions for fitting
*/

#ifndef STATSHELPERS_H
#define STATSHELPERS_H

namespace calib {

    /**
    * @brief Statistical error on an MPV.
    *
    * @param proj_y The histogram MPV is from
    * @param fitParams Fit params from MPV fit
    *
    */

    double pointError(TH1D *proj_y, fitLGParameters fitParams);

}

#endif
