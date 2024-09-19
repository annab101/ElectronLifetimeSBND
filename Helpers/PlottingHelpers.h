/*
* Helper functions for plotting
*/

#ifndef PLOTTINGHELPERS_H
#define PLOTTINGHELPERS_H

namespace calib {

    /**
    * @brief Set marker style.
    *
    * @param t The graph/histogram
    * @param color Marker colour
    * @param style Marker style
    * @param size Marker size
    *
    */

    template<class T> void setMarker(T *t, Color_t color, Style_t style, Size_t size);
    
    /**
    * @brief Save canvas in png, pdf and root format.
    *
    * @param c The canvas
    * @param plotName Name of plot
    *
    */

    void saveFig(TCanvas *c, std::string plotName);

    /**
    * @brief Open existing canvas and clear it.
    *
    * @param c The canvas
    *
    */

    void openAndClear(TCanvas *c);

    /**
    * @brief Draw lifetime stats box for two functions.
    *
    * @param pos Stats box position
    * @param track_count Number of AC tracks
    * @param f1 The first function
    * @param f2 The second function
    *
    */

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1, TF1* f2);
    
    /**
    * @brief Draw lifetime stats box for one function.
    *
    * @param pos Stats box position
    * @param track_count Number of AC tracks
    * @param f1 The function
    *
    */

    TPaveText *statsBox(std::vector<double> pos, int track_count, TF1* f1);
    
    /**
    * @brief Draw stats box with number of AC crossers.
    *
    * @param pos Stats box position
    * @param track_count Number of AC tracks
    *
    */

    TPaveText *statsBox(std::vector<double> pos, int track_count);

    /**
    * @brief Put labels on plot (e.g SBND Preliminary,).
    *
    * @param pos Stats box position
    * @param lines Lines for stats box
    * @param textColor Colour of text
    * @param font Font
    * @param size Font size
    */

    TPaveText *plotLabels(std::vector<double> pos, std::vector<std::string> lines, Color_t textColor, int font, int size);

    /**
    * @brief Set font size for a plot.
    *
    * @param t The thing being plotted
    * @param font Font
    * @param size Font size
    *
    */

    template<class T> void setFontSize(T *t, int font, int size);

    /**
    * @brief Set font size for the Z axis of a plot.
    *
    * @param h The histogram
    * @param font Font
    * @param size Font size
    *
    */

    void setFontSizeZ(TH2 *h, int font, int size);
}

#endif /* PLOTTINGHELPERS_H */