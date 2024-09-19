/*
* Constants for use in code
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {

    //Sheffield colour scheme
	Int_t deepViolet = TColor::GetColor("#440099");
	Int_t coral = TColor::GetColor("#E7004C");
	Int_t powderBlue = TColor::GetColor("#9ADBE8");
	Int_t flamingo = TColor::GetColor("#FF6371");
	Int_t mintGreen = TColor::GetColor("#00CE7C");
	double vDrift = 156.267;

	std::map<int, std::string> configMap{{0, "classic"}, {1, "multiwire"}, {2, "multihit"}, {3, "angle"}, {4, "wireAndAngle"}, {5, "octant"}};
	std::map<int, std::string> octantName{{1, "NW Upper"}, {2, "NE Upper"}, {3, "SW Upper"}, {4, "SE Upper"}, {5, "NW Lower"}, {6, "NE Lower"}, {7, "SW Lower"}, {8, "SE Lower"}};
	std::map<int, std::string> plotFitName{{0, "noFit"}, {1, "withFit"}};

}

#endif
