/* ephemeris.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for Ephemeris class
 * The class helps to store the ephemeris information of a satellite
 */

#include <iostream>
#include <unordered_map>
#include "ephemeris.h"


Ephemeris::Ephemeris() {
    constellation = ' ';
    satID = "";
    gtime.time = 0;
    gtime.sec = 0;
    gpsTime = 0;
    ephemerisData = std::vector<float>(MAX_NUM_EPHEMERIS_FIELDS, std::numeric_limits<float>::quiet_NaN());
}

bool Ephemeris::operator<(const Ephemeris &ephemeris) const {
    if (gpsTime == ephemeris.gpsTime) {
        return satID[0] == ephemeris.satID[0] ? satID < ephemeris.satID : orderOfConstellations.at(constellation) < orderOfConstellations.at(ephemeris.constellation);
    }
    return gpsTime < ephemeris.gpsTime;
}


