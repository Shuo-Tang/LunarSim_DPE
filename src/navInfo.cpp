/* navInfo.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for NavInfo class
 * The class helps to load and store the navigation information from a RINEX file
 */


#include "navInfo.h"

// NavInfo
// constructor
NavInfo::NavInfo() {
    this->numOfSatellites = 0;
}

// destructor
NavInfo::~NavInfo() = default;

// read navigation info from file
void NavInfo::readNavInfo(SettingsLoader& navSettings) {
    if (navSettings.enableGPS) {
        RinexReader gpsRinexReader(navSettings.gpsFilePath);
        gpsRinexReader.EphemerisParser(satEphemeris);
    }
    if (navSettings.enableGLONASS) {
        RinexReader gloRinexReader(navSettings.glonassFilePath);
        gloRinexReader.EphemerisParser(satEphemeris);
    }
    if (navSettings.enableGalileo) {
        RinexReader galRinexReader(navSettings.galileoFilePath);
        galRinexReader.EphemerisParser(satEphemeris);
    }
    if (navSettings.enableBeiDou) {
        RinexReader bdsRinexReader(navSettings.beidouFilePath);
        bdsRinexReader.EphemerisParser(satEphemeris);
    }

    numOfSatellites = satEphemeris.size();

}




