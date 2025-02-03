/* navInfo.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header filer for NavInfo class
 * This file contains the necessary includes and declarations for the navInfo.cpp file
 * The class helps to load and store the navigation information from a RINEX file
 */


#pragma once
#include <iostream>
#include <vector>
#include <string>

#include <yaml-cpp/yaml.h>
#include "rinexReader.h"
#include "settingsLoader.h"



class NavInfo {
  public:
    NavInfo();
    ~NavInfo();

    // read navigation info from file
    void readNavInfo(SettingsLoader& navSettings);

    // number of satellites
    int numOfSatellites;

    // all constellations and satellites
    // G: GPS, R: GLONASS, E: Galileo, C: BeiDou
    std::vector<Ephemeris> satEphemeris;


};