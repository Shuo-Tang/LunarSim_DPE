/* settingsLoader.h
*  Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for settingsLoader.cpp
 * This file contains the necessary includes and declarations for the settingsLoader.cpp file
 * The class helps to load the settings from a YAML file
 */

#pragma once
#include <iostream>
#include <string>
#include <yaml-cpp/yaml.h>
#include "timeFunctions.h"

struct LLA {
    double lat;
    double lon;
    double alt;
};

class SettingsLoader {
  public:
    SettingsLoader();
    ~SettingsLoader();

    void loadSettings(const std::string& fileName);

    std::string satPosVelPath;
    std::string satEphPath;
    std::string moonPosVelPath;

    bool enableGPS;
    bool enableGLONASS;
    bool enableGalileo;
    bool enableBeiDou;
    std::string gpsFilePath;
    std::string glonassFilePath;
    std::string galileoFilePath;
    std::string beidouFilePath;
    std::string gpsSignalPath;
    std::string gloSignalPath;
    std::string galSignalPath;
    std::vector<char> enabledConstellations;

    std::string timeEpoch;
    std::tm utcTimeEpoch;
    int timeDuration;
    int satUpdateFrequency;
    double satUpdateRate;

    LLA receiverPosLLA{};

    bool enableLunarSatellites;
    std::vector<LLA> lunarSatPosLLA;

    double mainLobeAngle;
    double sideLobeAngle;

    bool enableSignalOutput;
    double samplingFreq;

    bool enable2SP;
    bool enableDPE;
    int coherentIntegrationTime;


};
