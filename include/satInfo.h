/* satInfo.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for satInfo.cpp
 * This file contains the necessary includes and declarations for the satInfo.cpp file
 * The class helps to compute the satellite information from ephemeris
 */

#pragma once

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <array>
#include <Eigen/Dense>
#include <iomanip>
#include "ephemeris.h"
#include "navInfo.h"
#include "settingsLoader.h"

void displayProgressBar(int current, int total, int barWidth = 50);

class SatInfo {
  public:
    SatInfo();

    // load satellite ephemeris according to desired epoch (batch)
    void loadSatEph(const NavInfo& navInfo, const SettingsLoader& navSettings);

    // generate satellite position and velocity (batch) in the time duration
    void generateSatPosVel();

    // find the closest ephemeris to the desired epoch
    void findClosestEph(const std::vector<Ephemeris>& ephemeris);

    // Export satellite positions and velocities
    void exportPosVel(const std::string& filename) const;

    // Export selected ephemeris
    void exportEphemeris(const std::string& filename) const;

    int numOfSatellites;

    // desired epoch in GPS time
    double atTimeEpoch;
    // time duration
    int timeDuration{};
    // frequency of the satellite information update (position and velocity)
    int satUpdateFrequency;
    // rate of the satellite information update
    double satUpdateRate;
    // number of satellite information update epochs
    int numOfSatUpdates;
    // selected satellites at the desired epoch
    std::vector<Ephemeris> selectedSatellites;

    std::vector<std::vector<double>> satPosVel;
};

