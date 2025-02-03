/* geoFunctions.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for geoFunctions.cpp
 * This file contains the declarations of the geometry functions for GNSS
 */

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include "timeFunctions.h"
#include "ephemeris.h"

// some predefined constants for satellite orbit computation
constexpr double OMEGAE_DOT_GPS = 7.2921151467e-5;      // GPS     Angular velocity of the Earth rotation [rad/s]
constexpr double OMEGAE_DOT_GLO = 7.292115e-5;          // GLONASS Angular velocity of the Earth rotation [rad/s]
constexpr double OMEGAE_DOT_GAL = 7.2921151467e-5;      // Galileo Angular velocity of the Earth rotation [rad/s]
constexpr double OMEGAE_DOT_BDS = 7.292115e-5;          // BeiDou  Angular velocity of the Earth rotation [rad/s]

constexpr double GM_GPS = 3.986005e14;                  // Gravitational constant for GPS
constexpr double GM_GLO = 3.9860044e14;                 // Gravitational constant for GLONASS
constexpr double GM_GAL = 3.986004418e14;               // Gravitational constant for Galileo
constexpr double GM_BDS = 3.986004418e14;               // Gravitational constant for BDS

constexpr int ELL_A_GPS = 6378137;                      // GPS (WGS-84)      Ellipsoid semi-major axis [m]
constexpr int ELL_A_GLO = 6378136;                      // GLONASS (PZ-90)   Ellipsoid semi-major axis [m]
constexpr int ELL_A_GAL = 6378137;                      // Galileo (GTRF)    Ellipsoid semi-major axis [m]
constexpr int ELL_A_BDS = 6378136;                      // BeiDou (CGCS2000) Ellipsoid semi-major axis [m]

constexpr double J2_GLO = 1.0826257e-3;

constexpr double CIRCLE_RAD = 2 * M_PI;         // Circle in radians



// Compute the eccentric anomaly
std::pair<double, double> eccAnomaly(const double& time, const Ephemeris& eph);

// GLONASS satellite motion differential equations
std::vector<double> satelliteMotionDiff(
    const std::vector<double>& pos,
    const std::vector<double>& vel,
    const std::vector<double>& acc,
    double a,
    double GM,
    double J2,
    double Omegae_dot = 0.0
);

