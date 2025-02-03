/* geoFunctions.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for timeFunctions.cpp
 * This file contains the declarations of the geodetic-utility functions for GNSS
 */

#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cmath>

// Constants
constexpr double earthRadius = 6.356752314245179e6;
constexpr double moonRadius = 1737e3;
constexpr double inclineMoon2Earth = 21.9 * M_PI / 180.0;

// Convert ECI to ECEF (of the Earth)
Eigen::Vector3d eci2ecef(const Eigen::Vector3d &eci, const std::tm &time);

// Compute Greenwich Sidereal Time (GST) in radians
double getGST(const std::tm &time);

// Convert moon LLA to ECEF (of the Earth)
std::vector<double> moonLla2ecef(const std::tm &time, double lat, double lon, double h, std::vector<double> moonPosECI);

// Convert moon LLA to MCMF (of the moon)
Eigen::Vector3d moonLla2MCMF(double lat, double lon, double h);

// Converts MCMF to ECI (of the Earth)
Eigen::Vector3d mcmf2eci(const Eigen::Vector3d &posMCMF, const Eigen::Vector3d &translation, const Eigen::Matrix3d &rotation);

// Compute the rotation matrix from MCMF to ECI
Eigen::Matrix3d getRmcmf2eci(const Eigen::Vector3d &moonPosECI);
