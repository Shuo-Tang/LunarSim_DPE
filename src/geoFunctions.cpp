/* geoFunctions.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for timeFunctions.h
 * This file contains the declarations of the geodetic-utility functions for GNSS
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include "geoFunctions.h"

// Convert ECI to ECEF (of the Earth)
Eigen::Vector3d eci2ecef(const Eigen::Vector3d &eci, const std::tm &time) {
    // Compute Greenwich Sidereal Time (GST)
    double gst = getGST(time);

    // Define rotation matrix from GST
    Eigen::Matrix3d rotMatrix;
    rotMatrix << std::cos(gst), -std::sin(gst), 0,
                 std::sin(gst),  std::cos(gst), 0,
                          0,             0,     1;

    // Convert to ECI by multiplying the rotation matrix with the ECEF vector
    return rotMatrix * eci;
}

// Compute Greenwich Sidereal Time (GST) in radians
double getGST(const std::tm &time) {
    // Julian Date calculation
    int year = time.tm_year + 1900;
    int month = time.tm_mon + 1;
    if (month <= 2) {
        year--;
        month += 12;
    }
    int day = time.tm_mday;
    double hour = time.tm_hour + time.tm_min / 60.0 + time.tm_sec / 3600.0;

    int A = year / 100;
    int B = 2 - A + A / 4;
    double JD = static_cast<int>(365.25 * (year + 4716)) + static_cast<int>(30.6001 * (month + 1)) + day + hour / 24.0 + B - 1524.5;

    // Calculate GST
    double T = (JD - 2451545.0) / 36525.0; // Julian centuries since J2000.0
    double GST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + T * T * (0.000387933 - T / 38710000.0);
    GST = std::fmod(GST, 360.0); // Restrict to [0, 360)
    if (GST < 0) GST += 360.0;

    return GST * M_PI / 180.0; // Convert to radians
}


// Convert moon LLA to ECEF (of the Earth)
std::vector<double> moonLla2ecef(const std::tm &time, double lat, double lon, double h, std::vector<double> moonPosECI) {
    Eigen::Vector3d posMCMF = moonLla2MCMF(lat, lon, h);
    auto translation = Eigen::Vector3d(moonPosECI[0], moonPosECI[1], moonPosECI[2]);
    auto rotation = getRmcmf2eci(translation);
    Eigen::Vector3d posECI = mcmf2eci(posMCMF, translation, rotation);
    Eigen::Vector3d posECEF = eci2ecef(posECI, time);
    return {posECEF[0], posECEF[1], posECEF[2]};
}

// Convert moon LLA to MCMF (of the moon)
Eigen::Vector3d moonLla2MCMF(double lat, double lon, double h) {
    // Convert degrees to radians
    double latRad = lat * M_PI / 180.0;
    double lonRad = lon * M_PI / 180.0;

    // Compute Cartesian coordinates
    double r = moonRadius + h;
    double X = r * std::cos(lonRad) * std::cos(latRad);
    double Y = r * std::cos(lonRad) * std::sin(latRad);
    double Z = r * std::sin(lonRad);

    return {X, Y, Z};
}

// Converts MCMF to ECI (of the Earth)
Eigen::Vector3d mcmf2eci(const Eigen::Vector3d &posMCMF, const Eigen::Vector3d &translation, const Eigen::Matrix3d &rotation) {
    // Apply rotation and translation
    return translation + rotation * posMCMF;
}

Eigen::Matrix3d getRmcmf2eci(const Eigen::Vector3d &moonPosECI) {
    // Normalize MCMF X-axis
    Eigen::Vector3d xMCMF = - moonPosECI.normalized();

    // MCMF Z-axis: Moon's rotation axis inclined by 21.9° to Earth's equator
    Eigen::Vector3d zMCMF(0, std::sin(inclineMoon2Earth), std::cos(inclineMoon2Earth));

    // MCMF Y-axis: Right-hand rule
    Eigen::Vector3d yMCMF = zMCMF.cross(xMCMF).normalized();

    // Recompute MCMF Z-axis to ensure orthogonality
    zMCMF = xMCMF.cross(yMCMF);

    // Construct the rotation matrix
    Eigen::Matrix3d Rmcmf2eci;
    Rmcmf2eci.col(0) = xMCMF;
    Rmcmf2eci.col(1) = yMCMF;
    Rmcmf2eci.col(2) = zMCMF;

    return Rmcmf2eci;
}