/* timeFunctions.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for timeFunctions.h
 * This file contains the necessary functions for geometry computation
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include "orbitFunctions.h"

std::pair<double, double> eccAnomaly(const double& time, const Ephemeris& eph) {
    double GM;

    switch (eph.constellation) {
     case 'G':
      GM = GM_GPS;
     break;
     case 'R':
      GM = GM_GLO;
     break;
     case 'E':
      GM = GM_GAL;
     break;
     case 'C':
      GM = GM_BDS;
     break;
     default: {
      std::cerr << "ERROR: Invalid constellation!" << std::endl;
      return {0.0, 0.0};
     }
    }

    double M0 = eph.ephemerisData[6];
    double sqrtA = eph.ephemerisData[10];
    double deltaN = eph.ephemerisData[5];
    double ecc = eph.ephemerisData[8];
    double toe = eph.ephemerisData[11];
    int numOfWeek = static_cast<int>(eph.ephemerisData[21]);
    double timeEph = weekTow2gps(numOfWeek, toe, eph.constellation);

    double A = sqrtA * sqrtA;
    double tk = checkWeekCross(time - timeEph);
    double n0 = sqrt(GM / pow(A, 3));
    double n = n0 + deltaN;
    double Mk = M0 + n * tk;
    Mk = fmod(Mk + CIRCLE_RAD, CIRCLE_RAD);

    double Ek = Mk;
    constexpr int maxIter = 10;

    for (int i = 0; i < maxIter; ++i) {
        double Ek_old = Ek;
        Ek = Mk + ecc * sin(Ek);
        double dEk = fmod(Ek - Ek_old, CIRCLE_RAD);
        if (std::abs(dEk) < 1e-12) {
            break;
        }
    }

    Ek = fmod(Ek + CIRCLE_RAD, CIRCLE_RAD);
    return {Ek, n};
}

std::vector<double> satelliteMotionDiff(
    const std::vector<double>& pos,
    const std::vector<double>& vel,
    const std::vector<double>& acc,
    double a,
    double GM,
    double J2,
    double omegaeDot
) {
    std::vector<double> diff(6); // pos and vel differential

    // Rename variables for clarity
    double x = pos[0], y = pos[1], z = pos[2];
    double vx = vel[0], vy = vel[1];
    double ax = acc[0], ay = acc[1], az = acc[2];

    // Compute parameters
    double r = std::sqrt(x * x + y * y + z * z);
    double g = -GM / std::pow(r, 3);
    double h = J2 * 1.5 * std::pow(a / r, 2);
    double k = 5 * std::pow(z / r, 2);

    // Differential position
    diff[0] = vx;
    diff[1] = vy;
    diff[2] = vel[2];

    // Differential velocity
    diff[3] = g * x * (1 - h * (k - 1)) + ax + omegaeDot * omegaeDot * x + 2 * omegaeDot * vy;
    diff[4] = g * y * (1 - h * (k - 1)) + ay + omegaeDot * omegaeDot * y - 2 * omegaeDot * vx;
    diff[5] = g * z * (1 - h * (k - 3)) + az;

    return diff;
}
