/* satInfo.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for SatInfo class
 * The class helps to compute the satellite information from ephemeris
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include "satInfo.h"
#include "ephemeris.h"
#include "orbitFunctions.h"


// comparator for sorting at the closest time epoch
struct EphemerisComparatorAtEpoch {
    bool operator()(const Ephemeris& eph1, const Ephemeris& eph2) const {
        return eph1.satID[0] == eph2.satID[0] ? eph1.satID < eph2.satID : orderOfConstellations.at(eph1.constellation) < orderOfConstellations.at(eph2.constellation);
    }
};

void displayProgressBar(int current, int total, int barWidth) {
    double progress = static_cast<double>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "% (" << current << "/" << total << ")" << std::flush;
}



SatInfo::SatInfo() {
    this->atTimeEpoch = 0;
	this->numOfSatellites = 0;
    this->timeDuration = 0;
    this->satUpdateFrequency = 0;
    this->satUpdateRate = 0;
    this->numOfSatUpdates = 0;
}

void SatInfo::loadSatEph(const NavInfo& navInfo, const SettingsLoader& navSettings) {
    gtime_t txgTime = str2gtime(navSettings.timeEpoch.substr(2));
    atTimeEpoch = gtime2gps(txgTime);
    timeDuration = navSettings.timeDuration;
    satUpdateFrequency = navSettings.satUpdateFrequency;
    satUpdateRate = navSettings.satUpdateRate;
    numOfSatUpdates = timeDuration * satUpdateFrequency;
    this->satPosVel.reserve(numOfSatUpdates * numOfSatellites);

    std::vector<Ephemeris> gpsEph;
    std::vector<Ephemeris> gloEph;
    std::vector<Ephemeris> galEph;
    std::vector<Ephemeris> bdsEph;
    // group by constellation
    for (int iSat = 0; iSat < navInfo.numOfSatellites; iSat++) {
        if (navInfo.satEphemeris[iSat].constellation == 'G') {
            gpsEph.push_back(navInfo.satEphemeris[iSat]);
        } else if (navInfo.satEphemeris[iSat].constellation == 'R') {
            gloEph.push_back(navInfo.satEphemeris[iSat]);
        } else if (navInfo.satEphemeris[iSat].constellation == 'E') {
            galEph.push_back(navInfo.satEphemeris[iSat]);
        } else if (navInfo.satEphemeris[iSat].constellation == 'C') {
            bdsEph.push_back(navInfo.satEphemeris[iSat]);
        }
    }
    // pick by the closest time epoch
    if (navSettings.enableGPS == true) {
        findClosestEph(gpsEph);
    }
    if (navSettings.enableGLONASS == true) {
        findClosestEph(gloEph);
    }
    if (navSettings.enableGalileo == true) {
        findClosestEph(galEph);
    }
    if (navSettings.enableBeiDou == true) {
        findClosestEph(bdsEph);
    }
    std::sort(selectedSatellites.begin(), selectedSatellites.end(), EphemerisComparatorAtEpoch());
    numOfSatellites = selectedSatellites.size();
}

void SatInfo::findClosestEph(const std::vector<Ephemeris> &ephemeris) {
    // binary search (since the ephemeris is sorted by time)
    auto itInit = std::lower_bound(ephemeris.begin(), ephemeris.end(), atTimeEpoch, [](const Ephemeris& eph, const double& time) {
        return eph.gpsTime < time;
    });
    std::unordered_set<std::string> satIDs;
    if (itInit != ephemeris.end()) {
        selectedSatellites.push_back(*itInit);
        satIDs.insert(itInit->satID);
    }
    float dtMax = 0.0;
    if (itInit->constellation == 'G') {
        dtMax = itInit->ephemerisData[28] * 3600 / 2;
    }
    else if (itInit->constellation == 'R') {
        dtMax = 900;
    }
    else {
        dtMax = 7200;
    }
    // check the all closest ephemeris
    if (itInit == ephemeris.begin()) {
        auto itForward = itInit;
        while (itForward != ephemeris.end()) {
            ++itForward;
            if (std::abs(itForward->gpsTime - atTimeEpoch) > dtMax) {
                break;
            }
            if (satIDs.contains(itForward->satID)) {
                continue;
            }
            selectedSatellites.push_back(*itForward);
            satIDs.insert(itForward->satID);
        }
    } else if (itInit == ephemeris.end()) {
        auto itBackward = itInit;
        while (itBackward != ephemeris.begin()) {
            --itBackward;
            if (std::abs(itBackward->gpsTime - atTimeEpoch) > dtMax) {
                break;
            }
            if (satIDs.contains(itBackward->satID)) {
                continue;
            }
            selectedSatellites.push_back(*itBackward);
            satIDs.insert(itBackward->satID);
        }
    }
    else {
        auto itBackward = itInit;
        auto itForward = itInit;
        if (std::abs((itInit - 1)->gpsTime - atTimeEpoch) <= std::abs((itInit + 1)->gpsTime - atTimeEpoch)) {
            while (itBackward != ephemeris.begin()) {
                --itBackward;
                if (std::abs(itBackward->gpsTime - atTimeEpoch) > dtMax) {
                    break;
                }
                if (satIDs.contains(itBackward->satID)) {
                    continue;
                }
                selectedSatellites.push_back(*itBackward);
                satIDs.insert(itBackward->satID);
            }
            while (itForward != ephemeris.end()) {
                ++itForward;
                if (std::abs(itForward->gpsTime - atTimeEpoch) > dtMax) {
                    break;
                }
                if (satIDs.contains(itForward->satID)) {
                    continue;
                }
                selectedSatellites.push_back(*itForward);
                satIDs.insert(itForward->satID);
            }
        }
        else {
            while (itForward != ephemeris.end()) {
                ++itForward;
                if (std::abs(itForward->gpsTime - atTimeEpoch) > dtMax) {
                    break;
                }
                if (satIDs.contains(itForward->satID)) {
                    continue;
                }
                selectedSatellites.push_back(*itForward);
                satIDs.insert(itForward->satID);
            }
            while (itBackward != ephemeris.begin()) {
                --itBackward;
                if (std::abs(itBackward->gpsTime - atTimeEpoch) > dtMax) {
                    break;
                }
                if (satIDs.contains(itBackward->satID)) {
                    continue;
                }
                selectedSatellites.push_back(*itBackward);
                satIDs.insert(itBackward->satID);
            }
        }
    }
}


void SatInfo::generateSatPosVel() {
    double t = atTimeEpoch;
    for (int iEpoch = 0; iEpoch < numOfSatUpdates; iEpoch++) {
        // Display the progress bar
        displayProgressBar(iEpoch + 1, numOfSatUpdates);
        for (int iSat = 0; iSat < numOfSatellites; iSat++) {
            const Ephemeris eph = selectedSatellites[iSat];
            char constellation = eph.constellation;
            int prn = std::stoi(eph.satID.substr(1));
            std::vector<double> posVel(6, 0.0);
            double px, py, pz;
            double vx, vy, vz;
            // correct the time for BDS
            if (constellation == 'C') {
                t -= 14.0;
            }
            // determine the angular velocity of the Earth rotation
            double omegaeDot;
            switch (constellation) {
                case 'G':
                    omegaeDot = OMEGAE_DOT_GPS;
                break;
                case 'R':
                    omegaeDot = OMEGAE_DOT_GLO;
                break;
                case 'E':
                    omegaeDot = OMEGAE_DOT_GAL;
                break;
                case 'C':
                    omegaeDot = OMEGAE_DOT_BDS;
                break;
                default:
                    throw std::runtime_error("Unrecognized satellite system!");
            }
            // for GPS, Galileo and Beidou
            if (constellation != 'R') {
                double omegak;
                // load the ephemeris data
                double sqrtA = eph.ephemerisData[10];
                double ecc = eph.ephemerisData[8];
                double omega = eph.ephemerisData[17];
                double cuc = eph.ephemerisData[7];
                double cus = eph.ephemerisData[9];
                double crc = eph.ephemerisData[16];
                double crs = eph.ephemerisData[4];
                double i0 = eph.ephemerisData[15];
                double idot = eph.ephemerisData[19];
                double cic = eph.ephemerisData[12];
                double cis = eph.ephemerisData[14];
                double omega0 = eph.ephemerisData[13];
                double omegaDot = eph.ephemerisData[18];
                double toe = eph.ephemerisData[11];
                int numOfWeek = static_cast<int>(eph.ephemerisData[21]);
                if (constellation == 'E' and numOfWeek > 2500) {
                    numOfWeek -= 1024;
                }
                double timeEph = weekTow2gps(numOfWeek, toe, constellation);

                /**************************************/
                /*** compute the satellite position ***/
                /**************************************/
                double A = sqrtA * sqrtA;  // Semi-major axis
                double tk = checkWeekCross(t - timeEph);
                auto [Ek, n] = eccAnomaly(t, eph);
                double fk = atan2(sqrt(1 - ecc * ecc) * sin(Ek), cos(Ek) - ecc); // True anomaly
                double phik = fk + omega;
                phik = fmod(phik, CIRCLE_RAD);

                double uk = phik + cuc * cos(2 * phik) + cus * sin(2 * phik); // Corrected argument of latitude
                double rk = A * (1 - ecc * cos(Ek)) + crc * cos(2 * phik) + crs * sin(2 * phik); // Corrected radial distance
                double ik = i0 + idot * tk + cic * cos(2 * phik) + cis * sin(2 * phik); // Corrected inclination

                // Satellite positions in orbital plane
                double x1k = cos(uk) * rk;
                double y1k = sin(uk) * rk;

                // GPS/Galileo/Beidou MEO satellites
                if (constellation != 'C' or (constellation == 'C' and prn > 5)) {
                    omegak = omega0 + (omegaDot - omegaeDot) * tk - omegaeDot * toe;
                    omegak = fmod(omegak + CIRCLE_RAD, CIRCLE_RAD);

                    px = x1k * cos(omegak) - y1k * cos(ik) * sin(omegak);
                    py = x1k * sin(omegak) + y1k * cos(ik) * cos(omegak);
                    pz = y1k * sin(ik);
                }
                // Beidou GEO satellites
                else {
                    // If GEO BeiDou satellite (ranging code number <= 5)
                    omegak = omega0 + omegaDot * tk - omegaeDot * toe;
                    omegak = fmod(omegak + CIRCLE_RAD, CIRCLE_RAD);

                    double xgk = x1k * cos(omegak) - y1k * cos(ik) * sin(omegak);
                    double ygk = x1k * sin(omegak) + y1k * cos(ik) * cos(omegak);
                    double zgk = y1k * sin(ik);

                    // Inertial coordinates
                    Eigen::Matrix3d Rx;
                    Rx << 1, 0, 0,
                          0, cos(-5 * M_PI / 180.0), sin(-5 * M_PI / 180.0),
                          0, -sin(-5 * M_PI / 180.0), cos(-5 * M_PI / 180.0);

                    double oedt = omegaeDot * tk;
                    Eigen::Matrix3d Rz;
                    Rz << cos(oedt), sin(oedt), 0,
                          -sin(oedt), cos(oedt), 0,
                          0, 0, 1;
                    // Apply rotation
                    Eigen::Vector3d Xgk(xgk, ygk, zgk);
                    Eigen::Vector3d Xk = Rz * Rx * Xgk;
                    px = Xk(0);
                    py = Xk(1);
                    pz = Xk(2);
                }

                /**************************************/
                /*** compute the satellite velocity ***/
                /**************************************/
                double MkDot = n;
                double EkDot = MkDot / (1 - ecc * cos(Ek));
                double fkDot = sin(Ek) * EkDot * (1 + ecc * cos(fk)) / ((1 - cos(Ek) * ecc) * sin(fk));
                double phikDot = fkDot;
                double ukDot = phikDot + 2 * (cus * cos(2 * phik) - cuc * sin(2 * phik)) * phikDot;
                double rkDot = A * ecc * sin(Ek) * EkDot + 2 * (crs * cos(2 * phik) - crc * sin(2 * phik)) * phikDot;
                double ikDot = idot + 2 * (cis * cos(2 * phik) - cic * sin(2 * phik)) * phikDot;
                double omegakDot = omegaDot - omegaeDot;
                double x1kDot = rkDot * cos(uk) - y1k * ukDot;
                double y1kDot = rkDot * sin(uk) + x1k * ukDot;

                vx = x1kDot * cos(omegak) - y1kDot * cos(ik) * sin(omegak) + y1k * sin(ik) * sin(omegak) * ikDot - py *
                            omegakDot;
                vy = x1kDot * sin(omegak) + y1kDot * cos(ik) * cos(omegak) - y1k * sin(ik) * ikDot * cos(omegak) + px *
                            omegakDot;
                vz = y1kDot * sin(ik) + y1k * cos(ik) * ikDot;
            }
            // for GLONASS
            else {
                // load the ephemeris data
                auto [weekToe, toe] = gps2weekTow(eph.gpsTime);
                double timeEph = weekTow2gps(weekToe, toe, constellation);
                double x = eph.ephemerisData[3];
                double y = eph.ephemerisData[7];
                double z = eph.ephemerisData[11];
                double dx = eph.ephemerisData[4];
                double dy = eph.ephemerisData[8];
                double dz = eph.ephemerisData[12];
                double dx2 = eph.ephemerisData[5];
                double dy2 = eph.ephemerisData[9];
                double dz2 = eph.ephemerisData[13];
                // compute the satellite position
                double int_step = 60; // Integration step in seconds
                double tk = checkWeekCross(t - timeEph);
                // confirm the number of iterations
                int n = static_cast<int>(std::floor(std::abs(tk / int_step)));
                double int_step_res = std::fmod(tk, int_step);
                if (int_step_res != 0) {
                    n++;
                }

                std::vector<double> pos = {x, y, z};
                std::vector<double> vel = {dx, dy, dz};
                std::vector<double> acc = {dx2, dy2, dz2};

                for (int s = 0; s < n; ++s) {
                    double step = (s == n - 1 && int_step_res != 0) ? int_step_res : int_step;

                    // Runge-Kutta steps
                    // step 1
                    auto pos_vel1_dot = satelliteMotionDiff(pos, vel, acc, ELL_A_GLO, GM_GLO, J2_GLO, OMEGAE_DOT_GLO);

                    std::vector<double> pos2(3), vel2(3);
                    for (int i = 0; i < 3; ++i) {
                        pos2[i] = pos[i] + pos_vel1_dot[i] * step / 2;
                        vel2[i] = vel[i] + pos_vel1_dot[i + 3] * step / 2;
                    }

                    auto pos_vel2_dot = satelliteMotionDiff(pos2, vel2, acc, ELL_A_GLO, GM_GLO, J2_GLO, OMEGAE_DOT_GLO);

                    std::vector<double> pos3(3), vel3(3);
                    for (int i = 0; i < 3; ++i) {
                        pos3[i] = pos[i] + pos_vel2_dot[i] * step / 2;
                        vel3[i] = vel[i] + pos_vel2_dot[i + 3] * step / 2;
                    }

                    auto pos_vel3_dot = satelliteMotionDiff(pos3, vel3, acc, ELL_A_GLO, GM_GLO, J2_GLO, OMEGAE_DOT_GLO);

                    std::vector<double> pos4(3), vel4(3);
                    for (int i = 0; i < 3; ++i) {
                        pos4[i] = pos[i] + pos_vel3_dot[i] * step;
                        vel4[i] = vel[i] + pos_vel3_dot[i + 3] * step;
                    }

                    auto pos_vel4_dot = satelliteMotionDiff(pos4, vel4, acc, ELL_A_GLO, GM_GLO, J2_GLO, OMEGAE_DOT_GLO);

                    for (int i = 0; i < 3; ++i) {
                        pos[i] += (pos_vel1_dot[i] + 2 * pos_vel2_dot[i] + 2 * pos_vel3_dot[i] + pos_vel4_dot[i]) * step / 6;
                        vel[i] += (pos_vel1_dot[i + 3] + 2 * pos_vel2_dot[i + 3] + 2 * pos_vel3_dot[i + 3] + pos_vel4_dot[i + 3]) * step / 6;
                    }
                }
                // Transformation from PZ-90.02 to WGS-84 (G1150)
                px = pos[0] - 0.36;
                py = pos[1] + 0.08;
                pz = pos[2] + 0.18;

                vx = vel[0];
                vy = vel[1];
                vz = vel[2];
            }
            // format the satellite position and velocity
            posVel[0] = px;
            posVel[1] = py;
            posVel[2] = pz;
            posVel[3] = vx;
            posVel[4] = vy;
            posVel[5] = vz;
            // format the time and sat ID
            satPosVel.push_back(posVel);
        }
        // update GPS time
        t += satUpdateRate;
    }
    std::cout << std::endl;
}




void SatInfo::exportPosVel(const std::string &filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing.\n";
        return;
    }

    // Write header
    file << "Satellite,Position_X,Position_Y,Position_Z,Velocity_X,Velocity_Y,Velocity_Z\n";
    for (int iSat = 0; iSat < numOfSatellites; iSat++) {
        std::vector posVel = satPosVel[iSat];
        file << selectedSatellites[iSat].satID << ',' << std::fixed << std::setprecision(10) << posVel[0] << ','
        << std::fixed << std::setprecision(10) << posVel[1] << ','
        << std::fixed << std::setprecision(10) << posVel[2] << ','
        << std::fixed << std::setprecision(10) << posVel[3] << ','
        << std::fixed << std::setprecision(10) << posVel[4] << ','
        << std::fixed << std::setprecision(10) << posVel[5] << '\n';
    }
    file.close();
    std::cout << "Satellite positions and velocities exported to " << filename << " successfully.\n";
}

void SatInfo::exportEphemeris(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing.\n";
        return;
    }

    // Write header
    file << "Constellation,SatID,GpsTime,EphemerisData\n";
    for (const auto& eph : selectedSatellites) {
        file << eph.constellation << ',' << eph.satID << ','
             << std::fixed << std::setprecision(2) << eph.gpsTime; // Ensures full precision
        for (const auto& data : eph.ephemerisData) {
            file << ',' <<  std::fixed << std::setprecision(31) << data;
        }
        file << '\n';
    }

    file.close();
    std::cout << "Ephemeris data exported to " << filename << " successfully.\n";
}


