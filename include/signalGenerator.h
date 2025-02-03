/* signalGenerator.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for signalGenerator class
 * The class helps to generate the baseband received signal on the moon
 */

#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include <deque>
#include "satInfo.h"
#include "settingsLoader.h"
#include "signalFunctions.h"

constexpr double c = 299792458.0;
constexpr double fGpsL1 = 1575.42e6;
constexpr double fGalE1 = 1575.42e6;
const std::vector<double> fGloL1 = {
    1598.0625e6, // k = -7
    1598.6250e6, // k = -6
    1599.1875e6, // k = -5
    1599.7500e6, // k = -4
    1600.3125e6, // k = -3
    1600.8750e6, // k = -2
    1601.4375e6, // k = -1
    1602.0000e6, // k = 0
    1602.5625e6, // k = +1
    1603.1250e6, // k = +2
    1603.6875e6, // k = +3
    1604.2500e6, // k = +4
    1604.8125e6, // k = +5
    1605.3750e6  // k = +6
};
constexpr double codePeriodGps = 1.0e-3;
constexpr double codePeriodGal = 4.0e-3;
constexpr double codePeriodGlo = 1.0e-3;
constexpr double TcGps = 1.0 / 1.023e6;
constexpr double TcGal = 4.0 / 4.092e6;
constexpr double TcGlo = 1.0 / 511.0e6;
constexpr int caCodeLength = 1023;
constexpr int galE1CodeLength = 4092;
constexpr int gloL1CodeLength = 511;
constexpr int maxElements = 2046000000;

class SignalGenerator {
  public:
    SignalGenerator(const SatInfo& satInfo, const SettingsLoader& navSettings);

    void generateSignalAll(const SatInfo &satInfo, const SettingsLoader &navSettings);

    void generateSignalGps(const std::vector<int>& index, const SatInfo &satInfo, const SettingsLoader &navSettings);

    void generateSignalGNSS(char constellation, const SatInfo &satInfo, const SettingsLoader &navSettings);

    int lobeRange(const Eigen::Vector3d &satPos, const Eigen::Vector3d &receiverPos, const SettingsLoader &navSettings);

    void readMoonPosVel();

    void saveSignal(const char, const SettingsLoader &navSettings, const SatInfo &satInfo);


    int numOfSatellites;
    int numOfEpochs;
    int maxUpdates;//  this is the maximum number of updates which can be restored in one signal vector


    std::string moonPosVelPath;
    std::vector<double> receiverPosLLA;
    std::vector<std::vector<double>> receiverPosECEF;
    std::vector<std::vector<double>> receiverVelECEF;
    std::vector<std::vector<double>> moonPosVel;

    double fs;

    std::vector<std::vector<int>> validGpsSatellites;
    std::vector<std::vector<int>> validGalSatellites;
    std::vector<std::vector<int>> validGloSatellites;

    std::vector<int> gpsIndex;
    std::vector<int> glonassIndex;
    std::vector<int> galileoIndex;
    std::vector<int> beidouIndex;

    std::vector<std::complex<double>> gpsSignal;
    std::vector<std::complex<double>> galSignal;
    std::vector<std::vector<std::complex<double>>> gloSignal;


};


