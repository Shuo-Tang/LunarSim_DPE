/* twoStepSolver.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for signalGenerator class
 * The class helps to solve lunar position with traditional two-step method
 */

#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <complex>
#include <fftw3.h>
#include "satInfo.h"
#include "settingsLoader.h"
#include "signalGenerator.h"
#include "signalFunctions.h"

class TwoStepSolver {
  public:
    TwoStepSolver();

    void solvePVT(SatInfo& satInfo, SignalGenerator& signalGenerator, const SettingsLoader& navSettings);

    std::vector<int> getLocalReplica2SP(char constellation, int prn, double fs, int coherentIntegrationTime);

    std::pair<double, double> getFracPseudorange(const std::vector<std::complex<double>>& signalRx, const std::vector<int>& localCode, double fs);


    std::vector<std::vector<double>> pseudorange;

};