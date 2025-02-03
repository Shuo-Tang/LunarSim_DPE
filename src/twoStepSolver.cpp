/* twoStepSolver.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for signalGenerator class
 * The class helps to solve lunar position with traditional two-step method
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include "twoStepSolver.h"

void saveCAF2SP(const std::string& filename, const std::vector<std::vector<double>>& caf) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    for (const auto& row : caf) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ","; // Add comma between elements
            }
        }
        file << "\n"; // Newline after each row
    }

    file.close();
    std::cout << "Data saved to " << filename << "\n";
}


TwoStepSolver::TwoStepSolver()= default;

void TwoStepSolver::solvePVT(SatInfo &satInfo, SignalGenerator &signalGenerator, const SettingsLoader &navSettings) {
    double fs = navSettings.samplingFreq;
    int numOfMs = navSettings.timeDuration * 1000;
    int samplePerMs = fs / 1000;
    std::cout << "Generating pseudorange..." << std::endl;
    if (navSettings.enableGPS) {
        // get pseudorange and doppler frequency for each ms (with coherent integration time)
        for (int i = 0; i <= numOfMs - navSettings.coherentIntegrationTime + 1; i++) {
            // display progress
            displayProgressBar(i + 1, numOfMs - navSettings.coherentIntegrationTime + 1, 50);
            // extract the received signal in coherent integration time
            int left = i * samplePerMs;
            int right = (i + navSettings.coherentIntegrationTime) * samplePerMs - 1;
            std::vector<std::complex<double>> receivedSignal(navSettings.coherentIntegrationTime * samplePerMs, std::complex<double>(0.0, 0.0));
            for (int iSample = 0; iSample < navSettings.coherentIntegrationTime * samplePerMs; iSample++) {
                receivedSignal[iSample] = signalGenerator.gpsSignal[left + iSample];
            }
            // get update epoch
            int iUpdate = i / (navSettings.satUpdateRate * 1000);
            std::vector<int> validPRN = signalGenerator.validGpsSatellites[iUpdate];
            // loop for every satellite
            std::vector<double> prAtMs = std::vector(validPRN.size(), 0.0); // pseudorange at this ms
            for (int iSat = 0; iSat < validPRN.size(); iSat++) {
                int prn = validPRN[iSat];
                // generate the local replica signal
                std::vector localCode = getLocalReplica2SP('G', prn, fs, navSettings.coherentIntegrationTime);
                // get pseudorange measurement
                auto [pr, _] = getFracPseudorange(receivedSignal, localCode, fs);
                prAtMs[iSat] = pr;
            }
            pseudorange.push_back(prAtMs);
        }
    }
}

std::vector<int> TwoStepSolver::getLocalReplica2SP(char constellation, int prn, double fs, int coherentIntegrationTime) {
    int numberOfSamples = fs * coherentIntegrationTime / 1000;
    std::vector<int> localReplica = std::vector(numberOfSamples, 0);
    // prepare the signal parameters
    double Tc = 0.0;
    int codeLength = 0;
    std::vector<int> code;
    if (constellation == 'G') {
        code = genCAcode(prn);
        Tc = TcGps;
        codeLength = caCodeLength;
    }
    else if (constellation == 'E') {
        code = genGalE1(prn);
        Tc = TcGal;
        codeLength = galE1CodeLength;
    }
    else if (constellation == 'R') {
        code = genGloL1(prn);
        Tc = TcGlo;
        codeLength = gloL1CodeLength;
    }
    // apply sampling
    for (int i = 0; i < numberOfSamples; i++) {
        int index = positiveMod(static_cast<int>(std::round(static_cast<double>(i) / fs / Tc)), codeLength);
        // summation of the signal among all satellites
        localReplica[i] = code[index % codeLength];
    }
    return localReplica;
}


std::pair<double, double> TwoStepSolver::getFracPseudorange(const std::vector<std::complex<double>>& signalRx, const std::vector<int>& localCode, double fs) {
    std::vector<int> freqBin = getFreqBin(-20000, 20000, 100);
    int numberOfSample = signalRx.size();
    double maxCaf = -1.0;
    int maxFreqIndex = -1;
    int maxSampleIndex = -1;
    // std::vector<std::vector<double>> caf = std::vector(freqBin.size(), std::vector(numberOfSample, 0.0));
    for (int iFreq = 0; iFreq < freqBin.size(); iFreq++) {
        // get the frequency bin
        int fd = freqBin[iFreq];
        // compute phase
        std::vector<double> phase = std::vector(numberOfSample, 0.0);
        for (int iSample = 0; iSample < numberOfSample; iSample++) {
            phase[iSample] = 2.0 * M_PI * fd * (static_cast<double>(iSample) / fs);
        }
        // compute the signal with phase
        std::vector<std::complex<double>> localSignal = std::vector(numberOfSample, std::complex<double>(0.0, 0.0));
        for (int iSample = 0; iSample < numberOfSample; iSample++) {
            localSignal[iSample] =  static_cast<double>(localCode[iSample]) * std::polar(1.0, phase[iSample]);
        }
        // compute the ifft(fft(x) .* conj(fft(y)))
        std::vector<std::complex<double>> inverseFFT = shiftCorrelation(signalRx, localSignal);
        // compute the cross-ambiguity function
        for (int iSample = 0; iSample < numberOfSample; iSample++) {
            // caf[iFreq][iSample] = std::abs(inverseFFT[iSample]);
            double caf = std::abs(inverseFFT[iSample]);
            if (caf > maxCaf) {
                maxCaf = caf;
                maxFreqIndex = iFreq;
                maxSampleIndex = iSample;
            }
        }
    }
    // save the cross-ambiguity function
    // saveCAF2SP("../exportData/caf2sp.csv", caf);
    // prepare the output: pseudorange and doppler frequency
    double fd = freqBin[maxFreqIndex];
    double pr = static_cast<double>(maxSampleIndex) / fs * c;
    return std::make_pair(pr, fd);
}





