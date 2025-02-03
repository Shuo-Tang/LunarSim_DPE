/* signalGenerator.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for signalGenerator class
 * The class helps to generate the baseband received signal on the moon
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include "signalGenerator.h"
#include "geoFunctions.h"


SignalGenerator::SignalGenerator(const SatInfo &satInfo, const SettingsLoader &navSettings) {
    // load settings
    this->numOfSatellites = satInfo.numOfSatellites;
    this->numOfEpochs = satInfo.numOfSatUpdates;
    this->receiverPosLLA = {navSettings.receiverPosLLA.lat, navSettings.receiverPosLLA.lon, navSettings.receiverPosLLA.alt};
    this->moonPosVelPath = navSettings.moonPosVelPath;
    this->fs = navSettings.samplingFreq;

    // read moon position and velocity (pre-generated)
    readMoonPosVel();
    std::cout << "moon position: " << moonPosVel[0][0] << ", " << moonPosVel[0][1] << ", " << moonPosVel[0][2] << std::endl;
    std::cout << "receiver position: " << receiverPosLLA[0] << ", " << receiverPosLLA[1] << ", " << receiverPosLLA[2] << std::endl;
    std::cout << "UTC time: " << navSettings.utcTimeEpoch.tm_year << "-" << navSettings.utcTimeEpoch.tm_mon << "-" << navSettings.utcTimeEpoch.tm_mday << " "
              << navSettings.utcTimeEpoch.tm_hour << ":" << navSettings.utcTimeEpoch.tm_min << ":" << navSettings.utcTimeEpoch.tm_sec << std::endl;
    // compute receiver position in ECEF
    for (int iEpoch = 0; iEpoch < numOfEpochs; iEpoch++) {
        std::vector<double> posECEF = moonLla2ecef(navSettings.utcTimeEpoch, receiverPosLLA[0], receiverPosLLA[1], receiverPosLLA[2], moonPosVel[iEpoch]);
        this->receiverPosECEF.push_back(posECEF);
    }
    // estimate receiver velocity in ECEF
    std::vector<double> velECEF(3, 0.0);
    for (int iEpoch = 0; iEpoch < numOfEpochs - 1; iEpoch++) {
        for (int i = 0; i < 3; i++) {
            velECEF[i] = (this->receiverPosECEF[iEpoch + 1][i] - this->receiverPosECEF[iEpoch][i]) / navSettings.satUpdateRate;
        }
        this->receiverVelECEF.push_back(velECEF);
    }
    this->receiverVelECEF.push_back(velECEF);
}

void SignalGenerator::readMoonPosVel() {
    std::ifstream file(moonPosVelPath);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + moonPosVelPath);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<double> row;
        std::string value;

        while (std::getline(lineStream, value, ',')) {
            row.push_back(std::stod(value));
        }

        if (row.size() != 6) {
            throw std::runtime_error("Malformed line in file: " + line);
        }
        moonPosVel.push_back(row);
    }
    file.close();
}


void SignalGenerator::generateSignalAll(const SatInfo &satInfo, const SettingsLoader &navSettings) {
    // group the satellites by constellation
    for (int i = 0; i < numOfSatellites; i++) {
        switch (satInfo.selectedSatellites[i].constellation) {
            case 'G':
                gpsIndex.push_back(i);
                break;
            case 'R':
                glonassIndex.push_back(i);
                break;
            case 'E':
                galileoIndex.push_back(i);
                break;
            case 'C':
                beidouIndex.push_back(i);
                break;
            default:
                break;
        }
    }

    if (navSettings.enableGPS) {
        std::cout << "Generating and writing GPS signal..." << std::endl;
        generateSignalGNSS('G', satInfo, navSettings);
        std::cout << std::endl;
    }
    if (navSettings.enableGLONASS) {
        std::cout << "Generating and writing GLONASS signal..." << std::endl;
        generateSignalGNSS('R', satInfo, navSettings);
        std::cout << std::endl;
    }
    if (navSettings.enableGalileo) {
        std::cout << "Generating and writing Galileo signal..." << std::endl;
        generateSignalGNSS('E', satInfo, navSettings);
        std::cout << std::endl;
    }

}

void SignalGenerator::generateSignalGps(const std::vector<int>& gpsIndex, const SatInfo &satInfo, const SettingsLoader &navSettings) {
    int samplePerMs = fs / 1000; // number of samples in 1ms
    // Loop for every update. Each update is valid for satUpdateRate seconds.
    for (int iUpdate = 0; iUpdate < numOfEpochs; iUpdate++) {
        // loop for every satellite
        std::vector<std::complex<double>> signalTx1ms(samplePerMs, std::complex<double>(0.0, 0.0));
        for (int iSat : gpsIndex) {
            std::string satID = satInfo.selectedSatellites[iSat].satID;
            const int prn = std::stoi(satID.substr(1));
            int satIndex = numOfSatellites * iUpdate + iSat;
            // Generate signal for every 1ms. Concatenate the signal for the whole update (satUpdateRate * 1e3 ms)
            Eigen::Vector3d satPos(satInfo.satPosVel[satIndex][0], satInfo.satPosVel[satIndex][1], satInfo.satPosVel[satIndex][2]);
            Eigen::Vector3d satVel(satInfo.satPosVel[satIndex][3], satInfo.satPosVel[satIndex][4], satInfo.satPosVel[satIndex][5]);
            Eigen::Vector3d userPos(receiverPosECEF[iUpdate][0], receiverPosECEF[iUpdate][1], receiverPosECEF[iUpdate][2]);
            Eigen::Vector3d userVel(receiverVelECEF[iUpdate][0], receiverVelECEF[iUpdate][1], receiverVelECEF[iUpdate][2]);
            // Filter the valid satellites within the lobe range
            const int lobeFlag = lobeRange(satPos, userPos, navSettings);
            std::cout << "Satellite: " << satID << " Lobe: " << lobeFlag << std::endl;
            if (lobeFlag == 0) {
                continue;
            }
            // prepare the signal parameters
            std::vector<int> caCode = genCAcode(prn);
            double range = (satPos - userPos).norm(); // geometric range
            double delay = range / c; // delay
            Eigen::Vector3d unitVec = (satPos - userPos) / range; // los unit vector
            double doppler = - (satVel - userVel).dot(unitVec) * fGpsL1 / c; // doppler frequency
            double fracDelay = fmod(delay, codePeriodGps); // fractional delay
            double prevNCOIndex = -fracDelay / TcGps; // initial NCO index
            // Generate local replica signal
            for (int i = 0; i < samplePerMs; i++) {
                int index = positiveMod(static_cast<int>(std::round(prevNCOIndex + static_cast<double>(i) / fs / TcGps)), caCodeLength);
                double phase = 2.0 * M_PI * doppler * (static_cast<double>(i) / fs);
                // summation of the signal among all satellites
                signalTx1ms[i] += static_cast<double>(caCode[index % caCodeLength]) * std::polar(1.0, phase);
            }
        }
        // copy the signal for the whole update
        int numOfMs = navSettings.satUpdateRate * 1000;
        std::vector<std::complex<double>> signalTx = copySignal(signalTx1ms, numOfMs);
        // Append signalTx to gpsSignal
        gpsSignal.insert(gpsSignal.end(), signalTx.begin(), signalTx.end());
    }

}

void SignalGenerator::generateSignalGNSS(char constellation, const SatInfo &satInfo, const SettingsLoader &navSettings) {
    int numberOfSamples = 0;
    double fc = 0.0;
    double codePeriod = 0.0;
    double Tc = 0.0;
    int codeLength = 0;
    std::vector<int> groupIndex;
    if (constellation == 'G') {
        groupIndex = gpsIndex;
        numberOfSamples = fs * codePeriodGps;
        fc = fGpsL1;
        codePeriod = codePeriodGps;
        Tc = TcGps;
        codeLength = caCodeLength;
    } else if (constellation == 'E') {
        groupIndex = galileoIndex;
        numberOfSamples = fs * codePeriodGal;
        fc = fGalE1;
        codePeriod = codePeriodGal;
        Tc = TcGal;
        codeLength = galE1CodeLength;
    } else if (constellation == 'R') {
        groupIndex = glonassIndex;
        numberOfSamples = fs * codePeriodGlo;
        codePeriod = codePeriodGlo;
        Tc = TcGlo;
        gloSignal.resize(groupIndex.size());
        codeLength = gloL1CodeLength;
    }
    else {
        numberOfSamples = fs * codePeriodGps;
        fc = fGpsL1;
        codePeriod = codePeriodGps;
        Tc = TcGps;
        codeLength = caCodeLength;
    }
    int numOfPeriods = navSettings.satUpdateRate / codePeriod;
    maxUpdates = maxElements / (numOfPeriods * numberOfSamples);
    if (constellation == 'R') {
        maxUpdates /= groupIndex.size();
    }
    int nUpdates = 1;
    // Loop for every update. Each update is valid for satUpdateRate seconds.
    for (int iUpdate = 0; iUpdate < numOfEpochs; iUpdate++) {
        // Display the progress bar
        displayProgressBar(iUpdate + 1, numOfEpochs);
        // loop for every satellite
        // prepare the transmitted signal for every satellite
        std::vector<std::vector<std::complex<double>>> signalTxPeriod;
        if (constellation != 'R') {
            signalTxPeriod = std::vector(1, std::vector(numberOfSamples, std::complex<double>(0.0, 0.0)));
        }
        else {
            signalTxPeriod = std::vector(24, std::vector(numberOfSamples, std::complex<double>(0.0, 0.0)));
        }
        std::vector<int> validGpsPRN;
        std::vector<int> validGloPRN;
        std::vector<int> validGalPRN;

        for (int iSat : groupIndex) {
            std::string satID = satInfo.selectedSatellites[iSat].satID;
            const int prn = std::stoi(satID.substr(1));
            int satIndex = numOfSatellites * iUpdate + iSat;
            // Generate signal for every 1ms. Concatenate the signal for the whole update (satUpdateRate * 1e3 ms)
            Eigen::Vector3d satPos(satInfo.satPosVel[satIndex][0], satInfo.satPosVel[satIndex][1], satInfo.satPosVel[satIndex][2]);
            Eigen::Vector3d satVel(satInfo.satPosVel[satIndex][3], satInfo.satPosVel[satIndex][4], satInfo.satPosVel[satIndex][5]);
            Eigen::Vector3d userPos(receiverPosECEF[iUpdate][0], receiverPosECEF[iUpdate][1], receiverPosECEF[iUpdate][2]);
            Eigen::Vector3d userVel(receiverVelECEF[iUpdate][0], receiverVelECEF[iUpdate][1], receiverVelECEF[iUpdate][2]);
            // Filter the valid satellites within the lobe range
            const int lobeFlag = lobeRange(satPos, userPos, navSettings);
            // std::cout << "Satellite: " << satID << " Lobe: " << lobeFlag << std::endl;
            // std::cout << "sat X: " << satPos[0] << " Y: " << satPos[1] << " Z: " << satPos[2] << std::endl;
            // std::cout << "user X: " << userPos[0] << " Y: " << userPos[1] << " Z: " << userPos[2] << std::endl;
            if (lobeFlag == 0) {
                continue;
            }
            // record the valid satellites

            if (constellation == 'G') {
                validGpsPRN.push_back(prn);
            }
            else if (constellation == 'E') {
                validGalPRN.push_back(prn);
            }
            else if (constellation == 'R') {
                validGloPRN.push_back(prn);
            }
            // prepare the signal parameters
            std::vector<int> code;
            if (constellation == 'G') {
                code = genCAcode(prn);
            }
            else if (constellation == 'E') {
                code = genGalE1(prn);
            }
            else if (constellation == 'R') {
                code = genGloL1(prn);
                fc = fGloL1[(prn - 1) % 14];
            }
            double range = (satPos - userPos).norm(); // geometric range
            double delay = range / c; // delay
            Eigen::Vector3d unitVec = (satPos - userPos) / range;
            // los unit vector
            double doppler = - (satVel - userVel).dot(unitVec) * fc / c; // doppler frequency
            // std::cout << "Doppler: " << doppler << std::endl;
            double fracDelay = fmod(delay, codePeriod); // fractional delay
            double prevNCOIndex = -fracDelay / Tc; // initial NCO index
            // Generate local replica signal
            for (int i = 0; i < numberOfSamples; i++) {
                int index = positiveMod(static_cast<int>(std::round(prevNCOIndex + static_cast<double>(i) / fs / Tc)), codeLength);
                double phase = 2.0 * M_PI * doppler * (static_cast<double>(i) / fs);
                // summation of the signal among all satellites
                if (constellation != 'R') {
                    signalTxPeriod[0][i] += static_cast<double>(code[index % codeLength]) * std::polar(1.0, phase);
                }
                else {
                    signalTxPeriod[prn - 1][i] += static_cast<double>(code[index % codeLength]) * std::polar(1.0, phase);
                }
            }
        }
        // copy the signal for the all updates
        if (constellation != 'R') {
            std::vector<std::complex<double>> signalTx = copySignal(signalTxPeriod[0], numOfPeriods);
            if (constellation == 'G') {
                gpsSignal.insert(gpsSignal.end(), signalTx.begin(), signalTx.end());
                // Check if gpsSignal exceeds the maxElements limit
                if (iUpdate == nUpdates * maxUpdates) {
                    saveSignal(constellation, navSettings, satInfo); // Save the current data
                    gpsSignal.clear(); // Clear the vector for the next batch
                    nUpdates++;
                }
            }
            else if (constellation == 'E') {
                galSignal.insert(galSignal.end(), signalTx.begin(), signalTx.end());
                // Check if galSignal exceeds the maxElements limit
                if (iUpdate == nUpdates * maxUpdates) {
                    saveSignal(constellation, navSettings, satInfo); // Save the current data
                    galSignal.clear(); // Clear the vector for the next batch
                    nUpdates++;
                }
            }
        }
        else {
            for (int iSat = 0; iSat < groupIndex.size(); iSat++) {
                std::vector<std::complex<double>> signalTx = copySignal(signalTxPeriod[iSat], numOfPeriods);
                gloSignal[iSat].insert(gloSignal[iSat].end(), signalTx.begin(), signalTx.end());
            }
            // Check if gloSignal exceeds the maxElements limit
            if (iUpdate == nUpdates * maxUpdates) {
                saveSignal(constellation, navSettings, satInfo); // Save the current data
                gloSignal.clear(); // Clear the vector for the next batch
                nUpdates++;
            }
        }
        // save the valid satellites
        if (constellation == 'G') {
            validGpsSatellites.push_back(validGpsPRN);
        }
        else if (constellation == 'E') {
            validGalSatellites.push_back(validGalPRN);
        }
        else if (constellation == 'R') {
            validGloSatellites.push_back(validGloPRN);
        }
    }
    std::cout << std::endl;
    // Save any remaining data in signal
    std::cout << "Finishing signal data writing..." << std::endl;
    if (constellation == 'G' && !gpsSignal.empty()) {
        saveSignal(constellation, navSettings, satInfo);
    }
    if (constellation == 'E' && !galSignal.empty()) {
        saveSignal(constellation, navSettings, satInfo);
    }
    if (constellation == 'R' && !gloSignal.empty()) {
        saveSignal(constellation, navSettings, satInfo);
    }
    std::cout << "Signal data writing completed." << std::endl;
}


int SignalGenerator::lobeRange(const Eigen::Vector3d &satPos, const Eigen::Vector3d &receiverPos, const SettingsLoader &navSettings) {
    // compute distance from Earth's center to the LOS
    Eigen::Vector3d los = receiverPos - satPos;
    Eigen::Vector3d crossProd = satPos.cross(los);
    double d =  crossProd.norm() / los.norm();
    // if the satellite is covered by the Earth or facing to the Earth
    if (d <= earthRadius) {
        return 0;
    }

    // Compute the cosine of the angle between vectors a and b
    double cosTheta = los.dot(-satPos) / (los.norm() * (-satPos).norm());
    double theta = std::acos(cosTheta) * (180.0 / M_PI);

    if (theta <= navSettings.mainLobeAngle) {
        return 1;  // Inside main lobe
    } else if (theta <= navSettings.sideLobeAngle) {
        return 2;  // Between main and side lobe
    } else {
        return 0;  // Outside side lobe
    }

}


void SignalGenerator::saveSignal(const char constellation, const SettingsLoader &navSettings, const SatInfo &satInfo) {
    if (constellation == 'G') {
        // open the file in append mode
        std::ios_base::openmode mode = std::ios::binary;
        mode |= std::ios::app; // Open in append mode

        std::ofstream file(navSettings.gpsSignalPath + ".dat", mode);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << navSettings.gpsSignalPath << ".dat" << " for writing.\n";
            return;
        }

        // Write the complex signal to the binary file
        for (const auto& value : gpsSignal) {
            double real = value.real();
            double imag = value.imag();
            file.write(reinterpret_cast<const char*>(&real), sizeof(double));
            file.write(reinterpret_cast<const char*>(&imag), sizeof(double));
        }

        file.close();
    }

    if (constellation == 'R') {
        // Iterate through each PRN and corresponding signal
        for (int iSat = 0; iSat < gloSignal.size(); ++iSat) {
            const auto& signal = gloSignal[iSat];
            std::string satID = satInfo.selectedSatellites[glonassIndex[iSat]].satID;
            const int prn = std::stoi(satID.substr(1));

            // Create a filename with zero-padded PRN (e.g., PRN01, PRN02, ...)
            std::ostringstream filenameStream;
            filenameStream << navSettings.gloSignalPath << std::setw(2) << std::setfill('0') << prn << ".dat";
            std::string filename = filenameStream.str();

            // open the file in append mode
            std::ios_base::openmode mode = std::ios::binary;
            mode |= std::ios::app; // Open in append mode

            std::ofstream file(filename, mode);
            if (!file.is_open()) {
                std::cerr << "Error: Could not open file " << filename << " for writing.\n";
                return;
            }

            // Write the complex signal to the binary file
            for (const auto& value : signal) {
                double real = value.real();
                double imag = value.imag();
                file.write(reinterpret_cast<const char*>(&real), sizeof(double));
                file.write(reinterpret_cast<const char*>(&imag), sizeof(double));
            }

            file.close();
        }

    }

    if (constellation == 'E') {
        // open the file in append mode
        std::ios_base::openmode mode = std::ios::binary;
        mode |= std::ios::app; // Open in append mode

        std::ofstream file(navSettings.galSignalPath + ".dat", mode);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << navSettings.galSignalPath << ".dat" << " for writing.\n";
            return;
        }

        // Write the complex signal to the binary file
        for (const auto& value : galSignal) {
            double real = value.real();
            double imag = value.imag();
            file.write(reinterpret_cast<const char*>(&real), sizeof(double));
            file.write(reinterpret_cast<const char*>(&imag), sizeof(double));
        }

        file.close();
    }
}







