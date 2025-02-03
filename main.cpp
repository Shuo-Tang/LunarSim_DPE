#include <iostream>

#include "settingsLoader.h"
#include "navInfo.h"
#include "satInfo.h"
#include "timeFunctions.h"
#include "signalGenerator.h"
#include "twoStepSolver.h"



int main() {
    std::cout << "Hello, World!" << std::endl;
    /*=====================================
     === set input and output file path ===
     ======================================*/
    // load path: settings
    const std::string settingsPath = "../importData/navSettings.yaml";

    // load settings
    SettingsLoader settingsLoader;
    settingsLoader.loadSettings(settingsPath);

    // load settings to navi file information
    NavInfo navInfo;
    // read navi information from Rinex files
    navInfo.readNavInfo(settingsLoader);

    // load satellite ephemeris
    std::cout << "Start epoch: " << std::endl;
    std::cout << settingsLoader.timeEpoch << std::endl;
    SatInfo satInfo;
    // select satellites
    satInfo.loadSatEph(navInfo, settingsLoader);
    std::cout << "Selected satellites:" << std::endl;
    std::cout << "Total " << satInfo.numOfSatellites << " satellites" << std::endl;
    for (int i = 0; i < satInfo.numOfSatellites; i++) {
        std::cout << gtime2str(satInfo.selectedSatellites[i].gtime) << " " << satInfo.selectedSatellites[i].satID << std::endl;
    }

    // generate satellite position and velocity
    std::cout << "Generating satellites position and velocity for " << satInfo.timeDuration << " s" << std::endl;
    std::cout << "Update every " << satInfo.satUpdateRate << " s" << std::endl;
    satInfo.generateSatPosVel();
    satInfo.exportEphemeris(settingsLoader.satEphPath);
    satInfo.exportPosVel(settingsLoader.satPosVelPath);
    std::cout << "satellites Information is saved." << std::endl;

    // generate signal
    SignalGenerator signalGenerator(satInfo, settingsLoader);
    signalGenerator.generateSignalAll(satInfo, settingsLoader);
    std::cout << "Signal is generated and saved." << std::endl;
    // save signal


    // solve PVT
    // if (settingsLoader.enable2SP) {
    //     TwoStepSolver twoStepSolver;
    //     twoStepSolver.solvePVT(satInfo, signalGenerator, settingsLoader);
    //     for (int i = 0; i < settingsLoader.coherentIntegrationTime; i++) {
    //         std::cout << "Pseudorange at " << i << " ms:" << std::endl;
    //         for (int j = 0; j < satInfo.numOfSatellites; j++) {
    //             std::cout << satInfo.selectedSatellites[j].satID << ": " << twoStepSolver.pseudorange[i][j] << std::endl;
    //         }
    //     }
    // }
    // display each pseudorange


    while (true) {
        std::cout << "Press any key to exit." << std::endl;
        // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (std::cin.get()) {
            break;
        }
    }
    return 0;
}
