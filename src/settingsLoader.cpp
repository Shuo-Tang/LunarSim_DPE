/* settingsLoader.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for settingsLoader class
 * The class helps to load the settings from a YAML file
 */

#include "settingsLoader.h"
#include "timeFunctions.h"

SettingsLoader::SettingsLoader() {
    enableGPS = true;
    enableGLONASS = false;
    enableGalileo = false;
    enableBeiDou = false;
    gpsFilePath = "";
    glonassFilePath = "";
    galileoFilePath = "";
    beidouFilePath = "";
    timeEpoch = "";
    timeDuration = 0;
    receiverPosLLA = LLA{0, 0, 0};
    enableLunarSatellites = false;
    lunarSatPosLLA = std::vector<LLA>();
}

SettingsLoader::~SettingsLoader() = default;

void SettingsLoader::loadSettings(const std::string& fileName) {
    try {
        YAML::Node config = YAML::LoadFile(fileName);

        satPosVelPath = config["filePaths"]["satPosVelOutput"].as<std::string>();
        satEphPath = config["filePaths"]["satEphOutput"].as<std::string>();
        moonPosVelPath = config["filePaths"]["moonPosVel"].as<std::string>();

        enableGPS = config["enableGPS"].as<bool>();
        enableGLONASS = config["enableGLONASS"].as<bool>();
        enableGalileo = config["enableGalileo"].as<bool>();
        enableBeiDou = config["enableBeiDou"].as<bool>();

        if (enableGPS) {
            gpsFilePath = config["filePaths"]["GPS"].as<std::string>();
            gpsSignalPath = config["filePaths"]["gpsSignalOutput"].as<std::string>();
            enabledConstellations.push_back('G');
        }
        if (enableGLONASS) {
            glonassFilePath = config["filePaths"]["GLONASS"].as<std::string>();
            gloSignalPath = config["filePaths"]["glonassSignalOutput"].as<std::string>();
            enabledConstellations.push_back('R');
        }
        if (enableGalileo) {
            galileoFilePath = config["filePaths"]["Galileo"].as<std::string>();
            galSignalPath = config["filePaths"]["galileoSignalOutput"].as<std::string>();
            enabledConstellations.push_back('E');
        }
        if (enableBeiDou) {
            beidouFilePath = config["filePaths"]["BeiDou"].as<std::string>();
            enabledConstellations.push_back('C');
        }

        timeEpoch = config["timeEpoch"].as<std::string>();
        utcTimeEpoch = str2ctime(timeEpoch.substr(2));
        timeDuration = config["timeDuration"].as<int>();
        satUpdateFrequency = config["satFrequency"].as<double>();
        satUpdateRate = 1.0 / satUpdateFrequency;


        receiverPosLLA.lat = config["lunarReceiver"]["latitude"].as<double>();
        receiverPosLLA.lon = config["lunarReceiver"]["longitude"].as<double>();
        receiverPosLLA.alt = config["lunarReceiver"]["altitude"].as<double>();

        enableLunarSatellites = config["enableLunarSatellites"].as<bool>();
        if (enableLunarSatellites) {
            for (const auto& satPos : config["LunarSatellite"]) {
                LLA satLLA{};
                satLLA.lat = satPos["latitude"].as<double>();
                satLLA.lon = satPos["longitude"].as<double>();
                satLLA.alt = satPos["altitude"].as<double>();
                lunarSatPosLLA.push_back(satLLA);
            }
        }

        mainLobeAngle = config["lobeAngleLimit"]["mainLobe"].as<double>();
        sideLobeAngle = config["lobeAngleLimit"]["sideLobe"].as<double>();

        enableSignalOutput = config["enableSignalOutput"].as<bool>();
        samplingFreq = config["samplingFrequency"].as<double>();

        enable2SP = config["enable2SP"].as<bool>();
        enableDPE = config["enableDPE"].as<bool>();
        coherentIntegrationTime = config["coherentIntegrationTime"].as<int>();
    } catch (const YAML::Exception& e) {
        std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
    }
}