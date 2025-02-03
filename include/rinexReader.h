/* rinexReader.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header filer for rinexReader.cpp
 * This file contains the necessary includes and declarations for the rinexReader.cpp file
 * The class helps to read the navigation information from a RINEX file
 */

# pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <chrono>
#include <tuple>
#include <algorithm>
#include "ephemeris.h"
#include "timeFunctions.h"
#include <regex>

const std::unordered_map<char, int> numberOfExtraLines = {
    {'G', 7},
    {'R', 3},
    {'E', 7}
};

const std::unordered_map<char, int> numberOfEphemerisData = {
    {'G', 29},
    {'R', 15},
    {'E', 28} // for Rinex 3.xx
};

class RinexReader {
public:
    // constructor
    explicit RinexReader(const std::string& filePath);
    // destructor
    ~RinexReader();

    // parse navigation file
    void EphemerisParser(std::vector<Ephemeris> &satEphemeris);

    // parse the header of RINEX file
    void parseRinex2Header(std::ifstream& file);
    void parseRinex3Header(std::ifstream& file);

    // get version and type of RINEX file
    void getRinexInfo(std::ifstream &file);

    // get version of RINEX file
    static std::tuple<std::string, bool> getRinexVersion(const std::string& line);

    // read navigation file
    void readRinexNav(std::ifstream& file, std::vector<Ephemeris> &satEphemeris);

    // all constellations and satellites
    // G: GPS, R: GLONASS, E: Galileo, C: BeiDou
    std::vector<Ephemeris> satEphemeris;
    // file path
    std::string rinexPath;
    // header information
    std::unordered_map<std::string, std::string> header;
    // rinex info
    struct RinexInfo {
        std::string version; // version {2.XX, 3.XX}
        char fileType; // file type {O, C, N, G, E}
        std::string rinexType; // rinex type {nav, obs}
        char constellation; // constellation {G, R, E, C}
        // ionospheric corrections: ION ALPHA and ION BETA / GPSA and GPSB / GAL...
        std::unordered_map <std::string, std::vector<float>> ionoCorrections;
        // system time corrections: GPUT / GAUT... (not included for now)
        std::unordered_map <std::string, std::vector<float>> timeCorrections;
    }rinexInfo;
};

