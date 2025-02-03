/* rinexReader.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for rinexReader class
 * The class helps to read the navigation information from a RINEX file
 * From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 */

#include "rinexReader.h"
#include "timeFunctions.h"

// constructor
RinexReader::RinexReader(const std::string& rinexPath) {
    // initialize the file path
    this->rinexPath = rinexPath;
    // initialize the rinex information
    this->rinexInfo.version = "";
    this->rinexInfo.constellation = '\0';
    this->rinexInfo.fileType = '\0';
    this->rinexInfo.rinexType = "";

}

// destructor
RinexReader::~RinexReader() {

}

// remove the white spaces at the end of a string
std::string rtrim(const std::string& s) {
    const size_t end = s.find_last_not_of(" \t\n\r\f\v");
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// parse navigation file
void RinexReader::EphemerisParser(std::vector<Ephemeris> &satEphemeris) {
    // open the Rinex file
    std::ifstream file(rinexPath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open RINEX file: " + rinexPath);
    }
    // get Rinex information
    getRinexInfo(file);
    // read the Rinex file according to version
    if (rinexInfo.rinexType == "nav") {
        // parse and read the navigation file
        readRinexNav(file, satEphemeris);
    } else {
        throw std::runtime_error("The input file is not a navigation file.");
    }

}

// get information of RINEX file
void RinexReader::getRinexInfo(std::ifstream &file) {
    file.seekg(0, std::ios::beg);
    std::string line;
    getline(file, line);
    if (line[0] == '#' && (line[1] == 'a' || line[1] == 'c' || line[1] == 'd')) {
        rinexInfo.version = line.substr(1, 1);
        rinexInfo.rinexType = "sp3";
    }
    // get rinex information
    std::tuple<std::string, bool> versionTuple = getRinexVersion(line);
    rinexInfo.version = std::get<0>(versionTuple); // rinex version
    if (line[20] == 'O' || line[20] == 'C') { // rinex type
        rinexInfo.rinexType = "obs";
    } else if (line[20] == 'N' || line.substr(20, 20).find("NAV") != std::string::npos) {
        rinexInfo.rinexType = "nav";
    } else {
        rinexInfo.rinexType = line[20];
    }
    rinexInfo.fileType = line[20]; // file type
    if (rinexInfo.version[0] == '2') { // rinex constellation
        switch (rinexInfo.fileType) {
            case 'N':
                rinexInfo.constellation = 'G';
            break;
            case 'G':
                rinexInfo.constellation = 'R';
            break;
            case 'E':
                rinexInfo.constellation = 'E';
            break;
            default:
                rinexInfo.constellation = line[40];
            break;
        }
    } else {
        rinexInfo.constellation = line[40];
    }
}

// get version of RINEX file
std::tuple<std::string, bool> RinexReader::getRinexVersion(const std::string& s) {
    if (s.size() < 2) {
        throw std::invalid_argument("Cannot decode RINEX/SP3 version from line: " + s);
    }

    if (s[0] == '#') {
        if (s.size() < 2 || (s[1] != 'a' && s[1] != 'c' && s[1] != 'd')) {
            throw std::invalid_argument("SP3 versions currently handled: a, c, d");
        }
        return std::make_tuple("0.00", false);
    }

    if (s.size() >= 80 && s.substr(60, 20) != "RINEX VERSION / TYPE" && s.substr(60, 20) != "CRINEX VERS   / TYPE") {
        throw std::invalid_argument("The first line of the RINEX file header is corrupted.");
    }

    try {
        std::string version = s.substr(5, 4);
        bool isCrinex = s.substr(20, 20) == "COMPACT RINEX FORMAT";
        return std::make_tuple(version, isCrinex);
    } catch (const std::exception& e) {
        throw std::invalid_argument("Could not determine file version: " + std::string(e.what()));
    }
}


// read navigation file
void RinexReader::readRinexNav(std::ifstream& file, std::vector<Ephemeris> &satEphemeris) {
    // read Rinex data according to version
    if (rinexInfo.version[0] == '2') {
        // read the header of the RINEX file
        parseRinex2Header(file);
        // save the ionospheric corrections
        if (header.contains("ION ALPHA") and header.contains("ION BETA")) {
            std::vector<float> alpha;
            std::vector<float> beta;
            for (int i = 0; i < 4; i++) {
                std::string alphaStr = header["ION ALPHA"].substr(2 + i * 12, 12);
                std::string betaStr = header["ION BETA"].substr(2 + i * 12, 12);
                std::ranges::replace(alphaStr, 'D', 'E');
                std::ranges::replace(betaStr, 'D', 'E');
                alpha.push_back(std::stof(alphaStr));
                beta.push_back(std::stof(betaStr));
            }
            rinexInfo.ionoCorrections["alpha"] = alpha;
            rinexInfo.ionoCorrections["beta"] = beta;
        }
        // read the ephemeris data
        std::string line;
        int numberOfValues = numberOfEphemerisData.at(rinexInfo.constellation);
        if (rinexInfo.constellation == 'E') {
            numberOfValues--;
        }
        // read the raw data for one satellite
        while (std::getline(file, line)) {
            // Read raw data
            std::string raw = line.substr(22, 57); // From column 22 to 79
            for (int i = 0; i < numberOfExtraLines.at(rinexInfo.constellation); ++i) {
                std::string nextLine;
                if (std::getline(file, nextLine)) {
                    raw += nextLine.substr(3, 77);
                }
            }
            std::ranges::replace(raw, 'D', 'E');
            // Parse raw data
            std::string timeStr = line.substr(3, 19);
            Ephemeris ephemeris;
            ephemeris.constellation = rinexInfo.constellation;
            ephemeris.satID = rinexInfo.constellation + line.substr(0, 2);
            // convert to GPS time tag
            gtime_t timeAtEpoch = str2gtime(timeStr);
            ephemeris.gtime = timeAtEpoch;
            ephemeris.gpsTime = gtime2gps(timeAtEpoch);
            if (rinexInfo.constellation != 'R') {
                int nLeaps = countLeaps(ephemeris.gpsTime, "utc2gps");
                ephemeris.gpsTime = ephemeris.gpsTime - nLeaps;
            }
            for (int i = 0; i < numberOfValues; ++i) {
                ephemeris.ephemerisData[i] = std::stof(raw.substr(i * 19, 19));
            }
            // GLONASS uses kilometers to report its ephemeris. Convert to meters.
            if (rinexInfo.constellation == 'R') {
                ephemeris.ephemerisData[3] *= 1e3;  // X
                ephemeris.ephemerisData[4] *= 1e3;  // dX
                ephemeris.ephemerisData[5] *= 1e3;  // dX2
                ephemeris.ephemerisData[7] *= 1e3;  // Y
                ephemeris.ephemerisData[8] *= 1e3;  // dY
                ephemeris.ephemerisData[9] *= 1e3;  // dY2
                ephemeris.ephemerisData[11] *= 1e3; // Z
                ephemeris.ephemerisData[12] *= 1e3; // dZ
                ephemeris.ephemerisData[13] *= 1e3; // dZ2
                if (ephemeris.satID < "R10") {
                    ephemeris.satID = "R0" + ephemeris.satID.substr(2);
                }
            }
            // insert the ephemeris to the navigation data
            satEphemeris.push_back(ephemeris);
        }
    }
    else if (rinexInfo.version[0] == '3') {
        // read the header of the RINEX file
        parseRinex3Header(file);
        //save the ionospheric corrections
        if (header.contains("IONOSPHERIC CORR GPSA")) {
            std::vector<float> gpsa;
            std::vector<float> gpsb;
            for (int i = 0; i < 4; i++) {
                std::string gpsaStr = header["IONOSPHERIC CORR GPSA"].substr(1 + i * 12, 11);
                std::string gpsbStr = header["IONOSPHERIC CORR GPSB"].substr(1 + i * 12, 11);
                std::ranges::replace(gpsaStr, 'D', 'E');
                std::ranges::replace(gpsbStr, 'D', 'E');
                gpsa.push_back(std::stof(gpsaStr));
                gpsb.push_back(std::stof(gpsbStr));
            }
            rinexInfo.ionoCorrections["GPSA"] = gpsa;
            rinexInfo.ionoCorrections["GPSB"] = gpsb;
        }
        // read the ephemeris data
        std::string line;
        int numberOfValues = numberOfEphemerisData.at(rinexInfo.constellation);
        // read the raw data for one satellite
        while (std::getline(file, line)) {
            // Read raw data
            std::string raw = line.substr(23, 57); // From column 23 to 79
            for (int i = 0; i < numberOfExtraLines.at(rinexInfo.constellation); ++i) {
                std::string nextLine;
                if (std::getline(file, nextLine)) {
                    raw += nextLine.substr(4, 77);
                }
            }
            std::ranges::replace(raw, 'D', 'E');
            // Parse raw data
            std::string timeStr = line.substr(6, 17);
            Ephemeris ephemeris;
            ephemeris.constellation = rinexInfo.constellation;
            ephemeris.satID = line.substr(0, 3);
            // convert to GPS time tag
            gtime_t timeAtEpoch = str2gtime(timeStr);
            ephemeris.gtime = timeAtEpoch;
            ephemeris.gpsTime = gtime2gps(timeAtEpoch);
            if (rinexInfo.constellation != 'R') {
                int nLeaps = countLeaps(ephemeris.gpsTime, "utc2gps");
                ephemeris.gpsTime = ephemeris.gpsTime - nLeaps;
            }
            for (int i = 0; i < numberOfValues; ++i) {
                ephemeris.ephemerisData[i] = std::stof(raw.substr(i * 19, 19));
            }
            // GLONASS uses kilometers to report its ephemeris. Convert to meters.
            if (rinexInfo.constellation == 'R') {
                ephemeris.ephemerisData[3] *= 1e3;  // X
                ephemeris.ephemerisData[4] *= 1e3;  // dX
                ephemeris.ephemerisData[5] *= 1e3;  // dX2
                ephemeris.ephemerisData[7] *= 1e3;  // Y
                ephemeris.ephemerisData[8] *= 1e3;  // dY
                ephemeris.ephemerisData[9] *= 1e3;  // dY2
                ephemeris.ephemerisData[11] *= 1e3; // Z
                ephemeris.ephemerisData[12] *= 1e3; // dZ
                ephemeris.ephemerisData[13] *= 1e3; // dZ2
                if (ephemeris.satID < "R10") {
                    ephemeris.satID = "R0" + ephemeris.satID.substr(2);
                }
            }
            // insert the ephemeris to the navigation data
            satEphemeris.push_back(ephemeris);
        }
    }
    std::sort(satEphemeris.begin(), satEphemeris.end());
}


// parse the header of RINEX file 2.xx
void RinexReader::parseRinex2Header(std::ifstream &file) {
    std::string line;

    // Read header lines
    while (std::getline(file, line)) {
        // Process lines until "END OF HEADER"
        if (line.find("END OF HEADER") != std::string::npos) {
            break;
        }
        if (line.size() >= 60) {
            std::string key = line.substr(60);
            key = rtrim(key);
            std::string value = line.substr(0, 60);
            value = rtrim(value);
            header.insert({key, value});
        }
    }
}

// parse the header of RINEX file 3.xx
void RinexReader::parseRinex3Header(std::ifstream &file) {
    std::string line;

    // Read header lines
    while (std::getline(file, line)) {
        // Process lines until "END OF HEADER"
        if (line.find("END OF HEADER") != std::string::npos) {
            break;
        }
        if (line.size() >= 60) {
            std::string key = line.substr(60);
            key = rtrim(key);
            std::string value = line.substr(0, 60);
            value = rtrim(value);

            // ionospheric / time system corrections needs more detailed parsing
            if (key == "IONOSPHERIC CORR" || key == "TIME SYSTEM CORR") {
                key += " " + value.substr(0, 4);
                value = value.substr(5);
            }
            header.insert({key, value});
        }
    }
}







