/* ephemeris.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for Ephemeris.cpp
 * This file contains the necessary includes and declarations for the Ephemeris.cpp file
 * The class helps to store the ephemeris information of a satellite
 */

#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include "timeFunctions.h"

#define MAX_NUM_EPHEMERIS_FIELDS 31

const std::unordered_map<char, int> orderOfConstellations = {
    {'G', 0},
    {'R', 1},
    {'E', 2},
    {'C', 3}
};


class Ephemeris {
public:
    // constructor
    Ephemeris();

    // comparator for sorting
    bool operator<(const Ephemeris& ephemeris) const;

    char constellation;

    std::string satID;

    gtime_t gtime{};

    double gpsTime;

    std::vector<float> ephemerisData;

    /* Ephemeris data vector include different parameters for different constellations
     * GPS: 'G' 29 parameters
     *  - 0: satClockBias       - 1: satClockDrift      - 2: satClockDriftRate
     *  - 3: IODE               - 4: Crs                - 5: DeltaN                 - 6: M0
     *  - 7: Cuc                - 8: eccentricity       - 9: Cus                    - 10: sqrtA
     *  - 11: Toe               - 12: Cic               - 13: Omega0                - 14: Cis
     *  - 15: I0                - 16: Crc               - 17: Omega                 - 18: OmegaDot
     *  - 19: IDOT              - 20: codesL2           - 21: GPSWeek               - 22: L2PFlag
     *  - 23: satAcc            - 24: health            - 25: TGD                   - 26: IODC
     *  - 27: transTime         - 28: fitIntvl
     * GLONASS: 'R' 15 parameters
     * - 0: satClockBias        - 1: satRelFreqBias     - 2: messageFrameTime
     * - 3: X                   - 4: dX                 - 5: dX2                    - 6: health
     * - 7: Y                   - 8: dY                 - 9: dY2                    - 10: freqNum
     * - 11: Z                  - 12: dZ                - 13: dZ2                   - 14: ageOpInfo
     *
     * Galileo: 'E' 27 parameters   (Rinex 2.xx)
     * - 0: satClockBias        - 1: satClockDrift      - 2: satClockDriftRate
     * - 3: IODNav              - 4: Crs                - 5: DeltaN                 - 6: M0
     * - 7: Cuc                 - 8: eccentricity       - 9: Cus                    - 10: sqrtA
     * - 11: Toe                - 12: Cic               - 13: Omega0                - 14: Cis
     * - 15: I0                 - 16: Crc               - 17: Omega                 - 18: OmegaDot
     * - 19: IDOT               - 20: DataSrc           - 21: GALWeek               - 22: SISA
     * - 23: health             - 24: BGDe5a            - 25: BGDe5b                - 26: transTime
     *
     * Galileo: 'E' 28 parameters   (Rinex 3.xx)
     * - 0: satClockBias        - 1: satClockDrift      - 2: satClockDriftRate
     * - 3: IODNav              - 4: Crs                - 5: DeltaN                 - 6: M0
     * - 7: Cuc                 - 8: eccentricity       - 9: Cus                    - 10: sqrtA
     * - 11: Toe                - 12: Cic               - 13: Omega0                - 14: Cis
     * - 15: I0                 - 16: Crc               - 17: Omega                 - 18: OmegaDot
     * - 19: IDOT               - 20: DataSrc           - 21: GALWeek               - 22: spare0
     * - 23: SISA               - 24: health            - 25: BGDe5a                - 26: BGDe5b
     * - 27: transTime
     *
     * Beidou: 32 parameters        (RINEX 3.xx)
     * -  0: SatClockBias       -  1: SatClockDrift     -  2: SatClockDriftRate
     * -  3: AODE               -  4: Crs               -  5: DeltaN                -  6: M0
     * -  7: Cuc                -  8: Eccentricity      -  9: Cus                   - 10: sqrtA
     * - 11: Toe                - 12: Cic               - 13: Omega0                - 14: Cis
     * - 15: I0                 - 16: Crc               - 17: omega                 - 18: OmegaDot
     * - 19: IDOT               - 20: spare0            - 21: BDTWeek               - 22: spare1
     * - 23: satAcc             - 24: SatH1             - 25: TGD1                  - 26: TGD2
     * - 27: TransTime          - 28: AODC              - 29: spare2                - 30: spare3
     */




};

