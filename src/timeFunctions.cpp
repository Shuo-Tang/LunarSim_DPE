/* timeFunctions.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for timeFunctions.h
 * This file contains the necessary functions for GNSS time conversion and computation
 */

#include "timeFunctions.h"

// Define GPS leap seconds
std::vector<long> getLeaps() {
    return {46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
            346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
            599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017};
}

// Count number of leap seconds that have passed
int countLeaps(long gpsTime, const std::string& dirFlag) {
    std::vector<long> leaps = getLeaps();
    int nLeaps = 0;  // number of leap seconds prior to gpsTime

    for (size_t i = 0; i < leaps.size(); i++) {
        if (dirFlag == "utc2gps") {
            if (gpsTime >= leaps[i] - i) {
                nLeaps++;
            }
        } else if (dirFlag == "gps2utc") {
            if (gpsTime >= leaps[i]) {
                nLeaps++;
            }
        } else {
            std::cerr << "ERROR Invalid Flag!" << std::endl;
            return -1;
        }
    }
    return nLeaps;
}

long utc2gps(double utcTime) {
    // Add offset in seconds
    long gpsTime = static_cast<long>(utcTime) - GPS_OFFSET;
    int nLeaps = countLeaps(gpsTime, "utc2gps");
    gpsTime = gpsTime + nLeaps;
    return gpsTime;
}

// Convert GPS Time to Unix Time
double gps2utc(long gpsTime) {
    // Add offset in seconds
    double utcTime = static_cast<double>(gpsTime) + GPS_OFFSET;
    int nLeaps = countLeaps(gpsTime, "gps2utc");
    utcTime = utcTime - nLeaps;
    return utcTime;
}

// Convert gtime format to UTC time
double gtime2utc(const gtime_t& gtime) {
    return gtime.time + static_cast<double>(gtime.sec);
}

// Convert gtime format to GPS time
double gtime2gps(const gtime_t& gtime) {
    return utc2gps(gtime2utc(gtime));
}

double timeDiff(gtime_t &t1, gtime_t &t2) {
    return difftime(t1.time,t2.time) + t1.sec - t2.sec;
}

// Convert string to gtime format
gtime_t str2gtime(const std::string& timeStr) {
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    // Parse the input string
    const char* str = timeStr.c_str();
    char* end;
    int year = static_cast<int>(std::strtol(str, &end, 10));
    year += (year < 80) ? 2000 : 1900;
    int month = static_cast<int>(std::strtol(end, &end, 10));
    int day = static_cast<int>(std::strtol(end, &end, 10));
    int hour = static_cast<int>(std::strtol(end, &end, 10));
    int minute = static_cast<int>(std::strtol(end, &end, 10));
    auto second = static_cast<float>(std::strtod(end, nullptr));

    // Calculate the number of days since 1970-01-01
    int days = (year-1970)*365 + (year-1969)/4 + doy[month-1] + day - 2 + (year%4==0&&month>=3?1:0);
    gtime_t gtime;
    gtime.time = static_cast<time_t>(days)*86400 + hour*3600 + minute*60 + static_cast<int>(second);
    gtime.sec = second - static_cast<int>(second);

    return gtime;
}

std::string gtime2str(const gtime_t& gtime) {
    std::tm* utcTime = std::gmtime(&gtime.time);
    char buffer[20];
    std::strftime(buffer, sizeof(buffer), "%Y %m %d %H %M %S", utcTime);
    std::string timeStr = std::string(buffer) + std::to_string(gtime.sec);
    return timeStr;
}

std::tm str2ctime(const std::string& timeStr) {
    // Parse the input string
    const char* str = timeStr.c_str();
    char* end;
    int year = static_cast<int>(std::strtol(str, &end, 10));
    int month = static_cast<int>(std::strtol(end, &end, 10));
    int day = static_cast<int>(std::strtol(end, &end, 10));
    int hour = static_cast<int>(std::strtol(end, &end, 10));
    int minute = static_cast<int>(std::strtol(end, &end, 10));
    auto second = static_cast<float>(std::strtod(end, nullptr));

    // Create a tm struct and set the values
    std::tm tm_time = {};
    year += (year < 80) ? 2000 : 1900;
    tm_time.tm_year = year - 1900;  // tm_year is years since 1900
    tm_time.tm_mon = month - 1;     // tm_mon is 0-based (January is 0)
    tm_time.tm_mday = day;
    tm_time.tm_hour = hour;
    tm_time.tm_min = minute;
    tm_time.tm_sec = static_cast<int>(second);

    return tm_time;
}

double weekTow2gps(int& week, double& tow, char constellation) {
    double time = week * 7 * 86400.0 + tow;
    if (constellation == 'C') { // BeiDou system
        const int GPS_BDS_WEEK = 1356; // GPS week number on 1st January 2006
        time = (GPS_BDS_WEEK + week) * 7 * 86400.0 + tow;
    }
    return time;
}

std::pair<int, double> gps2weekTow(const double& gpsTime) {
    int weekNumber = static_cast<int>(gpsTime / SECONDS_IN_WEEK);
    double secondOfWeek = gpsTime - (weekNumber * SECONDS_IN_WEEK);
    return std::make_pair(weekNumber, secondOfWeek);
}

double checkWeekCross(const double& time) {
    constexpr double halfWeekSec = 302400.0;
    if (time > halfWeekSec) {
        return time - 2 * halfWeekSec;
    } else if (time < - halfWeekSec) {
        return time + 2 * halfWeekSec;
    }
    return time;
}


