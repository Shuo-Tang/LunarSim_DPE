/* timeFunctions.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for timeFunctions.cpp
 * This file contains the declarations of the time-utility functions for GNSS
 */

#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>


#define GPS_OFFSET 315964800.0 // Offset between Unix epoch and GPS epoch (1970-01-06 00:00:00 UTC)

// time struct defined as RTKLIB
// gtime is always expressed in UTC time
typedef struct {
 time_t time;   // time (s) expressed by standard time_t
 float sec;    // fractional part of the second (under 1 second)
} gtime_t;

constexpr int SECONDS_IN_WEEK = 604800;


int countLeaps(long gpsTime, const std::string& dirFlag);

long utc2gps(double utcTime);

double gps2utc(long gpsTime);

double gtime2utc(const gtime_t& gtime);

double gtime2gps(const gtime_t& gtime);

double timeDiff(gtime_t &t1, gtime_t &t2);

gtime_t str2gtime(const std::string& timeStr);

std::string gtime2str(const gtime_t& gtime);

std::tm str2ctime(const std::string& timeStr);

double weekTow2gps(int& week, double& tow, char constellation);

std::pair<int, double> gps2weekTow(const double& gpsTime);

double checkWeekCross(const double& time);