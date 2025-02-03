/* signalFunctions.h
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Header file for signalFunctions.cpp
 * This file contains the necessary functions for signal operation
 */

#pragma once

#include <cmath>

#include <iostream>
#include <vector>
#include <deque>
#include <complex>
#include <fftw3.h>

int positiveMod(int x, int y);

std::vector<std::complex<double> > copySignal(std::vector<std::complex<double>> signal, int nCopy);

std::vector<int> genCAcode(int prn);

std::vector<int> genGalE1(int prn);

std::vector<int> genGloL1(int prn);

std::vector<int> getFreqBin(int start, int end, int space);

std::vector<std::complex<double>> shiftCorrelation(const std::vector<std::complex<double>>& xDelay, const std::vector<std::complex<double>>& xLocal);


