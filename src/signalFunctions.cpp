/* signalFunctions.cpp
 * Serve for LunarSimulator From Shuo Tang, PhD student in SPIRAL, NEU 12/05/2024
 * Source file for signalFunctions.h
 * This file contains the necessary functions for signal operation
 */

#include "signalFunctions.h"

int positiveMod(int x, int y) {
    int ans = fmod(x, y);
    return ans >= 0 ? ans : ans + y;
}

std::vector<std::complex<double> > copySignal(std::vector<std::complex<double>> signal, int nCopy) {
    std::vector<std::complex<double>> signalCopy;
    signalCopy.reserve(nCopy * signal.size());
    for (int i = 0; i < nCopy; ++i) {
        signalCopy.insert(signalCopy.end(), signal.begin(), signal.end());
    }
    return signalCopy;
}

std::vector<int> genCAcode(int prn) {
    // Define the G2 shifts for each PRN
    const std::vector<int> g2s = {5, 6, 7, 8, 17, 18, 139, 140, 141, 251,
                                  252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
                                  473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
                                  861, 862, 145, 175, 52, 21, 237, 235, 886, 657,
                                  634, 762, 355, 1012, 176, 603, 130, 359, 595, 68,
                                  386};

    // Check if the PRN is valid
    if (prn < 1 || prn > g2s.size()) {
        throw std::invalid_argument("Invalid PRN number. Must be between 1 and " + std::to_string(g2s.size()));
    }

    // Get the G2 shift for the given PRN
    int g2shift = g2s[prn - 1];

    // Generate G1 code
    std::vector<int> g1(1023, -1);  // Initialize G1 with -1
    std::vector<int> reg1(10, -1);  // Shift register for G1

    for (int i = 0; i < 1023; i++) {
        g1[i] = reg1[9];
        int saveBit = reg1[2] * reg1[9];
        for (int j = 9; j > 0; j--) {
            reg1[j] = reg1[j - 1];
        }saveBit;
        reg1[0] = saveBit;
    }

    // Generate G2 code
    std::vector<int> g2(1023, -1);  // Initialize G2 with -1
    std::vector<int> reg2(10, -1);  // Shift register for G2

    for (int i = 0; i < 1023; i++) {
        g2[i] = reg2[9];
        int saveBit = reg2[1] * reg2[2] * reg2[5] * reg2[7] * reg2[8] * reg2[9];
        for (int j = 9; j > 0; j--) {
            reg2[j] = reg2[j - 1];
        }
        reg2[0] = saveBit;
    }

    // Shift G2 code
    std::deque<int> g2Deque(g2.begin(), g2.end());
    for (int i = 0; i < g2shift; i++) {
        g2Deque.push_front(g2Deque.back());
        g2Deque.pop_back();
    }
    std::vector<int> g2Shifted(g2Deque.begin(), g2Deque.end());

    // Form C/A code by multiplying G1 and G2
    std::vector<int> CAcode(1023);
    for (int i = 0; i < 1023; i++) {
        CAcode[i] = -(g1[i] * g2Shifted[i]);
    }

    return CAcode;
}

std::vector<int> genGalE1(int prn) {
    // Check if the PRN is valid
    if (prn < 1 || prn > 50) { // Galileo PRN range is typically 1-50
        throw std::invalid_argument("Invalid PRN number. Must be between 1 and 50.");
    }

    // Initialize registers for G1 and G2
    std::vector<int> g1(4092, -1); // Galileo E1 codes are 4092 chips long
    std::vector<int> g2(4092, -1);
    std::vector<int> reg1(13, -1); // 13-stage shift register for G1
    std::vector<int> reg2(13, -1); // 13-stage shift register for G2

    // Generate G1 code
    for (int i = 0; i < 4092; i++) {
        g1[i] = reg1[12];
        int saveBit = reg1[0] * reg1[3] * reg1[7] * reg1[12]; // Feedback taps
        for (int j = 12; j > 0; j--) {
            reg1[j] = reg1[j - 1];
        }
        reg1[0] = saveBit;
    }

    // Generate G2 code
    for (int i = 0; i < 4092; i++) {
        g2[i] = reg2[12];
        int saveBit = reg2[0] * reg2[1] * reg2[8] * reg2[12]; // Feedback taps
        for (int j = 12; j > 0; j--) {
            reg2[j] = reg2[j - 1];
        }
        reg2[0] = saveBit;
    }

    // Combine G1 and G2 to create E1 code
    std::vector<int> e1Code(4092);
    for (int i = 0; i < 4092; i++) {
        e1Code[i] = g1[i] * g2[(i + prn) % 4092]; // PRN determines shift
    }

    return e1Code;
}

std::vector<int> genGloL1(int prn) {
    // GLONASS PRN codes are 511 chips long
    const int codeLength = 511;

    // GLONASS PRN numbers typically range from 1 to 24
    if (prn < 1 || prn > 24) {
        throw std::invalid_argument("Invalid PRN number. Must be between 1 and 24.");
    }

    // Initialize shift registers
    std::vector<int> g1(codeLength, -1);
    std::vector<int> reg1(9, -1); // 9-stage shift register for G1

    // Generate G1 code (simple feedback shift register)
    for (int i = 0; i < codeLength; i++) {
        g1[i] = reg1[8]; // Output from the last register
        int saveBit = reg1[0] * reg1[5] * reg1[8]; // Feedback taps
        for (int j = 8; j > 0; j--) {
            reg1[j] = reg1[j - 1];
        }
        reg1[0] = saveBit;
    }

    // Return the generated GLONASS code (same for all satellites)
    return g1;
}

std::vector<int> getFreqBin(int start, int end, int space) {
    std::vector<int> freqBin;
    for (int i = start; i <= end; i += space) {
        freqBin.push_back(i);
    }
    return freqBin;
}

// Function to compute shift correlation
std::vector<std::complex<double>> shiftCorrelation(const std::vector<std::complex<double>>& xDelay, const std::vector<std::complex<double>>& xLocal) {
    int N = xDelay.size();

    if (N != xLocal.size()) {
        throw std::invalid_argument("x_delay and x_local must have the same size");
    }

    // Allocate FFTW input/output arrays
    fftw_complex *fft_x_delay = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *fft_x_local = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *ifft_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Prepare FFT plans
    fftw_plan fft_plan_delay = fftw_plan_dft_1d(N, fft_x_delay, fft_x_delay, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan_local = fftw_plan_dft_1d(N, fft_x_local, fft_x_local, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan ifft_plan = fftw_plan_dft_1d(N, ifft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Fill input arrays
    for (int i = 0; i < N; ++i) {
        fft_x_delay[i][0] = xDelay[i].real(); // Real part
        fft_x_delay[i][1] = xDelay[i].imag(); // Imaginary part
        fft_x_local[i][0] = xLocal[i].real(); // Real part
        fft_x_local[i][1] = xLocal[i].imag(); // Imaginary part
    }

    // Execute FFT for x_delay and x_local
    fftw_execute(fft_plan_delay);
    fftw_execute(fft_plan_local);

    // Perform element-wise multiplication with conjugate
    for (int i = 0; i < N; ++i) {
        double real_local = fft_x_local[i][0];
        double imag_local = fft_x_local[i][1];
        double real_delay = fft_x_delay[i][0];
        double imag_delay = fft_x_delay[i][1];

        // Conjugate of fft_x_local
        double conj_real_local = real_local;
        double conj_imag_local = -imag_local;

        // Complex multiplication
        ifft_result[i][0] = real_delay * conj_real_local - imag_delay * conj_imag_local; // Real part
        ifft_result[i][1] = real_delay * conj_imag_local + imag_delay * conj_real_local; // Imaginary part
    }

    // Perform IFFT
    fftw_execute(ifft_plan);

    // Normalize and store results
    std::vector<std::complex<double>> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = std::complex<double>(ifft_result[i][0] / N, ifft_result[i][1] / N);
    }

    // Cleanup FFTW resources
    fftw_destroy_plan(fft_plan_delay);
    fftw_destroy_plan(fft_plan_local);
    fftw_destroy_plan(ifft_plan);
    fftw_free(fft_x_delay);
    fftw_free(fft_x_local);
    fftw_free(ifft_result);

    return result;
}

