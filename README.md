# GNSSLunarSim_v0_2

GNSSLunarSim_v0_2 is a C++/MATLAB simulation tool designed to generating the GNSS received baseband signal on the Moon for lunar positioning algorithm validation.

## 📂 Project Structure
```
GNSSLunarSim_v0_2/
│── include/               # Header files
│── src/                   # Source files
│── external/              # External dependencies (Eigen, FFTW, YAML-CPP)
│── importData/            # Navigation files and YAML configuration settings
│── exportData/            # Output data including lunar and satellite information
│── CMakeLists.txt
│── README.md
│── .gitignore
```
## 🚀 Build Instructions

### Prerequisites
- CMake 3.29+
- C++20 compatible compiler (GCC, Clang, MSVC)
- Git installed for YAML-CPP integration

### Steps
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/GNSSLunarSim_v0_2.git
   cd GNSSLunarSim_v0_2
   ```
2. Create a build directory and configure:
   ```sh
   mkdir build && cd build
   cmake ..
   ```
3. Compile:
   ```sh
   cmake --build .
   ```
4. Run:
   ```sh
   ./GNSSLunarSim_v0_2
   ```

## 📌 Program Versions
GNSSLunarSim_v0_2 is available in two versions:
- **C++ Implementation**: Generates and samples the GNSS received signal on the Moon for analysis.
- **MATLAB Implementation**: Provides real-time signal generation and processing for lunar GNSS signals without saving signal records. This version allows users to develop and test positioning algorithms at the baseband signal level.
  - Our implementation includes **Two-Step Positioning (2SP)** and **Direct Position Estimation (DPE)** solutions.
  - **Important:** Before executing any code, the MATLAB script **`moon_position.m`** must be run to compute the Moon's position in Earth-Centered Inertial (ECI) coordinates.

## 📚 Citation
If you use this project in your research, please cite:
```
Tang, S., Li, H., and Closas, P.,
Direct Position Estimation Framework for Lunar Positioning, Navigation, and Timing. TBD...
```

## 📜 License
This project is licensed under the MIT License.

## ⚠️ Important Notes
- The CMake configuration can be modified as needed, but the following dependencies must be properly linked:
  - [yaml-cpp](https://github.com/jbeder/yaml-cpp)
  - [Eigen 3.4.0](https://eigen.tuxfamily.org/index.php?title=Main_Page)
  - [FFTW 3.3.5](https://www.fftw.org/install/windows.html)
- Different versions of these dependencies may be used, but adjustments to `#include<>` directives in the source code might be required to ensure compatibility.
