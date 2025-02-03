# LunarSim_v0_2

LunarSim_v0_2 is a C++ simulation tool for computing satellite positions using GNSS signals.

## 📂 Project Structure
```
LunarSim_v0_2/
│── include/               # Header files
│── src/                   # Source files
│── external/              # External dependencies (Eigen, FFTW, YAML-CPP)
│── importData/            # Navigation files and .yaml settings
│── exportData/            # Inforamtion of the Moon and satellites, output of the signal
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
   git clone https://github.com/yourusername/LunarSim_v0_2.git
   cd LunarSim_v0_2
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
   ./LunarSim_v0_2
   ```
### Programs
Two versions are provided for LunarSim_v0_2: 
	- The c++ codes are used to generate and sample the GNSS received signal on the Moon.
	- The MATLAB codes are used to generate and process the lunar signals in real-time without saving the signal records. You can implement your own positioning algorithms in baseband-signal level and get the positioning result.
	  Here, we provide the two-step (2SP) and Direct Position Estimation solutions.\
	- !!! The MATLAB script **`moon_position.m`** must be run to generate the Moon's position in ECI coordinates before any code works.

## 📚 Citation 
- If you use this project in your research, please cite:
     ```
     Tang, S., Li, H. and Closas, P.,
     Direct Position Estimation Framework for Lunar Positioning, Navigation, and Timing. TBD...
     ```
    

## 📜 License
MIT License

---

## ⚠️ Note
- You may modify the cmake file in your way to link the library and all files, but the followinfg packages of dependencies have to be linked to the project.
	- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
	- [eigen-3.4.0](https://eigen.tuxfamily.org/index.php?title=Main_Page)
	- [fftw-3.3.5](https://www.fftw.org/install/windows.html)
- You can use diffenent versions of the above packages, but it may need modifying the ``#include<>`` in the codes.
