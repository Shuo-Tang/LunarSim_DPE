# LunarSim_v0_2

LunarSim_v0_2 is a C++ simulation tool for computing satellite positions using GNSS signals.

## 📂 Project Structure
```
LunarSim_v0_2/
│── include/               # Header files
│── src/                   # Source files
│── external/              # External dependencies (Eigen, FFTW, YAML-CPP)
│── test/                  # Unit tests (if applicable)
│── data/                  # Navigation files, test data
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

## 📜 License
MIT License

---

## ⚠️ Note
Ensure `fftw3.dll` is available in the execution path on Windows.
