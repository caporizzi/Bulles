#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.vtk> <output.dat>" << std::endl;
        return 1;
    }

    std::string vtkFile = argv[1];
    std::string outFile = argv[2];

    std::ifstream file(vtkFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the .vtk file: " << vtkFile << std::endl;
        return 1;
    }

    int nx = 0, ny = 0, nz = 0;
    int pointDataCount = -1, cellDataCount = -1;
    bool foundDimensions = false, foundScalars = false;
    std::string line;
    std::vector<float> scalars;

    while (std::getline(file, line)) {
        if (line.rfind("DIMENSIONS", 0) == 0) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> nx >> ny >> nz;
            foundDimensions = true;
        }
        else if (line.rfind("POINT_DATA", 0) == 0) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> pointDataCount;
        }
        else if (line.rfind("CELL_DATA", 0) == 0) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> cellDataCount;
        }
        else if (line.rfind("LOOKUP_TABLE", 0) == 0) {
            foundScalars = true;
            continue;
        }
        else if (foundScalars) {
            std::istringstream iss(line);
            float val;
            while (iss >> val) {
                scalars.push_back(val);
            }
        }
    }
    file.close();

    if (!foundDimensions || scalars.empty()) {
        std::cerr << "Error: Missing DIMENSIONS or SCALARS in .vtk file!" << std::endl;
        return 1;
    }

    // Decide whether data is POINT_DATA or CELL_DATA
    if (pointDataCount > 0) {
        // should match nx*ny*nz
        if ((int)scalars.size() != pointDataCount) {
            std::cerr << "Warning: POINT_DATA count mismatch. Got " << scalars.size()
                      << " expected " << pointDataCount << std::endl;
        }
    }
    else if (cellDataCount > 0) {
        // adjust to cell-centered grid
        nx = std::max(1, nx);
        ny = std::max(1, ny);
        nz = std::max(1, nz);
        if ((int)scalars.size() != cellDataCount) {
            std::cerr << "Warning: CELL_DATA count mismatch. Got " << scalars.size()
                      << " expected " << cellDataCount << std::endl;
        }
    }
    else {
        std::cerr << "Error: No POINT_DATA or CELL_DATA found in file!" << std::endl;
        return 1;
    }

    // Final check
    if ((int)scalars.size() != nx * ny * nz) {
        std::cerr << "Error: scalar count (" << scalars.size()
                  << ") does not match expected size (" << nx*ny*nz << ")" << std::endl;
        return 1;
    }

    // Write binary .dat
    std::ofstream out(outFile, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file: " << outFile << std::endl;
        return 1;
    }

    out.write((char*)&nx, sizeof(nx));
    out.write((char*)&ny, sizeof(ny));
    out.write((char*)&nz, sizeof(nz));
    out.write((char*)scalars.data(), scalars.size() * sizeof(float));
    out.close();

    std::cout << "Data saved to " << outFile
              << " (" << nx << "x" << ny << "x" << nz
              << ", total " << scalars.size() << " floats)" << std::endl;

    return 0;
}