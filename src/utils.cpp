/**
 * @file utils.cpp
 * @brief Implementation of readParticlesFromFile. 
 *        Minimal I/O.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib> 

#include "utils.h"
#include "config.h"
#include "benchmark.h"

/**
 * @brief read file => lines => parse ptype x y z => create Particle => push in vector
 * @param filename path
 * @param particles out vector
 * @param box_dim dimension
 * @param pBench optional instrumentation
 */
void readParticlesFromFile(const std::string& filename, 
                           std::vector<Particle>& particles, 
                           f64 box_dim, 
                           Benchmark* pBench /*=nullptr*/)
{
    // measure time
    if (pBench) {
        pBench->startTimer("readFile");
    }

    particles.clear();
    particles.reserve(N_PARTICULES_TOTAL);

    std::ifstream inF(filename);
    if (!inF.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }

    // skip first line
    std::string line;
    std::getline(inF, line);

    i32 ptype;
    f64 xx, yy, zz;
    i32 count = 0;
    // approximate ops
    u64 locOps = 0ULL;

    while (std::getline(inF, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        if (!(iss >> ptype >> xx >> yy >> zz)) {
            continue;
        }
        // create Particle => pos= (xx,yy,zz), vel= (0,0,0), mass=1
        std::vector<f64> posV = {xx, yy, zz};
        std::vector<f64> velV = {0.0, 0.0, 0.0};

        Particle tmpPrt(posV, velV, 1.0);
        tmpPrt.setBoxSize(box_dim);
        particles.push_back(tmpPrt);
        count++;

        locOps += 5ULL; // guess: parse+ store

        if (count >= N_PARTICULES_TOTAL) break;
    }

    if ((i32)particles.size() < N_PARTICULES_TOTAL) {
        std::cerr << "Warn: less than " << N_PARTICULES_TOTAL 
                  << " part. read from file.\n";
    }
    inF.close();

    // instrumentation
    if (pBench) {
        pBench->addOps("readFile", locOps);
        pBench->endTimer("readFile");
    }
}
