/**
 * @file utils.h
 * @brief Basic I/O routines for reading Particles from file.
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include "types.h"
#include "particle.h"
#include "benchmark.h" 

/**
 * @brief Read Particles from a file => store in 'particles'.
 *        Format: skip first line, then lines with "ptype x y z".
 *        If boxDim is given, set for each Particle.
 * @param filename path
 * @param particles out vector
 * @param box_dim the box dimension
 * @param pBench optional instrumentation
 */
void readParticlesFromFile(const std::string& filename, 
                           std::vector<Particle>& particles, 
                           f64 box_dim, 
                           Benchmark* pBench=nullptr);

#endif // UTILS_H
