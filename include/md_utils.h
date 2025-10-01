/**
 * @file md_utils.h
 * @brief Various utility functions for simulation: kineticEnergy, temperature, velocity init, PDB saving.
 *
 * We can pass a Benchmark* if we want instrumentation (time, ops).
 */

#ifndef MD_UTILS_H
#define MD_UTILS_H

#include "system.h"
#include "types.h"
#include <string>
#include "benchmark.h" // optional

/**
 * @brief Compute total kinetic energy in the system:
 *        E_kin = 1/(2*convForce) * sum_i( mass_i * v^2 ).
 * @param sys The System (we read Particles).
 * @param pBench Optional instrumentation
 * @return f64 kinetic energy
 */
f64 computeKineticEnergy(const System& sys, Benchmark* pBench=nullptr);

/**
 * @brief Compute temperature from E_kin. 
 *        T = (1 / (Ndl * const_R)) * Ecin, with Ndl=3*N - 3.
 * @param eCin The kinetic energy
 * @param N number of particles
 * @return temperature
 */
f64 computeTemperature(f64 eCin, u32 N);

/**
 * @brief Initialize velocities:
 *        1) random uniform [-0.5,0.5]
 *        2) subtract center-of-mass velocity
 *        3) compute Ecin => rescale => T0
 * @param sys The System
 * @param T0 The target temperature
 * @param pBench Optional instrumentation
 */
void initializeVelocities(System& sys, f64 T0, Benchmark* pBench=nullptr);

/**
 * @brief Save the current state to PDB file (very naive format).
 * @param sys The System
 * @param iteration The iteration number
 * @param outFile The PDB file path
 * @param pBench Optional instrumentation
 */
void savePDBFrameExact(const System& sys, u32 iteration, 
                       const std::string& outFile, 
                       Benchmark* pBench=nullptr);

#endif // MD_UTILS_H
