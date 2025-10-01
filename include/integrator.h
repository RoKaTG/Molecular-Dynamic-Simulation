/**
 * @file integrator.h
 * @brief Declaration of base Integrator class (abstract). 
 *        Derive from this to implement specific integration schemes (e.g. Verlet).
 *
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "system.h"
#include "types.h"
#include "benchmark.h" 

/**
 * @class Integrator
 * @brief Abstract base class for MD integrators. 
 */
class Integrator {
public:
    /**
     * @brief Virtual destructor
     */
    virtual ~Integrator() = default;

    /**
     * @brief Integrate one step. 
     *        Derived classes (like VerletIntegrator) implement the actual logic.
     * @param system The System to integrate
     * @param dt The time step
     * @param pBench 
     */
    virtual void integrate(System& system, f64 dt, Benchmark* pBench=nullptr) = 0;
};

#endif // INTEGRATOR_H
