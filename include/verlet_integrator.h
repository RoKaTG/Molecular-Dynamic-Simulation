/**
 * @file verlet_integrator.h
 * @brief Declaration of VerletIntegrator class: 
 *        Integrator derived from Integrator base. 
 */

#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

#include "integrator.h"
#include "benchmark.h" // for optional pointer

/**
 * @class VerletIntegrator
 * @brief Velocity-Verlet scheme integrator: update pos => recalc force => update vel.
 */
class VerletIntegrator : public Integrator {
public:
    /**
     * @brief Integrate one step with velocity-Verlet:
     *        1) store old a(t), update x(t+dt)
     *        2) compute new force => a(t+dt)
     *        3) update v(t+dt)
     * @param system The system to integrate
     * @param dt Timestep
     * @param pBench Optional benchmark
     */
    void integrate(System& system, f64 dt, Benchmark* pBench=nullptr) override;
};

#endif // VERLET_INTEGRATOR_H
