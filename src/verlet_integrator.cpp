/**
 * @file verlet_integrator.cpp
 * @brief Implementation of VerletIntegrator class. 
 *        Performs a velocity-Verlet integration step: 
 *        (1) store old acceleration, update positions, 
 *        (2) compute new forces => new acceleration,
 *        (3) update velocities.
 *
 */

#include <omp.h>
#include <cmath>
#include <vector>

#include "verlet_integrator.h"
#include "system.h"
#include "config.h"     
#include "benchmark.h"

/**
 * @brief Integrate one step using velocity-Verlet scheme:
 *        x(t+dt)= x(t)+v(t)*dt+0.5*a_old*convFrc/m *dt^2
 *        Recompute force => a_new
 *        v(t+dt)= v(t)+0.5*(a_old+a_new)*convFrc *dt
 * 
 * @param system The System containing Particles, neighbor list, etc.
 * @param dt Timestep in (e.g.) femtoseconds
 * @param pBench Optional pointer to Benchmark for timing/ops
 */
void VerletIntegrator::integrate(System& system, f64 dt, Benchmark* pBench /*=nullptr*/)
{
    if (pBench) {
        pBench->startTimer("integrateVerlet");
    }

    auto& vecParts = system.getParticles();
    const f64 convFrc = CONVERSION_FORCE; // shorten

    // I count approximate ops
    // step(1) => ~ 10 ops * N
    // step(2) => done inside computeForcesWithVerlet => separate label
    // step(3) => ~ 5 ops * N
    u64 locOps = 0ULL;

    // 0) store old acceleration => a_old = a(t)
    //    a(t)=F(t)/m => already computed if not the first step
    std::vector<std::vector<f64>> oldAcc(vecParts.size(), std::vector<f64>(3,0.0));

#pragma omp parallel for
    for (u32 idxP = 0; idxP < (u32)vecParts.size(); idxP++) {
        oldAcc[idxP] = vecParts[idxP].getAcceleration();
    }

    // 1) update position => x(t+dt)= x(t)+v(t)*dt+0.5*a_old*convFrc/m*(dt^2)
#pragma omp parallel for
    for (u32 idxP = 0; idxP < (u32)vecParts.size(); ++idxP) {
        auto posTmp = vecParts[idxP].getPosition();
        auto velTmp = vecParts[idxP].getVelocity();
        const auto& accTmp = oldAcc[idxP]; 
        f64 mVal = vecParts[idxP].getMass();

        std::vector<f64> newPos(3);
        for (u32 jj = 0; jj < 3; ++jj) {
            newPos[jj] = posTmp[jj]
                         + velTmp[jj]*dt
                         + 0.5*(accTmp[jj]*convFrc)/mVal*(dt*dt);
        }
        vecParts[idxP].setPosition(newPos);
    }
    // approx ops => ~10*N
    locOps += 10ULL * (u64)vecParts.size();

    // 2) recompute forces => a(t+dt) => system calls neighbor approach
    //    internally, we do pBench->startTimer("forceCalc") => etc.
    //    So I won't double count ops here. 
    system.computeForcesWithVerlet(pBench);

    // 2b) store new acceleration => a_new = F(t+dt)/m
    std::vector<std::vector<f64>> newAcc(vecParts.size(), std::vector<f64>(3,0.0));

#pragma omp parallel for
    for (u32 idxP = 0; idxP < (u32)vecParts.size(); ++idxP) {
        auto frcTmp = vecParts[idxP].getForce();
        f64 mVal = vecParts[idxP].getMass();
        for (u32 jj=0; jj<3; jj++){
            newAcc[idxP][jj] = frcTmp[jj]/mVal;
        }
        // store in particle
        vecParts[idxP].setAcceleration(newAcc[idxP]);
    }

    // 3) update velocity => v(t+dt)= v(t)+0.5*(a_old+a_new)*convFrc*dt
#pragma omp parallel for
    for (u32 idxP = 0; idxP < (u32)vecParts.size(); ++idxP) {
        auto velTmp = vecParts[idxP].getVelocity();
        for (u32 jj=0; jj<3; jj++){
            velTmp[jj] += 0.5*(oldAcc[idxP][jj] + newAcc[idxP][jj])*convFrc*dt;
        }
        vecParts[idxP].setVelocity(velTmp);
    }
    // approx ops => ~5*N
    locOps += 5ULL * (u64)vecParts.size();

    // finalize instrumentation
    if (pBench) {
        pBench->addOps("integrateVerlet", locOps);
        pBench->endTimer("integrateVerlet");
    }
}
