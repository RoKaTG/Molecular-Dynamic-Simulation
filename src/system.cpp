/**
 * @file system.cpp
 * @brief Implementation of System class (particle container, neighbor building, forces, energies).
 *        Optimized & commented for better performance. Includes calls to Benchmark.
 * - I use a triple sum (i_sym, i, j) for periodic images (27 offsets).
 * - Lennard-Jones potential => short-range cut at R_CUT.
 * - For energy: E_tot = sum(kinetic) + sum(potential).
 * - For forces: F_ij from LJ formula => if (r < R_CUT).
 */

#include <omp.h>
#include <cmath>
#include <array>
#include <vector>
#include <iostream>

#include "system.h"
#include "types.h"
#include "benchmark.h"

/**
 * @brief Constructor for System. Potential is nullptr, boxDim=0.0 by default.
 * 
 * This class aggregates:
 *  - A list of Particles (vecParts)
 *  - A list of symmetry offsets (vecSym)
 *  - A neighbor list (vecNbrList) used in Verlet approach
 *  - A pointer to Potential
 *  - The box dimension (boxDim)
 */
System::System()
  : m_potential(nullptr),
    boxDim(0.0)
{}

/**
 * @brief Build symmetry offset vectors (e.g. -L, 0, +L) for periodic images.
 *        Fills vecSym with N_SYM=27 offsets.
 * 
 * I do a small triple loop for nx, ny, nz in {-1,0,1}, 
 * storing each offset in vecSym. 
 */
void System::initSymVectors() {
    vecSym.clear();
    vecSym.reserve(N_SYM);

    // nx, ny, nz in {-1,0,1}
    i32 offsets[3] = {-1, 0, 1};
    for (i32 offX : offsets) {
        for (i32 offY : offsets) {
            for (i32 offZ : offsets) {
                std::array<f64, 3> tmpArr = {
                    offX * L, 
                    offY * L, 
                    offZ * L
                };
                vecSym.push_back(tmpArr);
            }
        }
    }
}

/**
 * @brief Return the const ref to the symmetry offset vectors (27 images).
 * @return const std::vector<std::array<f64,3>>&
 */
const std::vector<std::array<f64,3>>& System::getSymVectors() const {
    return vecSym;
}

/**
 * @brief Add a new Particle to the system. 
 *        Sets boxDim for that particle too (for periodic usage).
 * @param inPart The incoming Particle.
 */
void System::addParticle(const Particle& inPart) {
    Particle tmpP = inPart;
    tmpP.setBoxSize(boxDim);
    vecParts.push_back(tmpP);
}

/**
 * @brief Assign the Potential pointer (e.g. Lennard-Jones).
 * @param pot Shared ptr to Potential object.
 */
void System::setPotential(std::shared_ptr<Potential> pot) {
    m_potential = pot;
}

/**
 * @brief Return reference to the internal vector of Particles.
 */
std::vector<Particle>& System::getParticles() {
    return vecParts;
}

/**
 * @brief Return const reference to the internal vector of Particles.
 */
const std::vector<Particle>& System::getParticles() const {
    return vecParts;
}

/**
 * @brief Set the box dimension for each Particle in the system.
 * @param newBox The new box dimension.
 */
void System::setBoxSize(f64 newBox) {
    boxDim = newBox;
    for (auto& pp : vecParts) {
        pp.setBoxSize(boxDim);
    }
}

/**
 * @brief Get current box dimension.
 * @return f64
 */
f64 System::getBoxSize() const {
    return boxDim;
}

/**
 * @brief Build neighbor list in O(N^2 * N_SYM).
 *        For each i_sym, for each (i<j), check distance < R_CUT^2, store (j, i_sym).
 * 
 *        loop for i_sym in [0..N_SYM-1]:
 *          loop for i in [0..N-1]:
 *            loop for j in [i+1..N-1]:
 *              => check (r2 < R_CUT^2)
 *              => if true => neighborList[i] += (j, i_sym)
 * 
 *        This approach is naive but straightforward.
 *
 * @param pBench Optional pointer to Benchmark for timing/ops instrumentation.
 */
void System::buildNeighborList(Benchmark* pBench /*=nullptr*/)
{
    // (A) Optionally measure time for building neighbors
    if (pBench) {
        pBench->startTimer("buildNbr");
    }

    // (B) Clear old list
    vecNbrList.clear();
    vecNbrList.resize(vecParts.size());

    // I'll count ops for approximate flops in neighbor building
    // For each pair => dx/dy/dz => ~ 6 ops + compare => let's guess 10 ops total
    u64 locOps = 0ULL; 

#pragma omp parallel
    {
        // localVec => to reduce dynamic memory usage in parallel
        std::vector<std::pair<u32,u32>> localVec;

#pragma omp for schedule(dynamic)
        for (i32 kkSym = 0; kkSym < N_SYM; kkSym++) {
            for (u64 idxI = 0; idxI < vecParts.size(); ++idxI) {
                const auto& posI = vecParts[idxI].getPosition();

                for (u64 idxJ = idxI + 1; idxJ < vecParts.size(); ++idxJ) {
                    const auto& posJ = vecParts[idxJ].getPosition();

                    f64 xLoc = posJ[0] + vecSym[kkSym][0];
                    f64 yLoc = posJ[1] + vecSym[kkSym][1];
                    f64 zLoc = posJ[2] + vecSym[kkSym][2];

                    f64 dX = posI[0] - xLoc;
                    f64 dY = posI[1] - yLoc;
                    f64 dZ = posI[2] - zLoc;

                    f64 dist2 = dX*dX + dY*dY + dZ*dZ;
                    locOps += 10ULL; // rough

                    if (dist2 < (R_CUT*R_CUT)) {
                        localVec.emplace_back( (u32)idxJ, (u32)kkSym );
                    }
                }

                // push localVec into vecNbrList[idxI]
                if (!localVec.empty()) {
#pragma omp critical
                    {
                        auto& refVec = vecNbrList[idxI];
                        refVec.insert(refVec.end(), localVec.begin(), localVec.end());
                    }
                    localVec.clear();
                }
            }
        }
    } // end omp parallel

    // (C) register ops in benchmark
    if (pBench) {
        pBench->addOps("buildNbr", locOps);
        pBench->endTimer("buildNbr");
    }
}

/**
 * @brief Compute forces using neighbor list.
 *        For each i, for each (j, i_sym) in neighborList[i], compute force if r < R_CUT.
 * 
 *   loop for i in [0..N-1]:
 *     for (j, iSym) in vecNbrList[i]:
 *       => compute dx/dy/dz => sqrt => if (r<R_CUT) => LJ force => addForce
 * 
 * @param pBench Optional pointer for instrumentation (time, ops).
 */
void System::computeForcesWithVerlet(Benchmark* pBench /*=nullptr*/)
{
    // measure time
    if (pBench) {
        pBench->startTimer("forceCalc");
    }

    // reset forces
#pragma omp parallel for
    for (u64 idxA = 0; idxA < vecParts.size(); ++idxA) {
        vecParts[idxA].resetForce();
    }

    // approximate flops
    u64 locOps = 0ULL; 
    // each neighbor => distance => sqrt => ~20 ops => plus potential overhead

#pragma omp parallel for schedule(dynamic) reduction(+:locOps)
    for (u64 idxI = 0; idxI < vecParts.size(); ++idxI) {
        const auto& posI = vecParts[idxI].getPosition();

        // neighbor list => pairs
        for (auto& nPair : vecNbrList[idxI]) {
            u32 idxJ   = nPair.first;
            u32 idxSym = nPair.second;

            const auto& posJ = vecParts[idxJ].getPosition();

            // 1) compute local coords
            f64 xLoc = posJ[0] + vecSym[idxSym][0];
            f64 yLoc = posJ[1] + vecSym[idxSym][1];
            f64 zLoc = posJ[2] + vecSym[idxSym][2];

            // 2) compute distance
            f64 dX = posI[0] - xLoc;
            f64 dY = posI[1] - yLoc;
            f64 dZ = posI[2] - zLoc;
            f64 dist = std::sqrt(dX*dX + dY*dY + dZ*dZ);
            locOps += 20ULL; // rough

            // 3) check cutoff => if valid => compute LJ force
            if (dist < R_CUT) {
                // local part => overhead but design
                Particle tmpPrt({xLoc, yLoc, zLoc}, {0, 0, 0}, 1.0);
                tmpPrt.setBoxSize(boxDim);

                auto valForce = m_potential->computeForce(vecParts[idxI], tmpPrt);
                locOps += 20ULL; // rough

                // critical => block concurrency => big overhead
#pragma omp critical
                {
                    vecParts[idxI].addForce(valForce);

                    std::vector<f64> negF(3);
                    negF[0] = -valForce[0];
                    negF[1] = -valForce[1];
                    negF[2] = -valForce[2];
                    vecParts[idxJ].addForce(negF);
                }
            }
        }
    }

    // end timer + ops
    if (pBench) {
        pBench->addOps("forceCalc", locOps);
        pBench->endTimer("forceCalc");
    }
}

/**
 * @brief Compute total energy = sum(kinetic) + sum(potential).
 *        Potential part is O(N^2 * N_SYM).
 *        I do a triple loop: i_sym, i, j => if r^2 < R_CUT^2 => LJ pot.
 * 
 *   E_kin = 0.5*m*(vx^2 + vy^2 + vz^2).
 *   E_pot = sum_{i_sym=1..N_SYM} sum_{i<j} LJ(...).
 * 
 * @param pBench Optional pointer for instrumentation
 * @return f64 The total energy.
 */
f64 System::computeTotalEnergy(Benchmark* pBench /*=nullptr*/) const
{
    // measure time
    if (pBench) {
        pBench->startTimer("energyCalc");
    }

    f64 accumE = 0.0;

    // (A) Kinetic => short loop => O(N)
    for (const auto& onePart : vecParts) {
        const auto& velRef = onePart.getVelocity();
        // KE = 0.5*m*(v^2)
        f64 kE = 0.5 * onePart.getMass() * 
                 (velRef[0]*velRef[0] + velRef[1]*velRef[1] + velRef[2]*velRef[2]);
        accumE += kE;
    }

    // (B) Potential => triple nested => O(N^2*N_SYM)
    // I'll approximate ops
    u64 locOps = 0ULL;

    // loop for i_sym in [0..N_SYM-1], i in [0..N-1], j in [i+1..N-1]
    for (i32 symI = 0; symI < N_SYM; symI++) {
        for (u64 idxA = 0; idxA < vecParts.size(); idxA++) {
            for (u64 idxB = idxA + 1; idxB < vecParts.size(); idxB++) {
                const auto& posA = vecParts[idxA].getPosition();
                const auto& posB = vecParts[idxB].getPosition();

                f64 xLocB = posB[0] + vecSym[symI][0];
                f64 yLocB = posB[1] + vecSym[symI][1];
                f64 zLocB = posB[2] + vecSym[symI][2];

                f64 ddX = posA[0] - xLocB;
                f64 ddY = posA[1] - yLocB;
                f64 ddZ = posA[2] - zLocB;

                f64 dist2 = ddX*ddX + ddY*ddY + ddZ*ddZ;
                // ~8 ops
                locOps += 8ULL;

                if (dist2 < (R_CUT*R_CUT)) {
                    // local part => overhead but design
                    Particle tmpPrt({xLocB, yLocB, zLocB}, {0,0,0}, 1.0);
                    tmpPrt.setBoxSize(boxDim);

                    f64 potVal = m_potential->computePotentialEnergy(vecParts[idxA], tmpPrt);
                    // ~15 ops
                    locOps += 15ULL;

                    accumE += potVal;
                }
            }
        }
    }

    // finalize
    if (pBench) {
        pBench->addOps("energyCalc", locOps);
        pBench->endTimer("energyCalc");
    }

    return accumE;
}
