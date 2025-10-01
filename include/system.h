/**
 * @file system.h
 * @brief Header for System class: particle container, neighbor list building, force & energy computations.
 * 
 * The System class holds:
 * - A vector of Particles (vecParts)
 * - A pointer to a Potential (pPot)
 * - A box dimension (boxDim)
 * - 27 offset vectors for periodic images (vecSym)
 * - A neighbor list (vecNbrList) used by Verlet approach
 * - Add particles
 * - Set potential
 * - Build neighbor list (O(N^2 * N_SYM))
 * - Compute forces using neighbor list
 * - Compute total energy (kinetic + potential)
 * - Manage box dimension
 * - Initialize symmetry vectors for periodic images
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <memory>
#include <array>

#include "particle.h"
#include "potential.h"
#include "types.h"
#include "config.h"

// Forward declaration if needed
class Benchmark;

/**
 * @class System
 * @brief Container for Particles, building neighbor list, computing forces & energies, etc.
 */
class System {
public:
    /**
     * @brief Construct a System with empty data. Potential= nullptr, boxDim=0.0.
     */
    System();

    /**
     * @brief Add a new Particle to internal vector. 
     *        Sets the box dimension for that particle too.
     * @param inPart The incoming Particle
     */
    void addParticle(const Particle& inPart);

    /**
     * @brief Set the Potential pointer (e.g. LennardJones).
     * @param pot Shared pointer to Potential object.
     */
    void setPotential(std::shared_ptr<Potential> pot);

    /**
     * @brief Return reference to vector of Particles (modifiable).
     */
    std::vector<Particle>& getParticles();

    /**
     * @brief Return const reference to vector of Particles.
     */
    const std::vector<Particle>& getParticles() const;

    /**
     * @brief Old function O(N^2 * 27). Left for reference (commented out in .cpp).
     *        Possibly not used if Verlet approach is chosen.
     * @note If re-enabled, we can add Benchmark* param as well.
     */
    void computeForces(); 

    /**
     * @brief Compute forces using neighbor list. 
     *        For each i, for each (j, i_sym) in vecNbrList[i], compute LJ if dist < R_CUT.
     * @param pBench Optional pointer for instrumentation (timing, ops).
     */
    void computeForcesWithVerlet(Benchmark* pBench=nullptr);

    /**
     * @brief Build neighbor list in O(N^2 * N_SYM).
     *        For each (i_sym, i<j), store (j, i_sym) in vecNbrList[i] if r^2 < R_CUT^2.
     * @param pBench Optional pointer for instrumentation (timing, ops).
     */
    void buildNeighborList(Benchmark* pBench=nullptr);

    /**
     * @brief Compute total energy = sum(kinetic) + sum(potential).
     *        Potential part done via triple loop i_sym, i, j => O(N^2*N_SYM).
     * @param pBench Optional pointer for instrumentation (timing, ops).
     * @return f64 The total energy (kin + pot).
     */
    f64 computeTotalEnergy(Benchmark* pBench=nullptr) const;

    /**
     * @brief Set the cubic box dimension. Also updates each Particle's boxDim.
     * @param newBox The new dimension.
     */
    void setBoxSize(f64 newBox);

    /**
     * @brief Get current box dimension.
     * @return f64
     */
    f64 getBoxSize() const;

    /**
     * @brief Initialize 27 offset vectors for periodic images (vecSym).
     */
    void initSymVectors();

    /**
     * @brief Return const reference to offset vectors (vecSym).
     */
    const std::vector<std::array<f64,3>>& getSymVectors() const;

private:
    // *** Data Members ***

    std::vector<Particle> vecParts; ///< Internal vector of Particles
    std::shared_ptr<Potential> m_potential; ///< Potential pointer
    f64 boxDim; ///< The box dimension

    // 27 offset vectors for periodic images
    std::vector<std::array<f64,3>> vecSym;

    // For each i => vector of (j, i_sym) pairs => neighbor list
    std::vector<std::vector<std::pair<u32,u32>>> vecNbrList;
};

#endif // SYSTEM_H
