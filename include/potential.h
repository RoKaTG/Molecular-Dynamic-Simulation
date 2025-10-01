/**
 * @file potential.h
 * @brief Base Potential class, with pure virtual computeForce/computePotentialEnergy.
 */

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include "particle.h"
#include "types.h"

/**
 * @class Potential
 * @brief Abstract base for all potentials. 
 *        Must implement force & potential energy for two Particles.
 */
class Potential {
public:
    virtual ~Potential() = default;

    /**
     * @brief Compute force vector for p1 & p2 => 3D
     */
    virtual std::vector<f64> computeForce(const Particle& p1, const Particle& p2) const = 0;

    /**
     * @brief Compute potential energy for p1 & p2 => scalar
     */
    virtual f64 computePotentialEnergy(const Particle& p1, const Particle& p2) const = 0;
};

#endif // POTENTIAL_H
