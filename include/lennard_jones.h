/**
 * @file lennard_jones.h
 * @brief Declaration of LennardJones potential class (inherited from Potential).
 *        Uses sigma, epsilon, optional rCut for LJ interactions.
 */

#ifndef LENNARD_JONES_H
#define LENNARD_JONES_H

#include "potential.h"

/**
 * @class LennardJones
 * @brief Implements a standard 12-6 Lennard-Jones potential:
 *        U(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6).
 *        Force is derivative of that expression.
 */
class LennardJones : public Potential {
public:
    /**
     * @brief Construct with given epsilon, sigma. 
     *        rCut set to 2.5*sigma by default.
     * @param epsilon E.g. 0.2
     * @param sigma   E.g. 3.0
     */
    LennardJones(f64 epsilon, f64 sigma);

    /**
     * @brief Compute LJ force between p1 & p2. 
     *        rCut check is done outside, so no cutoff check here.
     * @param p1 The first Particle
     * @param p2 The second Particle
     * @return 3D force vector
     */
    std::vector<f64> computeForce(const Particle& p1, const Particle& p2) const override;

    /**
     * @brief Compute potential energy for p1 & p2 under LJ. 
     *        rCut check is done outside, so no check here.
     * @param p1 The first Particle
     * @param p2 The second Particle
     * @return The LJ potential energy
     */
    f64 computePotentialEnergy(const Particle& p1, const Particle& p2) const override;

private:
    f64 epsVal;   ///< LJ epsilon
    f64 sigVal;   ///< LJ sigma
    f64 rCutVal;  ///< radius cutoff => 2.5*sigma
};

#endif // LENNARD_JONES_H
