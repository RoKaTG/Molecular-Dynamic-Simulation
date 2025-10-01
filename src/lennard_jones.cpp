/**
 * @file lennard_jones.cpp
 * @brief Implementation of LennardJones class. 
 */

#include <cmath>
#include "lennard_jones.h"
#include "types.h"

/**
 * @brief Construct with given epsilon & sigma. rCutVal=2.5*sigma.
 */
LennardJones::LennardJones(f64 epsilon, f64 sigma)
    : epsVal(epsilon), sigVal(sigma), rCutVal(2.5 * sigma)
{}

/**
 * @brief Compute LJ force for p1 & p2, ignoring any cutoff check.
 * 
 * Force = 24*eps * ( (sig/r)^12 - (sig/r)^6 ) / r^2  * (dx,dy,dz)
 * 
 * Implementation detail: I do manual exponent:
 *   s = sigVal / r
 *   s2 = s*s
 *   s4 = s2*s2
 *   s6 = s2*s4
 *   s12 = s6*s6
 * => fewer std::pow calls
 */
std::vector<f64> LennardJones::computeForce(const Particle& p1, const Particle& p2) const {
    std::vector<f64> outF(3, 0.0);

    const auto& posA = p1.getPosition();
    const auto& posB = p2.getPosition();

    // dx, dy, dz
    f64 dX = posA[0] - posB[0];
    f64 dY = posA[1] - posB[1];
    f64 dZ = posA[2] - posB[2];

    f64 dist2 = dX*dX + dY*dY + dZ*dZ;
    f64 dist  = std::sqrt(dist2);

    // s = sigVal / dist
    f64 s = sigVal / dist;
    f64 s2 = s*s;
    f64 s4 = s2*s2;
    f64 s6 = s2*s4;
    f64 s12= s6*s6;

    // coefficient => 24 * epsVal * (s^12 - s^6) / r^2
    f64 factor = 24.0 * epsVal * (s12 - s6) / dist2;

    // multiply factor by (dx,dy,dz)
    outF[0] = factor * dX;
    outF[1] = factor * dY;
    outF[2] = factor * dZ;

    return outF;
}

/**
 * @brief Compute LJ potential energy for p1 & p2, ignoring cutoff check.
 * 
 * U = 4*epsVal * ( (sig/r)^12 - (sig/r)^6 )
 * 
 * Similarly, do manual exponent. 
 */
f64 LennardJones::computePotentialEnergy(const Particle& p1, const Particle& p2) const {
    const auto& posA = p1.getPosition();
    const auto& posB = p2.getPosition();

    f64 dX = posA[0] - posB[0];
    f64 dY = posA[1] - posB[1];
    f64 dZ = posA[2] - posB[2];

    f64 dist2 = dX*dX + dY*dY + dZ*dZ;
    f64 dist  = std::sqrt(dist2);

    f64 s = sigVal / dist;
    f64 s2 = s*s;
    f64 s4 = s2*s2;
    f64 s6 = s2*s4;
    f64 s12= s6*s6;

    f64 potVal = 4.0 * epsVal * (s12 - 2.0 * s6); 
    // note: the usual LJ form => 4*eps*( (sig/r)^12 - (sig/r)^6 )
    // but we do -2*s^6 => eq. to s^12 - 2*s^6
    // (s12 - 2*s6) is the standard => ( (sig/r)^12 - 2*(sig/r)^6 )

    return potVal;
}
