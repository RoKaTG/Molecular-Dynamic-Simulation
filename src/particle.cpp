/**
 * @file particle.cpp
 * @brief Implementation of Particle class. 
 */

#include "particle.h"
#include "types.h"

/**
 * @brief Default ctor. pos=vel=frc=acc=0, mass=1.0, boxDim=0.
 */
Particle::Particle()
  : posVec(3, 0.0),
    velVec(3, 0.0),
    frcVec(3, 0.0),
    accVec(3, 0.0),
    massVal(1.0),
    boxDim(0.0)
{}

/**
 * @brief Construct with pos, vel, mass. frc/acc=0, boxDim=0
 */
Particle::Particle(const std::vector<f64>& position,
                   const std::vector<f64>& velocity,
                   f64 mass)
  : posVec(position),
    velVec(velocity),
    frcVec(3, 0.0),
    accVec(3, 0.0),
    massVal(mass),
    boxDim(0.0)
{}

/**
 * @brief Return const ref to posVec
 */
const std::vector<f64>& Particle::getPosition() const {
    return posVec;
}

/**
 * @brief Return const ref to velVec
 */
const std::vector<f64>& Particle::getVelocity() const {
    return velVec;
}

/**
 * @brief Return const ref to frcVec
 */
const std::vector<f64>& Particle::getForce() const {
    return frcVec;
}

/**
 * @brief Return const ref to accVec
 */
const std::vector<f64>& Particle::getAcceleration() const {
    return accVec;
}

/**
 * @brief Return massVal
 */
f64 Particle::getMass() const {
    return massVal;
}

/**
 * @brief Overwrite posVec
 */
void Particle::setPosition(const std::vector<f64>& position) {
    posVec = position;
}

/**
 * @brief Overwrite velVec
 */
void Particle::setVelocity(const std::vector<f64>& velocity) {
    velVec = velocity;
}

/**
 * @brief Overwrite frcVec
 */
void Particle::setForce(const std::vector<f64>& force) {
    frcVec = force;
}

/**
 * @brief Overwrite accVec
 */
void Particle::setAcceleration(const std::vector<f64>& accel) {
    accVec = accel;
}

/**
 * @brief Reset frcVec to zero
 */
void Particle::resetForce() {
    frcVec = std::vector<f64>(3, 0.0);
}

/**
 * @brief Add force to existing frcVec => frcVec += param
 */
void Particle::addForce(const std::vector<f64>& force) {
    // short loop => 3 comps
    for (u32 idx = 0; idx < frcVec.size(); ++idx) {
        frcVec[idx] += force[idx];
    }
}

/**
 * @brief Set boxDim
 */
void Particle::setBoxSize(f64 box_dim) {
    boxDim = box_dim;
}

/**
 * @brief Return boxDim
 */
f64 Particle::getBoxSize() const {
    return boxDim;
}
