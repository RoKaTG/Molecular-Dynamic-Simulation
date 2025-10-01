/**
 * @file particle.h
 * @brief Declaration of Particle class, storing position, velocity, force, acceleration, etc.
 *        Box dimension is also stored for potential boundary logic.
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "types.h"

/**
 * @class Particle
 * @brief Represents a single particle with position, velocity, force, etc.
 */
class Particle {
public:
    /**
     * @brief Default ctor: pos=0, vel=0, frc=0, acc=0, mass=1.0, boxDim=0.0
     */
    Particle();

    /**
     * @brief Construct with given pos, vel, mass. Force & accel set to zero. boxDim=0.0
     * @param position 3D coords
     * @param velocity 3D velocity
     * @param mass The mass
     */
    Particle(const std::vector<f64>& position, 
             const std::vector<f64>& velocity, 
             f64 mass);

    /**
     * @brief Get position vector (3D).
     * @return const reference to internal posVec
     */
    const std::vector<f64>& getPosition() const;

    /**
     * @brief Get velocity vector (3D).
     * @return const reference
     */
    const std::vector<f64>& getVelocity() const;

    /**
     * @brief Get force vector (3D).
     * @return const reference
     */
    const std::vector<f64>& getForce() const;

    /**
     * @brief Get acceleration vector (3D).
     * @return const reference
     */
    const std::vector<f64>& getAcceleration() const; 

    /**
     * @brief Return mass
     */
    f64 getMass() const;

    /**
     * @brief Set position to new 3D vector
     */
    void setPosition(const std::vector<f64>& position);

    /**
     * @brief Set velocity to new 3D vector
     */
    void setVelocity(const std::vector<f64>& velocity);

    /**
     * @brief Set force to new 3D vector
     */
    void setForce(const std::vector<f64>& force);

    /**
     * @brief Set acceleration to new 3D vector
     */
    void setAcceleration(const std::vector<f64>& accel); 

    /**
     * @brief Reset force to (0,0,0)
     */
    void resetForce();

    /**
     * @brief Add a force to existing force vector => frcVec += param
     */
    void addForce(const std::vector<f64>& force);

    /**
     * @brief Set box dimension (for boundary usage).
     */
    void setBoxSize(f64 box_dim);

    /**
     * @brief Get box dimension
     */
    f64 getBoxSize() const;

private:
    // renamed to shorter or more distinct
    std::vector<f64> posVec; ///< 3D position
    std::vector<f64> velVec; ///< 3D velocity
    std::vector<f64> frcVec; ///< 3D force
    std::vector<f64> accVec; ///< 3D acceleration
    f64 massVal;             ///< mass
    f64 boxDim;              ///< dimension of bounding box
};

#endif // PARTICLE_H
