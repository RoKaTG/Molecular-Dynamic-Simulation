/**
 * @file config.h
 * @brief Global constants and macros for the simulation
 *
 * - N_PARTICULES_TOTAL / N_PARTICULES_LOCAL
 * - N_SYM = 27 (for periodic images)
 * - L = 42.0 (box dimension)
 * - R_CUT = 10.0
 * - MASS_PARTICLE = 18.0
 * - CONVERSION_FORCE = 0.0001 * 4.186
 * - CONSTANTE_R = 0.00199
 */

#ifndef CONFIG_H
#define CONFIG_H

#include "types.h"

// total vs local
static const i32 N_PARTICULES_TOTAL = 1000;
static const i32 N_PARTICULES_LOCAL = 1000; 

// number of symmetry images
static const i32 N_SYM = 27;

// box dimension, LJ cutoff
static const f64 L = 42.0; 
static const f64 R_CUT = 10.0;

// default mass, force conversion, gas constant
static const f64 MASS_PARTICLE     = 18.0;      
static const f64 CONVERSION_FORCE  = 0.0001 * 4.186; // ~0.0004186
static const f64 CONSTANTE_R       = 0.00199;     

#endif // CONFIG_H
