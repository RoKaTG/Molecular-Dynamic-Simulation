/**
 * @file main.cpp
 * @brief Main entry for the MD simulation. 
 *        Prints config, runs velocity-Verlet, logs steps in a single table, 
 *        prints final benchmark summary.
 * 
 */

#include <cstdio>      // for std::remove
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>     // for file output

#include "types.h"
#include "config.h"

#include "particle.h"
#include "system.h"
#include "lennard_jones.h"
#include "verlet_integrator.h"
#include "md_utils.h"
#include "utils.h"
#include "benchmark.h"

int main(int /*argc*/, char** /*argv*/)
{
    // Create a Benchmark object
    Benchmark benchObj;
    // Start measuring total execution
    benchObj.startTimer("totalExec");

    // 1) Remove old PDB
    std::string outPdb = "out.pdb";
    std::remove(outPdb.c_str());

    // 2) Prepare system
    System sysObj;
    f64 boxDim = 42.0;
    sysObj.setBoxSize(boxDim);

    // Read file => store Particles
    std::vector<Particle> inpParts;
    readParticlesFromFile("../particule.xyz", inpParts, boxDim, &benchObj);

    f64 Tinit = 300.0;

    // Print config as ASCII table
    std::cout << "\n+----------------------------------------+\n";
    std::cout << "|            CONFIGURATION               |\n";
    std::cout << "+----------------------------------------+\n";
    std::cout << "| Number of particles read : " << inpParts.size() << "\n";
    std::cout << "| N_PARTICULES_TOTAL      : " << N_PARTICULES_TOTAL << "\n";
    std::cout << "| N_PARTICULES_LOCAL      : " << N_PARTICULES_LOCAL << "\n";
    std::cout << "| BoxDim                  : " << boxDim << "\n";
    std::cout << "| R_CUT                   : " << R_CUT << "\n";
    std::cout << "| N_SYM                   : " << N_SYM << "\n";
    std::cout << "| T0 (K)                  : " << Tinit << "\n";
    std::cout << "+----------------------------------------+\n\n";

    // Add only N_PARTICULES_LOCAL
    for (u32 idx = 0; idx < N_PARTICULES_LOCAL && idx < (u32)inpParts.size(); ++idx) {
        sysObj.addParticle(inpParts[idx]);
    }

    // init 27 images
    sysObj.initSymVectors();

    // create LJ potential
    auto potPtr = std::make_shared<LennardJones>(0.2, 3.0);
    sysObj.setPotential(potPtr);

    // Build neighbor => O(N^2*N_SYM)
    sysObj.buildNeighborList(&benchObj);

    // Compute initial forces
    sysObj.computeForcesWithVerlet(&benchObj);

    // Store a(t)=F/m
    {
        auto& pVec = sysObj.getParticles();
        for (u32 i = 0; i < pVec.size(); i++) {
            const auto& frcRef = pVec[i].getForce();
            f64 mVal = pVec[i].getMass();
            std::vector<f64> a0(3);
            for (u32 j = 0; j < 3; j++){
                a0[j] = frcRef[j]/mVal;
            }
            pVec[i].setAcceleration(a0);
        }
    }

    // 5) init velocities => T0
    initializeVelocities(sysObj, Tinit, &benchObj);

    // 6) MD loop => velocity-Verlet
    VerletIntegrator integrObj; 
    u32 stepCount = 500;
    f64 dtVal     = 1.0; // 1 fs

    // Print header for the step logs table
    std::cout << "+-----------------------------------------------------------------------------------------------------+\n";
    std::cout << "| Step  |    E_pot    |  E_cin  |    E_tot    |    T    |         sumF(x, y, z)                |\n";
    std::cout << "+-----------------------------------------------------------------------------------------------------+\n";


    {
        std::ofstream ftemp("../data/temp.dat");
        ftemp << "# step  T\n";
    }
    {
        std::ofstream fkin("../data/kinetic_energy.dat");
        fkin << "# step  E_kin\n";
    }
    {
        std::ofstream fpot("../data/potential_energy.dat");
        fpot << "# step  E_pot\n";
    }
    {
        std::ofstream ftot("../data/total_energy.dat");
        ftot << "# step  E_tot\n";
    }

    for (u32 step = 0; step < stepCount; step++){
        // integrate
        integrObj.integrate(sysObj, dtVal, &benchObj);

        // rebuild neighbor every 20 steps
        if (step % 20 == 0 && step > 0) {
            sysObj.buildNeighborList(&benchObj);
        }

        // recalc forces
        sysObj.computeForcesWithVerlet(&benchObj);

        // Ecin
        f64 EcinVal = computeKineticEnergy(sysObj, &benchObj);

        // E_tot
        f64 Eall = sysObj.computeTotalEnergy(&benchObj);
        f64 Epot = Eall - EcinVal;

        // T
        f64 Tval = computeTemperature(EcinVal, (u32)sysObj.getParticles().size());

        // => Append each metric to its own file
        {
            std::ofstream ftemp("../data/temp.dat", std::ios::app);
            ftemp << step << " " << Tval << "\n";
        }
        {
            std::ofstream fkin("../data/kinetic_energy.dat", std::ios::app);
            fkin << step << " " << EcinVal << "\n";
        }
        {
            std::ofstream fpot("../data/potential_energy.dat", std::ios::app);
            fpot << step << " " << Epot << "\n";
        }
        {
            std::ofstream ftot("../data/total_energy.dat", std::ios::app);
            ftot << step << " " << Eall << "\n";
        }

        // sumF => debug
        std::vector<f64> sumForce(3, 0.0);
        {
            const auto& pVec= sysObj.getParticles();
            for (u32 i=0; i<pVec.size(); i++){
                const auto& frcRef= pVec[i].getForce();
                sumForce[0]+= frcRef[0];
                sumForce[1]+= frcRef[1];
                sumForce[2]+= frcRef[2];
            }
        }

        // save PDB every 5 steps
        if (step % 5 == 0) {
            savePDBFrameExact(sysObj, step, outPdb, &benchObj);
        }

        // print step info every 50 steps
        if (step % 50 == 0) {
            std::ostringstream rowStr;
            rowStr << "| " << std::setw(4) << step << " | "
                   << std::setw(10) << std::fixed << std::setprecision(2) << Epot << " | "
                   << std::setw(10) << std::fixed << std::setprecision(2) << EcinVal << " | "
                   << std::setw(10) << std::fixed << std::setprecision(2) << Eall << " | "
                   << std::setw(6) << std::fixed << std::setprecision(2) << Tval << " | ";

            // sumF => 3 scientific
            std::ostringstream sumFcol;
            sumFcol << std::scientific << std::setprecision(3)
                    << sumForce[0] << ", "
                    << sumForce[1] << ", "
                    << sumForce[2];
            std::string sumFStr = sumFcol.str();

            if (sumFStr.size() < 48) {
                sumFStr.append(48 - sumFStr.size(),' ');
            } else if (sumFStr.size() > 48) {
                sumFStr = sumFStr.substr(0,48);
            }

            rowStr << sumFStr << "          |\n";

            std::cout << rowStr.str();

            // line separator
            std::cout << "+-----------------------------------------------------------------------------------------------------+\n";
        }
    }

    std::cout << "\n=== Simulation finished ===\n";

    // stop totalExec
    benchObj.endTimer("totalExec");

    std::cout << "\n--- Benchmark Summary ---\n\n";
    // print final summary in ASCII table
    benchObj.printFinalSummary();

    return 0;
}
