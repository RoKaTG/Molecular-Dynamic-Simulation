/**
 * @file md_utils.cpp
 * @brief Implementation of MD utilities: kineticEnergy, temperature, velocity init, PDB output.
 *
 */

#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <cstring>

#include "md_utils.h"
#include "config.h"  
#include "benchmark.h" 

/**
 * @brief E_kin = 1/(2*CONVERSION_FORCE)* sum(m_i*v^2). 
 *        We do a short loop over Particles.
 */
f64 computeKineticEnergy(const System& sys, Benchmark* pBench /*=nullptr*/) {
    if (pBench) {
        pBench->startTimer("kinEnergy");
    }

    // ~5 ops per particle => (v^2 => 3 mul + sum + mass?)
    u64 locOps = 0ULL;

    const auto& vecP = sys.getParticles();
    f64 accum = 0.0;
    for (const auto& oneP : vecP) {
        const auto& velV = oneP.getVelocity();
        f64 vsq = velV[0]*velV[0] + velV[1]*velV[1] + velV[2]*velV[2];
        accum += (oneP.getMass() * vsq);
        locOps += 5ULL; // approx
    }

    // multiply by 1/(2*CONVERSION_FORCE)
    accum *= (1.0 / (2.0 * CONVERSION_FORCE));

    if (pBench) {
        pBench->addOps("kinEnergy", locOps);
        pBench->endTimer("kinEnergy");
    }
    return accum;
}

/**
 * @brief T= (1/(Ndl*CONSTANTE_R))* Ecin. Ndl=3*N -3
 */
f64 computeTemperature(f64 eCin, u32 N) {
    u32 Ndl = 3*N - 3;
    f64 T = (1.0 / (Ndl * CONSTANTE_R)) * eCin;
    return T;
}

/**
 * @brief 1) random in [-0.5,0.5], 2) subtract center-of-mass, 
 *        3) compute Ecin => rescale => T0
 */
void initializeVelocities(System& sys, f64 T0, Benchmark* pBench /*=nullptr*/) {
    if (pBench) {
        pBench->startTimer("initVel");
    }

    auto& vecP = sys.getParticles();
    u32 N = (u32)vecP.size();
    if (N==0) {
        if (pBench) {
            pBench->endTimer("initVel");
        }
        return;
    }

    // ~some ops
    u64 locOps = 0ULL;

    // 1) random in [-0.5,0.5]
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<f64> distR(-0.5, 0.5);

    for (auto &pp : vecP) {
        std::vector<f64> vtmp(3);
        for (u32 i=0; i<3; i++){
            vtmp[i] = distR(gen);
        }
        pp.setVelocity(vtmp);
        locOps += 3ULL; // 3 random
    }

    // 2) subtract center-of-mass velocity
    std::vector<f64> sumV(3, 0.0);
    for (auto &pp : vecP) {
        const auto& vRef = pp.getVelocity();
        sumV[0] += vRef[0];
        sumV[1] += vRef[1];
        sumV[2] += vRef[2];
        locOps += 3ULL; 
    }
    for (u32 i = 0; i < 3; i++) {
        sumV[i] /= (f64)N;
    }
    locOps += 3ULL; // dividing

    for (auto &pp : vecP) {
        auto vRef = pp.getVelocity();
        vRef[0] -= sumV[0];
        vRef[1] -= sumV[1];
        vRef[2] -= sumV[2];
        pp.setVelocity(vRef);
        locOps += 3ULL; 
    }

    // 3) compute Ecin => rescale => T0
    f64 eCinInit = computeKineticEnergy(sys, pBench);
    // do not double count ops from computeKineticEnergy

    u32 Ndl = 3*N - 3;
    f64 ratio = (Ndl * CONSTANTE_R * T0) / eCinInit;
    f64 scale = std::sqrt(ratio);

    for (auto &pp : vecP) {
        auto vRef = pp.getVelocity();
        vRef[0] *= scale;
        vRef[1] *= scale;
        vRef[2] *= scale;
        pp.setVelocity(vRef);
        locOps += 3ULL;
    }

    if (pBench) {
        pBench->addOps("initVel", locOps);
        pBench->endTimer("initVel");
    }
}

/**
 * @brief Save system state to PDB file. 
 *        We do naive alignment, with 80-col lines.
 */
void savePDBFrameExact(const System& sys, 
                       u32 iteration, 
                       const std::string& outFile,
                       Benchmark* pBench /*=nullptr*/)
{
    if (pBench) {
        pBench->startTimer("savePDB");
    }

    std::ofstream ofs(outFile, std::ios::app);
    if (!ofs.is_open()) {
        std::cerr << "Err: cannot open " << outFile << " in append.\n";
        if (pBench) pBench->endTimer("savePDB");
        return;
    }

    // We skip counting ops => mostly I/O
    // Just do time measure

    // 1) CRYST1 ...
    {
        f64 boxL = sys.getBoxSize();
        int Lint = (int)std::round(boxL);
        char cbuf[128];
        std::snprintf(cbuf, sizeof(cbuf),
            "CRYST1   %2d   %2d   %2d  90.00  90.00  90.00 P",
            Lint, Lint, Lint
        );
        ofs << cbuf << "\n";
    }

    // 2) MODEL
    {
        char cbufM[128];
        std::snprintf(cbufM, sizeof(cbufM), "MODEL        %u", iteration);
        ofs << cbufM << "\n";
    }

    // 3) lines ATOM ...
    const auto& vecP = sys.getParticles();
    for (u32 i=0; i<vecP.size(); i++){
        f64 x = vecP[i].getPosition()[0];
        f64 y = vecP[i].getPosition()[1];
        f64 z = vecP[i].getPosition()[2];

        char lineA[128];
        std::memset(lineA, ' ', 80);
        lineA[80] = '\0';

        // col1..4 => "ATOM"
        std::memcpy(&lineA[0], "ATOM", 4);

        // col10 => i => can be 1 or 2 digits
        {
            char tmp[8];
            std::snprintf(tmp, sizeof(tmp), "%u", i);
            lineA[9] = tmp[0];
            if (tmp[1] != '\0') {
                lineA[10] = tmp[1];
            }
        }

        lineA[12] = 'C';   // col13 => 'C'
        lineA[22] = '0';   // col23 => '0'

        // X => col31..37 => %7.3f
        {
            char bx[16];
            std::snprintf(bx, sizeof(bx), "%7.3f", x);
            std::memcpy(&lineA[30], bx, std::strlen(bx));
        }

        lineA[37] = ' '; // space
        // Y => col39..45
        {
            char by[16];
            std::snprintf(by, sizeof(by), "%7.3f", y);
            std::memcpy(&lineA[38], by, std::strlen(by));
        }

        lineA[45] = ' ';
        // Z => col47..53
        {
            char bz[16];
            std::snprintf(bz, sizeof(bz), "%7.3f", z);
            std::memcpy(&lineA[46], bz, std::strlen(bz));
        }

        // col70..73 => "MRES"
        std::memcpy(&lineA[69], "MRES", 4);

        // col73 => '\0'
        lineA[73] = '\0';

        ofs << lineA << "\n";
    }

    // 4) TER
    ofs << "TER\n";

    // 5) ENDMDL
    ofs << "ENDMDL\n";

    ofs.close();

    if (pBench) {
        pBench->endTimer("savePDB");
    }
}
