#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include "types.h"

// Pour stocker les métriques d'un pas de simulation
struct StepMetrics {
    u32 step;
    f64 E_pot;
    f64 E_cin;
    f64 E_tot;
    f64 temperature;
    f64 sumFx;
    f64 sumFy;
    f64 sumFz;
    f64 gflops;  // gf lops mesuré ou estimé
};

class Benchmark {
public:
    Benchmark();
    ~Benchmark() = default;

    // 1) Timers
    void startTimer(const std::string& label);
    void endTimer(const std::string& label);

    // 2) Comptage d'opérations "manuel"
    void addOps(const std::string& label, u64 count);
    f64  getTotalOps(const std::string& label) const;

    // 3) Méthode autoCountFlopsInSource => plus fiable (ignore commentaires, strings)
    u64  autoCountFlopsInSource(const std::string& code);

    // 4) Affichage config
    void printSystemConfiguration(u32 nParticles, f64 rCut, f64 boxSize, u32 nSym);

    // 5) Enregistrement / affichage périodique
    void recordStepMetrics(u32 step, f64 E_pot, f64 E_cin, f64 E_tot,
                           f64 temperature, f64 sumFx, f64 sumFy, f64 sumFz,
                           f64 gflops);
    void printPeriodicMetrics(u32 step, u32 freq);

    // 6) Résumé final
    void printFinalSummary();

private:
    // removeComments => supprime // et /* */
    std::string removeComments(const std::string& code);

private:
    // Chrono
    std::unordered_map<std::string, std::chrono::time_point<std::chrono::high_resolution_clock>> m_startPoints;
    std::unordered_map<std::string, std::vector<f64>> m_durationsMs;

    // Comptage ops
    std::unordered_map<std::string, u64> m_opsCount;

    // Liste des StepMetrics
    std::vector<StepMetrics> m_stepMetrics;
};

#endif // BENCHMARK_H
