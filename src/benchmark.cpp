/**
 * @file benchmark.cpp
 * @brief Implementation of Benchmark class with ASCII table summary.
 */

#include "benchmark.h"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>

// constructor
Benchmark::Benchmark()
{
}

// startTimer
void Benchmark::startTimer(const std::string& label)
{
    m_startPoints[label] = std::chrono::high_resolution_clock::now();
}

// endTimer
void Benchmark::endTimer(const std::string& label)
{
    auto endT = std::chrono::high_resolution_clock::now();
    auto it = m_startPoints.find(label);
    if (it==m_startPoints.end()) {
        return;
    }
    f64 durMs = std::chrono::duration<f64,std::milli>(endT - it->second).count();
    m_durationsMs[label].push_back(durMs);
}

// addOps
void Benchmark::addOps(const std::string& label, u64 count)
{
    m_opsCount[label]+= count;
}

// getTotalOps
f64 Benchmark::getTotalOps(const std::string& label) const
{
    auto it= m_opsCount.find(label);
    if (it==m_opsCount.end()) {
        return 0.0;
    }
    return (f64)it->second;
}

// recordStepMetrics
void Benchmark::recordStepMetrics(u32 step, f64 E_pot, f64 E_cin, f64 E_tot,
                                  f64 temperature, f64 sumFx, f64 sumFy, f64 sumFz,
                                  f64 gflops)
{
    StepMetrics sm;
    sm.step= step;
    sm.E_pot= E_pot;
    sm.E_cin= E_cin;
    sm.E_tot= E_tot;
    sm.temperature= temperature;
    sm.sumFx= sumFx;
    sm.sumFy= sumFy;
    sm.sumFz= sumFz;
    sm.gflops= gflops;
    m_stepMetrics.push_back(sm);
}

// printPeriodicMetrics
void Benchmark::printPeriodicMetrics(u32 step, u32 freq)
{
    if (step%freq!=0) return;
    if (m_stepMetrics.empty()) return;

    const auto& last= m_stepMetrics.back();
    std::cout << "[Metrics] Step " << last.step
              << " E_pot=" << last.E_pot
              << " E_cin=" << last.E_cin
              << " E_tot=" << last.E_tot
              << " T=" << last.temperature
              << " sumF=("<< last.sumFx<<","<< last.sumFy<<","<< last.sumFz<<")"
              << " GFLOPS="<< last.gflops
              << std::endl;
}

/**
 * @brief Print a final summary in ASCII table format:
 *        - Timers => label, count, avg, stdev
 *        - Ops => label, total ops, GFLOPS
 *        - Global stats => totalExec time => totalOps => total GFLOPS
 */
void Benchmark::printFinalSummary()
{
    std::cout << "=== Final Benchmark Summary ===\n\n";

    // Timers table
    // columns: label(18), count(6), avg(10), stdev(10)
    // plus separators => 
    // e.g. => +---------------------------------------------------------------+
    std::cout << "+---------------------------------------------------------------+\n";
    std::cout << "| Timers (ms)       | count |   avg(ms) |  stdev(ms) |\n";
    std::cout << "+---------------------------------------------------------------+\n";

    f64 totalExecMs=0.0;
    bool haveTotalExec= false;

    // gather timer labels
    std::vector<std::string> tLabels;
    for (auto& kv : m_durationsMs){
        tLabels.push_back(kv.first);
    }
    std::sort(tLabels.begin(), tLabels.end());

    for (auto& label : tLabels) {
        const auto& arr= m_durationsMs[label];
        if (arr.empty()) continue;

        f64 sum=0.0;
        for (auto val: arr) sum+= val;
        f64 mean= sum/ arr.size();

        f64 accum=0.0;
        for (auto val: arr){
            f64 d= (val-mean);
            accum+= d*d;
        }
        f64 stdev= (arr.size()>1)? std::sqrt(accum/ arr.size()): 0.0;

        if (label=="totalExec"){
            totalExecMs= sum;
            haveTotalExec= true;
        }

        // label => 16 chars
        // count => 5
        // avg => 9
        // stdev => 9
        std::cout << "| " 
                  << std::left << std::setw(16) << label << " | "
                  << std::right<< std::setw(5) << arr.size() << " | "
                  << std::setw(8) << std::fixed << std::setprecision(5) << mean << " | "
                  << std::setw(9) << std::fixed << std::setprecision(5) << stdev 
                  << " |\n";

        std::cout << "+---------------------------------------------------------------+\n";
    }

    // Ops table
    std::cout << "\n+---------------------------------------------------------------+\n";
    std::cout << "| Ops count & GFLOPS                                            |\n";
    std::cout << "+---------------------------------------------------------------+\n";
    std::cout << "| label            |       ops   |   GFLOPS   |\n";
    std::cout << "+---------------------------------------------------------------+\n";

    std::vector<std::string> opLabels;
    for (auto& kv : m_opsCount){
        opLabels.push_back(kv.first);
    }
    std::sort(opLabels.begin(), opLabels.end());

    u64 totalOps=0ULL;

    for (auto& lb : opLabels){
        u64 ops= m_opsCount[lb];

        // sum durations => to get GF
        f64 sumMs=0.0;
        auto it= m_durationsMs.find(lb);
        if (it!=m_durationsMs.end()){
            for (auto v: it->second) sumMs+= v;
        }
        f64 gf=0.0;
        f64 sumSec= sumMs/1000.0;
        if (sumSec>0.0){
            gf= ((f64)ops /1e9)/ sumSec;
        }

        totalOps+= ops;

        // row => label(14), ops(10), gf(8)
        std::cout << "| " 
                  << std::left << std::setw(14) << lb << " | "
                  << std::right<< std::setw(10)<< ops << " | "
                  << std::setw(9)<< std::fixed << std::setprecision(6)<< gf
                  <<" |\n";
        std::cout << "+---------------------------------------------------------------+\n";
    }

    // final global stats
    std::cout << "\n+---------------------------------------------------------------+\n";
    std::cout << "|                   GLOBAL STATS                                |\n";
    std::cout << "+---------------------------------------------------------------+\n";

    if (haveTotalExec){
        f64 totalSec= totalExecMs/1000.0;
        f64 gfGlobal=0.0;
        if (totalSec>0.0){
            gfGlobal= ((f64)totalOps /1e9)/ totalSec;
        }

        std::cout<<"| totalExec time (ms)= "<< totalExecMs 
                 <<" => "<< totalSec <<" s\n"
                 <<"| totalOps= "<< totalOps <<"\n"
                 <<"| => total GFLOPS= "<< gfGlobal <<"\n";
    } else {
        std::cout<<"| No 'totalExec' label => can't compute global time/GFLOPS.\n";
    }
    std::cout<<"+---------------------------------------------------------------+\n\n";
}
