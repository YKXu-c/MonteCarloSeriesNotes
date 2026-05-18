/**
 * @file metropolis.cpp
 * @brief Metropolis single-spin-flip Monte Carlo for the 2D Ising model.
 *
 * One "sweep" = N random single-spin-flip attempts.  For each attempt:
 *   1. Pick a random site i
 *   2. Compute ΔE = E(-S_i) - E(S_i) = 2J Σ_{nn} S_i S_j + 2Jp Σ_{nnn} S_i S_k
 *   3. Accept with probability min(1, exp(-β ΔE))
 *
 * @ref  Metropolis et al., J. Chem. Phys. 21, 1087 (1953)
 * @complexity O(N) per sweep (each site visited once on average).
 *
 * Usage:
 *   metropolis --L 16 --J 1.0 --T 2.5 --sweeps 10000 --therm 1000 [--Jp 0.2] [--seed 42]
 */

#include "../include/ising_model.hpp"
#include "../include/mc_base.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <print>
#include <string>
#include <string_view>

using namespace mc;

// ---------------------------------------------------------------------------
// MetropolisSweep — CRTP update rule
// ---------------------------------------------------------------------------

/**
 * @brief Metropolis single-spin-flip sweep for IsingModel.
 *
 * @ref Metropolis et al., J. Chem. Phys. 21, 1087 (1953)
 * @complexity O(N) per sweep.
 */
class MetropolisSweep : public MCSamplerCRTP<MetropolisSweep> {
public:
    /**
     * @brief Construct a Metropolis sweep bound to an IsingModel.
     * @param model  Reference to the model (must outlive this object).
     */
    explicit MetropolisSweep(IsingModel& model) : model_(model) {}

    /**
     * @brief Perform one full sweep: N random single-spin-flip attempts.
     *
     * For each attempt:
     *   - Pick random site i
     *   - Compute ΔE (O(1) via neighbor lookup)
     *   - Accept with probability min(1, exp(-β ΔE))
     *
     * This satisfies detailed balance (see slide Frame 3.3):
     *   W(S→S')/W(S'→S) = min(1, e^{-βΔE}) / min(1, e^{βΔE}) = e^{-βΔE}
     */
    void sweep_impl() {
        std::uniform_int_distribution<int> site_dist(0, model_.N() - 1);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const double beta = model_.beta();

        for (int step = 0; step < model_.N(); ++step) {
            // Step 1: Pick a random site
            int site = site_dist(model_.rng());

            // Step 2: Compute energy change from flipping this spin
            //   ΔE = 2J Σ_{nn} S_i S_j  +  2Jp Σ_{nnn} S_i S_k
            double deltaE = model_.energyChange(site);

            // Step 3: Accept with probability min(1, exp(-β ΔE))
            if (deltaE <= 0.0 || uniform(model_.rng()) < std::exp(-beta * deltaE)) {
                model_.flipSpin(site);
                ++accepted_;
            }
            ++attempted_;
        }
    }

    /** @brief Algorithm name for output headers. */
    [[nodiscard]] static std::string name_impl() { return "Metropolis"; }

    /** @brief Acceptance ratio since last reset. */
    [[nodiscard]] double acceptanceRatio() const {
        return attempted_ > 0
            ? static_cast<double>(accepted_) / static_cast<double>(attempted_)
            : 0.0;
    }

    /** @brief Reset acceptance counters. */
    void resetCounters() { accepted_ = 0; attempted_ = 0; }

private:
    IsingModel& model_;
    long long accepted_ = 0;
    long long attempted_ = 0;
};

// ---------------------------------------------------------------------------
// CLI argument parsing (simple)
// ---------------------------------------------------------------------------

struct CLIArgs {
    int L = 16;
    double J = 1.0;
    double Jp = 0.0;
    double T = 2.269;      // Default: near Onsager Tc
    int therm_sweeps = 1000;
    int measure_sweeps = 10000;
    uint64_t seed = 0;
    bool init_random = true;
    bool auto_therm = false;  // Auto-detect thermalization
    bool time_series = false; // Output per-sweep time series
};

/**
 * @brief Parse command-line arguments.
 * @return Parsed arguments with defaults.
 */
CLIArgs parseArgs(int argc, char* argv[]) {
    CLIArgs args;
    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);
        if (arg == "--L" && i + 1 < argc)       args.L = std::atoi(argv[++i]);
        else if (arg == "--J" && i + 1 < argc)   args.J = std::atof(argv[++i]);
        else if (arg == "--Jp" && i + 1 < argc)  args.Jp = std::atof(argv[++i]);
        else if (arg == "--T" && i + 1 < argc)   args.T = std::atof(argv[++i]);
        else if (arg == "--therm" && i + 1 < argc) args.therm_sweeps = std::atoi(argv[++i]);
        else if (arg == "--sweeps" && i + 1 < argc) args.measure_sweeps = std::atoi(argv[++i]);
        else if (arg == "--seed" && i + 1 < argc) args.seed = std::strtoull(argv[++i], nullptr, 10);
        else if (arg == "--all-up")              args.init_random = false;
        else if (arg == "--auto-therm")          args.auto_therm = true;
        else if (arg == "--ts")                  args.time_series = true;
    }
    return args;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    auto args = parseArgs(argc, argv);

    // Create model: L×L Ising with NN coupling J and optional NNN Jp
    IsingModel model(args.L, args.J, args.Jp, args.seed, !args.init_random);
    model.setTemperature(args.T);

    // Create Metropolis sweep and simulation
    MetropolisSweep sweep(model);
    MCSimulation<IsingModel, MetropolisSweep> sim(model, sweep);

    // Thermalization
    int actual_therm = 0;
    if (args.auto_therm) {
        ThermConfig cfg;
        auto therm = sim.autoThermalize("abs_magnetization", cfg);
        actual_therm = therm.sweeps_used;
        std::cerr << "# auto_therm: converged=" << (therm.converged ? "true" : "false")
                  << " sweeps=" << therm.sweeps_used
                  << " mean_1st=" << therm.mean_first_half
                  << " mean_2nd=" << therm.mean_second_half << "\n";
    } else {
        for (int i = 0; i < args.therm_sweeps; ++i) {
            sweep.sweep();
        }
        actual_therm = args.therm_sweeps;
    }

    // Measurement
    auto result = sim.run(0, args.measure_sweeps,
                          /*measure_interval=*/1, args.time_series);

    // Output results to stdout (CSV with JSON header)
    sim.writeResults(result, std::cout);
    std::cerr << "# thermalization_sweeps_used: " << actual_therm << "\n";
    if (args.time_series) {
        sim.writeTimeSeries(result, std::cout);
    }

    return 0;
}
