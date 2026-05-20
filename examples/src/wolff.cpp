/**
 * @file wolff.cpp
 * @brief Wolff single-cluster algorithm for the 2D ferromagnetic Ising model.
 *
 * One "sweep":
 *   1. Pick a random seed site i0.
 *   2. Grow a single cluster via BFS: for each unvisited neighbor j of
 *      the cluster frontier, if S_j = S_{i0}, add j with probability
 *      p = 1 - exp(-2 β J).
 *   3. Flip the entire cluster.
 *
 * Satisfies detailed balance for the pure ferromagnetic Hamiltonian
 * H = -J Σ_{<ij>} S_i S_j (Jp must be zero).
 *
 * @ref  U. Wolff, Phys. Rev. Lett. 62, 361 (1989)
 * @complexity O(<|C|>) per sweep, where |C| is the grown cluster size.
 *
 * Usage:
 *   wolff --L 16 --J 1.0 --T 2.5 --sweeps 10000 --therm 1000 [--seed 42]
 */

#include "../include/ising_model.hpp"
#include "../include/mc_base.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

using namespace mc;

// ---------------------------------------------------------------------------
// WolffSweep — CRTP update rule
// ---------------------------------------------------------------------------

/**
 * @brief Wolff single-cluster sweep for IsingModel (Jp = 0 only).
 *
 * Grows one cluster per sweep via BFS with FK bond activation probability
 * p = 1 - exp(-2βJ), then flips the entire cluster.
 *
 * Unlike Swendsen-Wang, no Union-Find is needed — a single visited marker
 * suffices since only one cluster is built per sweep.
 *
 * @ref  U. Wolff, Phys. Rev. Lett. 62, 361 (1989)
 * @complexity O(<|C|>) per sweep, where |C| is the cluster size.
 */
class WolffSweep : public MCSamplerCRTP<WolffSweep> {
public:
    /**
     * @brief Construct a Wolff sweep bound to an IsingModel.
     * @param model  Reference to the model (must outlive this object).
     */
    explicit WolffSweep(IsingModel& model)
        : model_(model),
          in_cluster_(model.N(), false),
          stack_(),
          bond_prob_cache_(0.0) {
        stack_.reserve(model.N());
        updateBondProbability();
    }

    /**
     * @brief Perform one Wolff sweep: seed → BFS grow → flip cluster.
     *
     * Steps:
     *   1. Pick random seed site i0.
     *   2. BFS: for each neighbor of the frontier, if aligned and
     *      uniform < p, add to cluster.
     *   3. Flip all spins in the cluster.
     *
     * Detailed balance (Frame 3.6–3.9):
     *   W(S→S') = (|C|/N) · ∏_{I(C)} p · ∏_{∂C, aligned} (1-p)
     *   Ratio cancels internal bonds; boundary bonds give:
     *     W(S→S') / W(S'→S) = exp(-β ΔE) = π(S') / π(S)
     *
     * Key identity: (1-p) = exp(-2βJ)
     */
    void sweep_impl() {
        const int N = model_.N();
        updateBondProbability();

        // Clear cluster markers
        std::fill(in_cluster_.begin(), in_cluster_.end(), false);
        stack_.clear();

        // Step 1: Pick random seed
        std::uniform_int_distribution<int> site_dist(0, N - 1);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        auto& rng = model_.rng();

        const int seed = site_dist(rng);
        const int seed_spin = model_.spin(seed);
        in_cluster_[seed] = true;
        stack_.push_back(seed);

        // Step 2: BFS cluster growth
        int idx = 0;
        while (idx < static_cast<int>(stack_.size())) {
            const int site = stack_[idx];
            ++idx;

            auto neighbors = model_.nearestNeighbors(site);
            for (int nb : neighbors) {
                if (in_cluster_[nb]) continue;
                if (model_.spin(nb) != seed_spin) continue;
                if (uniform(rng) < bond_prob_cache_) {
                    in_cluster_[nb] = true;
                    stack_.push_back(nb);
                }
            }
        }

        // Step 3: Flip the entire cluster
        for (int site : stack_) {
            model_.flipSpin(site);
        }

        // Track diagnostics
        cluster_size_last_ = static_cast<int>(stack_.size());
        ++total_sweeps_;
    }

    /** @brief Algorithm name for output headers. */
    [[nodiscard]] static std::string name_impl() { return "Wolff"; }

    /** @brief Size of the cluster in the last sweep. */
    [[nodiscard]] int clusterSize() const { return cluster_size_last_; }

    /** @brief Total sweeps performed. */
    [[nodiscard]] int totalSweeps() const { return total_sweeps_; }

private:
    IsingModel& model_;
    std::vector<bool> in_cluster_;
    std::vector<int> stack_;
    double bond_prob_cache_;
    int cluster_size_last_ = 0;
    int total_sweeps_ = 0;

    /** @brief Update the cached FK bond activation probability. */
    void updateBondProbability() {
        bond_prob_cache_ = 1.0 - std::exp(-2.0 * model_.beta() * model_.J());
    }
};

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------

struct CLIArgs {
    int L = 16;
    double J = 1.0;
    double Jp = 0.0;
    double T = 2.269;
    int therm_sweeps = 1000;
    int measure_sweeps = 10000;
    uint64_t seed = 0;
    bool init_random = true;
    bool auto_therm = false;
    bool time_series = false;
};

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

    // Wolff only supports pure ferromagnetic (Jp = 0)
    if (args.Jp != 0.0) {
        std::cerr << "Error: Wolff requires Jp = 0 (pure ferromagnet).\n"
                  << "  The FK bond activation p = 1 - exp(-2βJ) is derived for\n"
                  << "  the NN-only Hamiltonian H = -J Σ S_i S_j.\n"
                  << "  For frustrated systems (Jp > 0), use Metropolis instead.\n";
        return 1;
    }

    IsingModel model(args.L, args.J, args.Jp, args.seed, !args.init_random);
    model.setTemperature(args.T);

    WolffSweep sweep(model);
    MCSimulation<IsingModel, WolffSweep> sim(model, sweep);

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

    // Output
    sim.writeResults(result, std::cout);
    std::cerr << "# thermalization_sweeps_used: " << actual_therm << "\n";
    std::cerr << "# cluster_size_last: " << sweep.clusterSize()
              << " / " << model.N() << " ("
              << 100.0 * sweep.clusterSize() / model.N() << "%)\n";
    if (args.time_series) {
        sim.writeTimeSeries(result, std::cout);
    }

    return 0;
}
