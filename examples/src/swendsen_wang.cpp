/**
 * @file swendsen_wang.cpp
 * @brief Swendsen-Wang cluster algorithm for the 2D ferromagnetic Ising model.
 *
 * One "sweep":
 *   1. For each pair of aligned nearest neighbors, activate an FK bond
 *      with probability p = 1 - exp(-2 β J).
 *   2. Identify connected clusters via union-find.
 *   3. Flip each cluster independently with probability 1/2.
 *
 * Satisfies detailed balance for the pure ferromagnetic Hamiltonian
 * H = -J Σ_{<ij>} S_i S_j (Jp must be zero).
 *
 * @ref  Swendsen & Wang, Phys. Rev. Lett. 58, 86 (1987)
 * @complexity O(N α(N)) per sweep, where α is the inverse Ackermann function.
 *
 * Usage:
 *   swendsen_wang --L 16 --J 1.0 --T 2.5 --sweeps 10000 --therm 1000 [--seed 42]
 */

#include "../include/ising_model.hpp"
#include "../include/mc_base.hpp"
#include "../include/union_find.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace mc;

// ---------------------------------------------------------------------------
// SwendsenWangSweep — CRTP update rule
// ---------------------------------------------------------------------------

/**
 * @brief Swendsen-Wang multi-cluster sweep for IsingModel (Jp = 0 only).
 *
 * @ref  Swendsen & Wang, Phys. Rev. Lett. 58, 86 (1987)
 * @complexity O(N α(N)) per sweep.
 */
class SwendsenWangSweep : public MCSamplerCRTP<SwendsenWangSweep> {
public:
    /**
     * @brief Construct a SW sweep bound to an IsingModel.
     * @param model  Reference to the model (must outlive this object).
     */
    explicit SwendsenWangSweep(IsingModel& model)
        : model_(model),
          uf_(model.N()),
          cluster_id_(model.N(), -1),
          flip_cluster_(model.N(), false),
          bond_prob_cache_(0.0) {
        updateBondProbability();
    }

    /**
     * @brief Perform one full SW sweep: bond activation → cluster ID → flip.
     *
     * Steps:
     *   1. Reset union-find; precompute p = 1 - exp(-2βJ).
     *   2. Scan all NN bonds (right + up to avoid double-counting).
     *      For aligned pairs, activate bond with probability p.
     *   3. Identify cluster roots, assign sequential IDs.
     *   4. For each cluster, decide flip with probability 1/2.
     *   5. Bulk-flip spins belonging to flipped clusters.
     *
     * Detailed balance proof (Frame 3.5):
     *   π(S)·W(S→S') = Σ_B π(S)·P(B|S)·(1/2)^{N_C}
     *                 = Σ_B π(S')·P(B|S')·(1/2)^{N_C}   (key lemma)
     *                 = π(S')·W(S'→S)
     *
     * Key identity: e^{βJ}(1-p) = e^{-βJ}  ⟹  p = 1 - e^{-2βJ}
     */
    void sweep_impl() {
        const int N = model_.N();
        const int L = model_.L();
        updateBondProbability();

        // Step 1: Reset union-find
        uf_.reset();

        // Step 2: FK bond activation on all NN pairs
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        auto& rng = model_.rng();

        for (int site = 0; site < N; ++site) {
            const int x = site % L;
            const int y = site / L;
            const int si = model_.spin(site);

            // Right neighbor (avoid double-counting: only scan right and up)
            const int right = y * L + (x + 1 == L ? 0 : x + 1);
            if (si == model_.spin(right)) {
                if (uniform(rng) < bond_prob_cache_) {
                    uf_.unite(site, right);
                }
            }

            // Up neighbor
            const int up = ((y + 1 == L ? 0 : y + 1)) * L + x;
            if (si == model_.spin(up)) {
                if (uniform(rng) < bond_prob_cache_) {
                    uf_.unite(site, up);
                }
            }
        }

        // Step 3: Assign sequential cluster IDs from roots
        int num_clusters = 0;
        std::fill(cluster_id_.begin(), cluster_id_.end(), -1);

        for (int site = 0; site < N; ++site) {
            int root = uf_.find(site);
            if (cluster_id_[root] == -1) {
                cluster_id_[root] = num_clusters++;
            }
            // Also assign the site's cluster ID (will be resolved below)
        }

        // Resolve all cluster IDs (root → sequential ID)
        for (int site = 0; site < N; ++site) {
            int root = uf_.find(site);
            cluster_id_[site] = cluster_id_[root];
        }

        // Step 4: Decide which clusters to flip (each independently with p=1/2)
        flip_cluster_.assign(num_clusters, false);

        for (int c = 0; c < num_clusters; ++c) {
            flip_cluster_[c] = (uniform(rng) < 0.5);
        }

        // Step 5: Bulk-flip spins in flipped clusters
        auto& spins = model_.spins();
        for (int site = 0; site < N; ++site) {
            if (flip_cluster_[cluster_id_[site]]) {
                spins[site] = -spins[site];
            }
        }

        // Track diagnostics
        num_clusters_last_ = num_clusters;
        std::vector<int> cluster_sizes(num_clusters, 0);
        for (int site = 0; site < N; ++site) {
            ++cluster_sizes[cluster_id_[site]];
        }
        largest_cluster_size_last_ = *std::max_element(cluster_sizes.begin(), cluster_sizes.end());

        ++total_sweeps_;
    }

    /** @brief Algorithm name for output headers. */
    [[nodiscard]] static std::string name_impl() { return "SwendsenWang"; }

    /** @brief Number of clusters in the last sweep. */
    [[nodiscard]] int numClusters() const { return num_clusters_last_; }

    /** @brief Size of the largest cluster in the last sweep. */
    [[nodiscard]] int largestClusterSize() const { return largest_cluster_size_last_; }

    /** @brief Total sweeps performed. */
    [[nodiscard]] int totalSweeps() const { return total_sweeps_; }

private:
    IsingModel& model_;
    UnionFind uf_;
    std::vector<int> cluster_id_;
    std::vector<bool> flip_cluster_;
    double bond_prob_cache_;
    int num_clusters_last_ = 0;
    int largest_cluster_size_last_ = 0;
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

    // SW only supports pure ferromagnetic (Jp = 0)
    if (args.Jp != 0.0) {
        std::cerr << "Error: Swendsen-Wang requires Jp = 0 (pure ferromagnet).\n"
                  << "  The FK bond activation p = 1 - exp(-2βJ) is derived for\n"
                  << "  the NN-only Hamiltonian H = -J Σ S_i S_j.\n"
                  << "  For frustrated systems (Jp > 0), use Metropolis instead.\n";
        return 1;
    }

    IsingModel model(args.L, args.J, args.Jp, args.seed, !args.init_random);
    model.setTemperature(args.T);

    SwendsenWangSweep sweep(model);
    MCSimulation<IsingModel, SwendsenWangSweep> sim(model, sweep);

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
    std::cerr << "# clusters_last: " << sweep.numClusters()
              << " largest_cluster: " << sweep.largestClusterSize()
              << " / " << model.N() << " ("
              << 100.0 * sweep.largestClusterSize() / model.N() << "%)\n";
    if (args.time_series) {
        sim.writeTimeSeries(result, std::cout);
    }

    return 0;
}
