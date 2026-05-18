#pragma once
/**
 * @file ising_model.hpp
 * @brief 2D Ising model on an L×L square lattice with periodic boundaries.
 *
 * Supports both nearest-neighbor (NN) ferromagnetic coupling J and
 * next-nearest-neighbor (NNN) frustration Jp:
 *
 *   H = -J Σ_{<ij>} S_i S_j  +  Jp Σ_{<<ik>>} S_i S_k
 *
 * where S_i = ±1, J > 0 is ferromagnetic, Jp > 0 introduces frustration.
 *
 * The class integrates with mc::ObservableRegistry for extensible measurements
 * and provides the parameterJson() method for reproducible output headers.
 *
 * @ref  Onsager, Phys. Rev. 65, 117 (1944) — exact solution for 2D square.
 * @complexity O(N) for energy computation, O(1) for single-spin energy change.
 */

#include "mc_base.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace mc {

/**
 * @brief 2D Ising model configuration with NN and optional NNN coupling.
 *
 * Stores spins as a flat vector of ±1 (row-major, index = y * L + x).
 * Provides neighbor lookup, energy computation, and observable registration.
 */
class IsingModel {
public:
    /**
     * @brief Construct an L×L Ising model.
     * @param L             Linear lattice size (L × L spins).
     * @param J             NN coupling constant (J > 0: ferromagnetic).
     * @param Jp            NNN coupling constant (Jp > 0: frustrating AFM).
     * @param seed          RNG seed (0 = use std::random_device).
     * @param init_all_up   If true, start with all spins +1; else random.
     *
     * @note Periodic boundary conditions are used in both directions.
     */
    IsingModel(int L, double J, double Jp = 0.0,
               uint64_t seed = 0, bool init_all_up = true)
        : L_(L), N_(L * L), J_(J), Jp_(Jp),
          rng_seed_(seed),
          spins_(N_, 1),  // default: all up
          rng_(seed == 0 ? std::random_device{}() : seed) {

        assert(L > 0 && "Lattice size must be positive");
        assert(J != 0.0 && "Coupling constant J must be nonzero");

        if (!init_all_up) {
            std::uniform_int_distribution<int> dist(0, 1);
            for (auto& s : spins_) {
                s = 2 * dist(rng_) - 1;  // ±1 with equal probability
            }
        }

        registerDefaultObservables();
    }

    // --- Accessors ---

    [[nodiscard]] int L() const { return L_; }
    [[nodiscard]] int N() const { return N_; }
    [[nodiscard]] double J() const { return J_; }
    [[nodiscard]] double Jp() const { return Jp_; }
    [[nodiscard]] uint64_t seed() const { return rng_seed_; }
    [[nodiscard]] double beta() const { return beta_; }

    /** @brief Set inverse temperature β = 1/(k_B T). */
    void setBeta(double beta) { beta_ = beta; }

    /** @brief Set temperature T (sets β = 1/T in natural units). */
    void setTemperature(double T) { beta_ = 1.0 / T; }

    /** @brief Access spin at site i (0-indexed, row-major). */
    [[nodiscard]] int spin(int i) const { return spins_[i]; }
    [[nodiscard]] int& spin(int i) { return spins_[i]; }

    /** @brief Access spin at 2D position (x, y) with periodic BC. */
    [[nodiscard]] int spin(int x, int y) const {
        return spins_[wrap(y) * L_ + wrap(x)];
    }
    [[nodiscard]] int& spin(int x, int y) {
        return spins_[wrap(y) * L_ + wrap(x)];
    }

    /** @brief Direct access to the spin array. */
    [[nodiscard]] const std::vector<int>& spins() const { return spins_; }
    [[nodiscard]] std::vector<int>& spins() { return spins_; }

    /** @brief Access the observable registry. */
    [[nodiscard]] ObservableRegistry& observables() { return observables_; }
    [[nodiscard]] const ObservableRegistry& observables() const { return observables_; }

    /** @brief Access the RNG (for sweep implementations). */
    [[nodiscard]] std::mt19937_64& rng() { return rng_; }

    // --- Neighbor queries ---

    /**
     * @brief Get the 4 nearest-neighbor site indices (with periodic BC).
     * @param site  Site index (row-major).
     * @return Array of 4 neighbor indices: right, left, up, down.
     */
    [[nodiscard]] std::array<int, 4> nearestNeighbors(int site) const {
        const int x = site % L_;
        const int y = site / L_;
        const int xr = wrap(x + 1);
        const int xl = wrap(x - 1);
        const int yu = wrap(y + 1);
        const int yd = wrap(y - 1);
        return {y * L_ + xr, y * L_ + xl, yu * L_ + x, yd * L_ + x};
    }

    /**
     * @brief Get the 4 next-nearest-neighbor site indices (diagonal).
     * @param site  Site index (row-major).
     * @return Array of 4 NNN indices: NE, NW, SE, SW.
     */
    [[nodiscard]] std::array<int, 4> nextNearestNeighbors(int site) const {
        const int x = site % L_;
        const int y = site / L_;
        return {
            wrap(y + 1) * L_ + wrap(x + 1),  // NE
            wrap(y + 1) * L_ + wrap(x - 1),  // NW
            wrap(y - 1) * L_ + wrap(x + 1),  // SE
            wrap(y - 1) * L_ + wrap(x - 1),  // SW
        };
    }

    // --- Energy computation ---

    /**
     * @brief Compute the total energy of the current configuration.
     * @return Total energy E = -J Σ_{NN} S_iS_j + Jp Σ_{NNN} S_iS_k.
     * @complexity O(N). Each bond counted once.
     */
    [[nodiscard]] double computeEnergy() const {
        double energy = 0.0;
        for (int i = 0; i < N_; ++i) {
            const int si = spins_[i];
            // NN bonds: only count right and up to avoid double-counting
            const auto [right, left, up, down] = nearestNeighbors(i);
            energy -= J_ * si * spins_[right];
            energy -= J_ * si * spins_[up];

            // NNN bonds: count NE and NW to avoid double-counting
            if (Jp_ != 0.0) {
                const auto [ne, nw, se, sw] = nextNearestNeighbors(i);
                energy += Jp_ * si * spins_[ne];
                energy += Jp_ * si * spins_[nw];
            }
        }
        return energy;
    }

    /**
     * @brief Compute energy change ΔE from flipping a single spin.
     * @param site  The site to flip.
     * @return Energy change (always an even multiple of J or Jp for Ising).
     * @complexity O(1). Only depends on the 4+4 neighbors.
     *
     * @note Flipping S_i → -S_i changes each bond contribution by ±2J.
     *       ΔE = 2J Σ_{nn} S_i S_j  +  2Jp Σ_{nnn} S_i S_k
     *       (positive means the flip costs energy).
     */
    [[nodiscard]] double energyChange(int site) const {
        const int si = spins_[site];
        double sum_nn = 0;
        for (int nb : nearestNeighbors(site)) {
            sum_nn += si * spins_[nb];
        }

        double sum_nnn = 0;
        if (Jp_ != 0.0) {
            for (int nb : nextNearestNeighbors(site)) {
                sum_nnn += si * spins_[nb];
            }
        }

        return 2.0 * J_ * sum_nn + 2.0 * Jp_ * sum_nnn;
    }

    /** @brief Flip spin at site i. */
    void flipSpin(int i) { spins_[i] = -spins_[i]; }

    // --- Observables ---

    /**
     * @brief Compute magnetization per spin.
     * @return m = (1/N) Σ_i S_i, range [-1, 1].
     */
    [[nodiscard]] double magnetization() const {
        int sum = 0;
        for (int s : spins_) sum += s;
        return static_cast<double>(sum) / static_cast<double>(N_);
    }

    /**
     * @brief Compute absolute magnetization per spin.
     * @return |m| — used to avoid symmetry-broken cancellation.
     */
    [[nodiscard]] double absMagnetization() const {
        return std::abs(magnetization());
    }

    /**
     * @brief Compute energy per spin.
     * @return e = E / N.
     */
    [[nodiscard]] double energyPerSpin() const {
        return computeEnergy() / static_cast<double>(N_);
    }

    // --- Output helpers ---

    /**
     * @brief Serialize model parameters as a JSON string.
     *
     * Used as the header comment in output CSV files for reproducibility.
     */
    [[nodiscard]] std::string parameterJson() const {
        std::ostringstream oss;
        oss << "{\"L\":" << L_
            << ",\"J\":" << J_
            << ",\"Jp\":" << Jp_
            << ",\"beta\":" << beta_
            << ",\"N\":" << N_
            << "}";
        return oss.str();
    }

private:
    int L_;                   ///< Linear lattice size
    int N_;                   ///< Total spins = L × L
    double J_;                ///< NN coupling (J > 0: ferromagnetic)
    double Jp_;               ///< NNN coupling (Jp > 0: frustrating)
    double beta_ = 1.0;       ///< Inverse temperature (default β = 1)
    uint64_t rng_seed_;       ///< Stored for reproducibility output
    std::vector<int> spins_;  ///< Spin array: each entry is ±1
    std::mt19937_64 rng_;     ///< Per-model RNG (injectable, not global)
    ObservableRegistry observables_;

    /** @brief Wrap coordinate into [0, L_) with periodic BC. */
    [[nodiscard]] int wrap(int coord) const {
        int r = coord % L_;
        return r < 0 ? r + L_ : r;
    }

    /** @brief Register the default set of observables (m, |m|, e). */
    void registerDefaultObservables() {
        observables_.registerObservable("magnetization",
            [](const void* ptr) -> double {
                return static_cast<const IsingModel*>(ptr)->magnetization();
            });
        observables_.registerObservable("abs_magnetization",
            [](const void* ptr) -> double {
                return static_cast<const IsingModel*>(ptr)->absMagnetization();
            });
        observables_.registerObservable("energy_per_spin",
            [](const void* ptr) -> double {
                return static_cast<const IsingModel*>(ptr)->energyPerSpin();
            });
    }
};

}  // namespace mc
