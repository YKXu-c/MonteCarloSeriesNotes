# Monte Carlo Examples

C++23 implementations of classical Monte Carlo algorithms for the 2D Ising model.

## Build Requirements

- **Compiler:** GCC 13+ or Clang 17+ (C++23 support required)
- **Build tools:** CMake 3.25+ (optional, direct compilation also works)
- **Dependencies (via spack):**
  - `googletest` — unit testing
  - `googlebenchmark` — performance benchmarks
  - `cppcheck` — static analysis

## Quick Build (without CMake)

```bash
# Single-file compilation (header-only framework)
g++ -std=c++23 -Wall -Wextra -I include -o build/metropolis src/metropolis.cpp

# Run a quick test
./build/metropolis --L 16 --J 1.0 --T 2.5 --therm 2000 --sweeps 5000
```

## Build with CMake

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Command-Line Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `--L` | 16 | Linear lattice size (L × L) |
| `--J` | 1.0 | NN coupling (J > 0: ferromagnetic) |
| `--Jp` | 0.0 | NNN coupling (Jp > 0: frustrating) |
| `--T` | 2.269 | Temperature |
| `--therm` | 1000 | Thermalization sweeps |
| `--sweeps` | 10000 | Measurement sweeps |
| `--seed` | 0 | RNG seed (0 = random) |
| `--all-up` | — | Start from all spins up |

## Output Format

Tab-separated with JSON parameter header:
```
# {"L":16,"J":1,"Jp":0,"beta":0.4,"N":256}
# algorithm: Metropolis
# thermalization_sweeps: 2000
# measurement_sweeps: 5000
# seed: 12397089208
observable	mean	variance	mean2
magnetization	0.52	0.03	0.271
```

## Spack Environment Setup

```bash
# Install dependencies
spack install googletest
spack install googlebenchmark
spack install cppcheck

# Load into environment
spack load googletest googlebenchmark
```
