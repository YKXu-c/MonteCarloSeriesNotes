# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A lecture series titled "MonteCarlo: From classic to quantum" — Beamer slides progressing from classical spin-system Monte Carlo algorithms to quantum impurity QMC, with connections to DMFT and dynamical-system analysis of algorithms. The project plan is maintained in Chinese at `编写计划与原则.md`.

**Lecture series scope** (7 lectures planned):
1. Classical MC Algorithms (Sections 1–3 written; Sec 4 C++ code in progress — Metropolis + SW + Wolff done; Sections 5–7 placeholders)
2. MC Algorithm as a Dynamic System
3a. Application in Kondo-Ising model (DFT+U fitting, BO, CMA-ES)
3b. Monte Carlo Results
4a. QMC for Impurity Problems (Metropolis)
4b. VMC, dQMC, DMC, FCIQMC

## Structure

- `MonteCarlo_From_classic_to_quantum_1.tex` — Beamer slides for lecture 1 (classical MC algorithms, 30 frames)
- `MonteCarlo_From_classic_to_quantum_2.tex` — Placeholder for lecture 2 (empty until lecture 1 is finished)
- `generatePic.py` — Generates figures to `pics/` at 200 DPI. Two categories: (1) drawing functions (e.g. `kagome_kondo_ising()`) for schematic diagrams, (2) data plotting infrastructure (`parse_mc_output`, `run_mc_binary`, `temperature_sweep`, `onsager_exact_M`, `plot_metropolis_mt`, `plot_thermalization`, `plot_sw_mt`, `plot_wolff_mt`, `plot_dynamic_exponent` (τ_int vs L, log-log fit for z)) that runs C++ binaries and plots numerical results. Add new functions and call them from `__main__`. Each function saves to `pics/$name.png` via `fig.savefig()`.
- `pics/` — Generated images and reference images. Generated figures go here; reference images (e.g., `TbCo2*.png`, `SW_example*.png`) are manually placed.
- `refers/` — Reference papers, named `[arXivID]Title.pdf` or `[legacyID]Title.pdf`.
- `examples/` — C++ example code. Infrastructure: `include/mc_base.hpp` (CRTP framework + `ObservableRegistry` + auto-thermalization via `autoThermalize()` with sliding-window plateau detection), `include/ising_model.hpp` (2D Ising model), `include/union_find.hpp` (path compression + rank union-find for cluster algorithms), `src/metropolis.cpp` (Metropolis implementation with `--auto-therm`, `--ts` time series output), `src/swendsen_wang.cpp` (SW cluster algorithm with FK bond activation), `src/wolff.cpp` (Wolff single-cluster algorithm). Uses `spack` for package management and CMake for builds.
- `编写计划与原则.md` — Project plan in Chinese. Progress tracking table at the bottom records each iteration (date, content, notes). This table is the authoritative record of what has been completed.

## Commands

**Python environment:** `mc_learn` virtual environment (Python 3.12), located at `ykxu/uenv/mc_learn`. Requires `matplotlib` and `numpy`. Activate before running `generatePic.py`.

**Spack environment (for C++ dependencies):**
```
spack install googletest googlebenchmark cppcheck
spack load googletest googlebenchmark
```

**Compile slides:**
```
xelatex MonteCarlo_From_classic_to_quantum_1.tex
```
Build produces auxiliary files (`.aux`, `.log`, `.nav`, `.out`, `.snm`, `.toc`) — these can be safely deleted.

**Build C++ examples:**
```
cd examples && mkdir -p build && cd build && cmake .. && make
```
Run with `./metropolis --L 16 --T 2.0 --sweeps 10000` or `./swendsen_wang` (same flags). For the full CLI argument table (`--auto-therm`, `--ts`, `--all-up`, `--seed`, `--Jp`, `--therm`), see `examples/README.md`.

**Generate figures:**
```
python generatePic.py
```
Note: `generatePic.py` runs C++ binaries from `examples/build/` via `run_mc_binary()`, so build the C++ examples first.

## Beamer Conventions

**Theme:** Madrid theme, whale color theme, custom blue palette:
- `mainblue` #004983, `darkblue` #002F5F, `accentblue` #4682B4, `lightgray` #F0F0F0

**Additional packages:** `algorithm2e`, `bm`. Graphics path set to `pics/` via `\graphicspath{{pics/}}`.

**TikZ libraries:** `arrows.meta`, `positioning`, `calc`, `decorations.pathreplacing`. Flowcharts use mermaid.

**Lecture layout:** 1. Motivation/model → 2. Derivation → 3. Algorithm examples → 4. Discussion. Section numbers are topic labels, not page numbers. Keep main slide text minimal — use flowcharts and key equations. Detailed derivations go on backup slides after the main page.

**Lecture 1 `.tex` section status:**
- Sec 1 (Motivation: TbCo₂ & Kondo-Ising) — **complete** (6 frames)
- Sec 2 (Monte Carlo Fundamentals) — **complete** (4 frames + backup)
- Sec 3 (Classical Algorithms: Metropolis, SW, Wolff) — **complete** (13 frames + backup derivation slides, comparison table, dynamic exponent z measurement)
- Sec 4 (Heisenberg generalization) — **C++ code in progress** (base classes + Metropolis + SW + Wolff done)
- Sec 5 (Wolff on Infinite Lattice) — **placeholder**
- Sec 6 (Loop & Worm) — **placeholder**
- Sec 7 (Summary & Outlook) — **placeholder**

## Key Principles

- All content except the plan document is written in **English**.
- Academic rigor: every claim must have a citable reference; include mathematical formulations where applicable.
- References go in `refers/` as `[arXivID]Title.pdf`. Prefer arXiv if behind paywalls.
- Web searches use the **MiniMax MCP tool** (`mcp__MiniMax__web_search`). Download found references to `refers/`.
- Each lecture plan section `0.` lists figures and references to prepare before writing slides.
- Parentheses `（）` in the plan denote operational instructions, not content.

## Development Workflow

**Small-step iteration:** Complete one self-contained unit (a figure, a slide section, or example code) → verify → stop → record progress in the table at the bottom of `编写计划与原则.md` (date, completed content, notes).

**Verification:** After each compile, use MiniMax image understanding (`mcp__MiniMax__understand_image`) to check the rendered PDF for text overlap, content overflow beyond slide/block boundaries, and visual correctness. This has caught numerous issues in past iterations.

**Cycle:** After each unit, run `/clear` + `/init`. The first action after `/init` should ask whether the previous iteration had any issues or corrections.

**Mathematical accuracy:** Past iterations have required corrections to detailed-balance proofs (SW and Wolff). When adding or modifying proofs, double-check: (1) transfer probability formulas match the algorithm's actual mechanism, (2) FK identity $e^{\beta J}(1-p) = e^{-\beta J}$ is applied correctly, (3) boundary bond contributions are accounted for properly.

**MCP tools:** MiniMax image understanding (`mcp__MiniMax__understand_image`) and web search (`mcp__MiniMax__web_search`) are pre-allowed in `.claude/settings.local.json`.

## C++ Example Code Conventions

- **Directory structure:** `include/` (headers), `src/` (implementations), `tests/` (GoogleTest), `benchmarks/` (Google Benchmark), `tools/` (Python utilities)
- **Naming:** PascalCase classes, camelCase functions, snake_case variables, trailing `_` for members, UPPER_SNAKE constants
- **Comments:** Doxygen-compatible (`@brief`, `@param`, `@return`, `@ref`, `@complexity`). Every class and public function documented.
- **Architecture:** CRTP for compile-time polymorphism. `MCSimulation<ModelType, UpdateRule>` decouples algorithms from models. `ObservableRegistry` for extensible measurements.
- **Reproducibility:** `std::mt19937_64` RNG, JSON parameter headers in output, `--checkpoint`/`--resume` support
- **Error handling:** `assert()` for preconditions, tolerance-based float comparison, verbose logging (`-v`, `-vv`, `-vvv`)
- **Performance:** Pre-allocated buffers in hot loops, SoA layout consideration, union-find with path compression + rank
- **Testing:** GoogleTest unit tests, regression tests against known solutions (Onsager), Valgrind/ASan
- **Build:** CMake or direct g++ (Apple Clang 21), C++23, spack for dependencies (Eigen, GoogleTest, Google Benchmark, cppcheck)
- Code syntax, algorithms, and tool usage are documented in `examples/techTutorial.md` for user learning. Claude generates this file but does not need to read it.