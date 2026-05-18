# Baby QA — MonteCarloSeriesNotes Teaching Session

**Date:** 2026-05-15 ~ 2026-05-18
**Mode:** Baby (保姆级)
**Files taught:** mc_base.hpp, ising_model.hpp, metropolis.cpp, generatePic.py
**Languages compared:** C++, Python, Rust (user not familiar with Rust syntax but interested in comparisons)

---

## Teaching Progress

### mc_base.hpp (Framework)

| Unit | Topic | Key Concepts |
|------|-------|-------------|
| 1 | File header + includes + namespace | `#pragma once` vs `#ifndef`/`#define`/`#endif`, namespace ≈ Python module |
| 2 | ObservableRegistry | `std::function`, `void*` type erasure, `std::move` (ownership transfer), `const` placements, `std::map` vs hash table, `[[nodiscard]]` |
| 3 | MCSamplerCRTP | CRTP = Curiously Recurring Template Pattern, `static_cast<Derived*>(this)`, virtual vs CRTP, `[[nodiscard]]` |
| 4 | MCResult | `struct` vs `class`, `<O^2> - <O>^2` variance, `uint64_t` |
| 5 | ThermConfig + ThermResult | Auto-thermalization: sliding window, first/second half comparison |
| 6 | MCSimulation constructor | Template with two params (ModelType, UpdateRule), initializer list, reference members (borrowing, not owning) |
| 7 | MCSimulation::run() | Thermalization loop → measurement loop → accumulate mean/mean2 → divide by count |
| 8 | writeResults + writeTimeSeries | `std::ostream` polymorphism, tab-separated output, `# time_series_begin/end` markers |
| 9 | autoThermalize() | `.reserve()`, batch loop, sliding window convergence, early return on convergence |

### ising_model.hpp (Model)

| Unit | Topic | Key Concepts |
|------|-------|-------------|
| 10 | Constructor + accessors | Member init list, ternary operator, `std::random_device`, lambda for observable registration, `const` vs non-const overloads, `explicit` |
| 11 | Neighbors + energy | `std::array<int,4>`, periodic boundary conditions (PBC) via `wrap()`, avoiding double-counting bonds, `energyChange` O(1) formula: Delta_E = 2J sum_nn + 2Jp sum_nnn |
| 12 | Observables + output | `std::ostringstream`, `parameterJson()` for reproducibility, `std::mt19937_64` |

### metropolis.cpp (Algorithm + main)

| Unit | Topic | Key Concepts |
|------|-------|-------------|
| 13 | MetropolisSweep | CRTP subclass, `explicit` constructor, Metropolis acceptance criterion, `uniform_int_distribution` / `uniform_real_distribution`, acceptance ratio |
| 14 | CLI + main | `argc`/`argv`, `std::atoi`/`std::atof`, `auto`, `std::cout` vs `std::cerr`, return code |

### generatePic.py (Visualization)

| Unit | Topic | Key Concepts |
|------|-------|-------------|
| 15 | kagome_kondo_ising() | matplotlib basics, kagome lattice construction from triangular lattice basis, numpy array operations |
| 16 | MC data parsing | `subprocess.run()`, JSON parsing, tab-separated parsing, `**kwargs`, compiled binary runs without environment |
| 17 | Onsager exact + plotting | `onsager_exact_M()`, temperature sweep, `np.sinh`, `np.where`, `np.errstate`, matplotlib LaTeX rendering |

---

## Key Questions & Corrections

### Q: `#pragma once` vs traditional include guard?
**A:** Both prevent double inclusion. `#pragma once` is one line, compiler extension. Traditional is `#ifndef`/`#define`/`#endif`, standard-compliant. Project uses `#pragma once`.

### Q: `namespace mc` ≈ Python module?
**A:** Yes. `namespace mc { }` adds a `mc::` prefix. `using namespace mc` ≈ `from module import *`. C++ namespaces are not tied to files (unlike Python).

### Q: Template equivalents in Rust and Fortran?
**A:** Rust has generics (compile-time monomorphization, like C++). Fortran has limited parameterized derived types (kind/len only, cannot parameterize behavior). CRTP is not possible in Fortran.

### Q: What is `trait` in Rust?
**A:** A list of capabilities a type must implement. ≈ Python ABC, C++20 concepts, Java interfaces. Constrains generics: `fn max<T: Comparable>(a: T, b: T)`.

### Q: What is `impl`?
**A:** Rust keyword for "implement". Two uses: `impl Type { methods }` (inherent methods ≈ Python class methods), `impl Trait for Type { }` (satisfy a trait).

### Q: Difference between sweep and step?
**A:** Step = one single spin-flip attempt. Sweep = N steps (N = total lattice sites). Sweep is the standard MC time unit so different system sizes are comparable.

### Q: How is sweep count N determined?
**A:** N = L × L. No physics derivation, just convention: "visit every site once on average."

### Q: Do SW and Wolff also have sweep?
**A:** SW: one sweep = build all clusters, flip all (updates all N sites). Wolff: no natural sweep definition; one step = flip one cluster. Can define "sweep" as N/<|C|> steps for normalization.

### Q: `std::move` — why not reference?
**A:** Reference borrows (temporary). Map needs to own the function permanently. If the original is destroyed, the reference becomes dangling. `std::move` transfers ownership — map now owns it.

### Q: Ownership concept?
**A:** Who is responsible for freeing memory. Copy = duplicate (2 owners). Move = transfer (1 owner). Reference = borrow (0 new owners, original responsible). Rust enforces this at compile time; C++ relies on programmer discipline.

### Q: CRTP full name? Difference from normal inheritance?
**A:** Curiously Recurring Template Pattern. Normal inheritance uses virtual functions (runtime vtable lookup). CRTP uses templates (compile-time binding, zero overhead). Trade-off: CRTP cannot store heterogeneous types in one container.

### Q: Why `explicit`?
**A:** Prevents implicit conversions. Without `explicit`, single-arg constructors allow the compiler to silently construct objects when you pass the wrong type. With `explicit`, you must explicitly construct. Rule: always add `explicit` to single-argument constructors.

### Q: Can compiled C++ binary run without environment?
**A:** Mostly yes — it's machine code, no interpreter needed. But depends on platform (macOS binary won't run on Linux) and dynamically linked libraries (libc, libc++). Python needs its interpreter + all pip packages.

---

## Cross-Language Reference Table

| Concept | C++ | Python | Rust |
|---------|-----|--------|------|
| Include/import | `#include` (text paste) | `import` (module) | `use` (module) |
| Include guard | `#pragma once` | Not needed | Not needed |
| Namespace | `namespace mc { }` | File = module | `mod mc { }` |
| Class | `class` / `struct` | `class` | `struct` + `impl` |
| Inheritance | `class A : public B` | `class A(B)` | No inheritance, use traits |
| Virtual | `virtual` + vtable | All methods virtual | `dyn Trait` |
| Static dispatch | CRTP | N/A | `impl Trait` / generics |
| Ownership | Manual / `std::move` | GC handles it | Compile-time enforced |
| Borrowing | `const &` / `&` | All variables are references | `&T` / `&mut T` |
| Template/generic | `template<typename T>` | Duck typing | `<T: Trait>` |
| Null safety | No (undefined behavior) | No (runtime exceptions) | Yes (Option<T>) |
| Error handling | Exceptions / assert | Exceptions | Result<T, E> |
| Header files | `.hpp` / `.h` | None | None |
| Compiled binary | Direct execution | Needs interpreter | Direct execution |

---

## Math Concepts Referenced

- **Detailed balance:** Metropolis acceptance min(1, exp(-βΔE))
- **Onsager exact solution:** M(T) = (1 - sinh^-4(2J/T))^(1/8) for T < Tc
- **Critical temperature:** Tc = 2/ln(1+sqrt(2)) ≈ 2.269 (2D Ising square, NN only)
- **Variance:** Var(O) = <O²> - <O>² (stored as mean2 - mean²)
- **Percolation:** Bond percolation with probability p; percolation phase transition at pc
- **FK bond representation:** p = 1 - exp(-2βJ), connects cluster algorithms to percolation
- **Periodic boundary conditions:** `wrap()` with modular arithmetic, simulates infinite system
