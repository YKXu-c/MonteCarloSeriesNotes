# Tech Tutorial — C++23 MC Framework

Auto-generated reference for C++ syntax, patterns, and tools used in this project.
(For self-study; Claude writes this but does not read it back.)

---

## 1. C++23 Features Used

### `std::print` / `std::println`
```cpp
#include <print>
std::print("L = {}, T = {}\n", L, T);
```
Replaces `std::cout` and `printf` with type-safe, format-string output.
Available since C++23. Requires `<print>` header.

### Structured bindings with `auto& [name, value]`
```cpp
for (auto& [name, fn] : observables_) { ... }
```
Decomposes `std::pair` or `std::tuple` elements. The `auto&` ensures
we iterate by reference (no copy).

### `[[nodiscard]]` attribute
```cpp
[[nodiscard]] int L() const { return L_; }
```
Tells the compiler to warn if the return value is ignored.
Prevents accidental `model.L();` when you meant `int x = model.L();`.

### `[[nodiscard]] std::string` on `static` functions
```cpp
[[nodiscard]] static std::string name_impl() { return "Metropolis"; }
```
`static` member functions don't need an object instance. Used in CRTP
where the derived class provides its name via a static method.

---

## 2. CRTP Pattern (Curiously Recurring Template Pattern)

```cpp
template <typename Derived>
class Base {
public:
    void doSomething() {
        static_cast<Derived*>(this)->doSomething_impl();
    }
};

class Concrete : public Base<Concrete> {
public:
    void doSomething_impl() { /* actual work */ }
};
```

**Why:** Zero-overhead polymorphism. No virtual table, no indirect call.
The compiler sees exactly which `impl` function to call at compile time.

**Cost:** Same performance as hand-writing the loop without any abstraction.

**Trade-off:** Cannot store heterogeneous containers of different Derived types
(without `std::variant` or type erasure).

---

## 3. Random Number Generation in Modern C++

```cpp
#include <random>

// 64-bit Mersenne Twister — high-quality, long period (2^19937 - 1)
std::mt19937_64 rng(seed);

// Uniform integer in [0, N-1]
std::uniform_int_distribution<int> site_dist(0, N - 1);
int site = site_dist(rng);

// Uniform real in [0, 1)
std::uniform_real_distribution<double> uniform(0.0, 1.0);
double r = uniform(rng);
```

**Key points:**
- `std::mt19937_64` is the recommended RNG for scientific computing
- Distribution objects are lightweight — create locally in hot loops
- Pass `rng` by reference (not by value!) to avoid reseeding
- `std::random_device{}()` gets a true random seed from the OS

---

## 4. Compile and Run

```bash
# Compile (header-only, no linking needed)
g++ -std=c++23 -Wall -Wextra -I include -o build/metropolis src/metropolis.cpp

# Run
./build/metropolis --L 16 --J 1.0 --T 2.5 --therm 2000 --sweeps 5000 --seed 42
```

**Flag explanations:**
- `-std=c++23` — enable C++23 features
- `-Wall -Wextra` — enable most compiler warnings
- `-I include` — add `include/` to header search path
- `-o build/metropolis` — output binary location

---

## 5. Design Pattern: Observable Registry

```cpp
// Register a named measurement function
observables_.registerObservable("magnetization",
    [](const void* ptr) -> double {
        return static_cast<const IsingModel*>(ptr)->magnetization();
    });

// Measure all registered observables at once
auto results = observables_.measureAll(&model);
```

**Why `const void*`:** The registry doesn't know the model type at
registration time. The lambda captures the concrete type via `static_cast`.
This allows adding new observables (e.g. "cluster_overlap" in Lecture 2)
without modifying the framework code.

---

*More entries will be appended as new algorithms and techniques are introduced.*
