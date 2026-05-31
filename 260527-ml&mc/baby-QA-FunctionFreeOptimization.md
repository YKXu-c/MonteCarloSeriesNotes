# Baby-QA: Function-free Optimization via Comparison Oracles (arXiv 2604.26867v2)

## Paper Overview
- **Authors:** Katya Scheinberg, Zikai Xiong (Georgia Tech)
- **Topic:** Optimization using only a comparison oracle (no function values, no gradients)
- **Key idea:** Replace function-based analysis with preference-relation geometry (level sets)

---

## Terminology (CN-EN)

| English | 中文 |
|---------|------|
| Comparison oracle | 比较预言机 |
| Preference relation | 偏好关系 |
| Function-free optimization | 无函数优化 |
| Level set | 水平集 |
| Sublevel set | 次水平集 |
| Normal direction | 法线方向 |
| Tangent direction | 切线方向 |
| Orthogonal complement | 正交补 |
| Regularity radius | 正则半径 |
| Level-set optimality gap | 水平集最优性差距 |
| Plateau-free | 无平台 |
| Convex preference | 凸偏好 |
| Regular | 正则 |
| Step size | 步长 |
| Comparison radius | 比较半径 |
| Line search | 线搜索 |
| Normal Direction Descent (NDD) | 法线方向下降法 |
| Normalized gradient descent | 归一化梯度下降 |
| Cutting-plane method | 割平面法 |
| Lower bound | 下界 |
| Nearly optimal | 几乎最优 |
| Quasiconvex | 拟凸 |
| Preference-equivalent | 偏好等价 |
| Representative function | 代表函数 |
| Best-so-far | 最好到目前为止 |
| Two-case guarantee | 双保证 |
| Zero chain | 零链 |
| Complete and transitive | 完备且传递 |
| Strict preference | 严格偏好 |
| Indifference | 无差异 |
| Feasible set | 可行集 |
| Optimal solution set | 最优解集 |
| Stationary point | 驻点/平稳点 |

---

## Unit-by-Unit Q&A

### Unit 0: Big Picture
- **Q:** 传统优化坐滑梯、找梯度。比较 oracle 只有差然后指出方向
- **A:** 正确。"坐滑梯"是 gradient descent，比较 oracle 只有 {−1, 0, +1} 的信息
- Three core questions: Q1 (framework), Q2 (normal direction estimation), Q3 (optimization algorithm)

### Unit 1: Preference Relation & Comparison Oracle
- **Key:** x ≼ y means "x is no worse than y" (x preferred)
- Completeness: any two points can be compared
- Transitivity: no cycles (A>B>C>A)
- Oracle: 𝒞_cmp(x,y) ∈ {−1, 0, +1}
- **Q:** f(x) = x₁² + x₂², X* = ? 𝒞_cmp((1,1),(0,2)) = ?
- **A:** X* = {(0,0)}, 𝒞_cmp = +1 (since (1,1) ≺ (0,2))

### Unit 2: Level Sets & Sublevel Sets
- **Key definitions:**
  - Sₓ = {y : y ≼ x} (sublevel set, all points no worse than x)
  - Lₓ = {y : y ∼ x} (level set, all points equally good as x)
- **Example:** f(x) = ‖x‖₂, x = (1,0): Sₓ = unit disk, Lₓ = unit circle
- **Three properties:**
  - Plateau-free: Lₓ has no interior (dim < d)
  - Convex: every Sₓ is convex
  - Regular: Lₓ is smooth near x (unique normal direction exists)
- Regularity guarantees nₓ exists; Convexity guarantees global optimality; Plateau-free ensures Lₓ is "thin"
- **Q:** Can plateau-free be understood as dim(L) < dim(X)?
- **A:** Essentially yes — Lₓ is typically a (d−1)-dimensional hypersurface with no interior

### Unit 3: Optimality Measures
- **Δ_LS(x) = dist(X*, Lₓ):** distance from level set to optimal set
- **rₓ:** regularity radius (how far the smooth region extends)
- Δ_LS ≤ dist(x, X*) always (because x ∈ Lₓ, picking the best point in Lₓ can only help)
- **Key insight:** Δ_LS is an analytical tool, not computed by the algorithm
- rₓ small ↔ near-stationary (analogous to small gradient)
- **Q:** How can we know Δ_LS without knowing X*?
- **A:** We don't need to. Δ_LS is only used in convergence proofs, not during algorithm execution
- **Q:** Why is Δ_LS more loose than dist(x, X*)?
- **A:** Lₓ = "all points as good as me" — picking the closest one to X* is never farther than using just x

### Unit 4: Relation to Function-based Optimization
- Preference-equivalent functions: f and g produce same comparisons ⟺ g = φ∘f for strictly increasing φ
- Example: f₁=‖x‖, f₂=log(‖x‖+1), f₃ (piecewise) — all same level sets (spheres), same comparisons
- **Key philosophy:** All assumptions stated in terms of preference relation, not any specific f
- **Q:** "I optimized a non-convex function using comparison oracle" — response?
- **A:** You may have chosen a non-convex representative f, but the preference relation itself might be convex (quasiconvex). The statement is ill-defined without specifying which f.

### Unit 5: Normal Direction Estimation (Algorithm 1 & 2)
- **Core method:** Orthogonal-complement construction
  1. Find d−1 tangent directions (where 𝒞_cmp = 0)
  2. Take orthogonal complement → normal direction
- Binary search to find each tangent direction: O(log(d/ε)) comparisons
- Total: O(d·log(d/ε)) comparisons — nearly optimal (lower bound: Ω(d·log(1/ε)))
- **Q:** Why not use "direction of largest improvement"?
- **A:** Comparison oracle only gives {−1, 0, +1} — no magnitude information. Can't tell "how much better"

### Unit 6: NDD Algorithm
- **Update:** xₖ₊₁ = Π_C(xₖ − η·n̂ₖ) (projected normalized descent)
- Same as normalized GD when f exists and is differentiable
- **Two-case guarantee:** (i) Δ_LS small, or (ii) rₓ small (near-stationary)
- **Complexity:** Õ(d·D²/ε²) comparisons, D = dist(x₁, X*)
- **Lower bound:** Ω(D²/ε²) normal direction steps — NDD is nearly optimal
- **adaNDD:** Adaptive step size (line search) + adaptive comparison radius, same complexity without knowing parameters
- **Cutting plane vs NDD:** Cutting plane Õ(d²·log(R/ε)) better for small d/high accuracy; NDD Õ(d·D²/ε²) better for large d

### Key Discussion Points

#### h vs rₓ
- rₓ: intrinsic property of preference at x (radius of smooth region)
- h: algorithm's probing distance (comparison radius)
- Need h < rₓ for reliable normal estimation
- Adaptive: check comparison consistency; if inconsistent, shrink h

#### Step size η vs Comparison radius h
- η: "how far to walk" along −n̂ₖ (controlled by line search in adaNDD)
- h: "how far to probe" when estimating normal (controlled by adaptive mechanism)

#### Convexity requirement
- Paper's "convex preference" = all Sₓ are convex = quasiconvex
- Weaker than convex function: f₂ = log(‖x‖+1) is not convex but Sₓ is convex
- Two pits (multiple local minima) → violates convexity (Sₓ is dumbbell-shaped)
- Cannot fix by changing representative function (level set geometry is invariant)
- Same limitation as classical convex optimization, but with weaker oracle

#### Connection to spin frustrated systems
- Spin frustrated systems: discrete configuration space, exponentially many local minima
- Does NOT satisfy convex or continuous assumptions
- Paper's framework cannot directly apply
- Methodology (level-set geometry analysis) may inspire new approaches

---

## Paper Structure Summary

| Section | Content | Key Result |
|---------|---------|------------|
| Sec 1 | Introduction, motivation, 3 questions | Framework overview |
| Sec 2 | Function-free optimization framework | Δ_LS, rₓ, level set geometry |
| Sec 3 | Normal direction estimation | O(d·log(d/ε)) comparisons |
| Sec 4 | NDD + adaNDD + lower bounds | Õ(dD²/ε²) comparisons |
| Sec 5 | Summary, future directions | Non-convex extension open |

## Complexity Summary

| Quantity | Upper bound | Lower bound |
|----------|-------------|-------------|
| Normal direction estimation (comparisons) | O(d·log(d/ε)) | Ω(d·log(1/ε)) |
| NDD normal direction steps | O(D²/ε²) | Ω(D²/ε²) |
| NDD total comparisons | Õ(d·D²/ε²) | — |
