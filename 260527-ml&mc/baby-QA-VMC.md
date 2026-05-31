# Baby-QA: Variational Monte Carlo (VMC) + Optimization Methods

## Overview
- **Topic**: VMC 方法、优化算法、FFO-CO 集成方案
- **Connection**: FFO-CO 论文 (arXiv 2604.26867v2) 的应用场景——VMC 参数优化中 MC 噪声使梯度不可靠，但 sign(ΔE) 更稳健

---

## Terminology (CN-EN)

| English | 中文 |
|---------|------|
| Variational Monte Carlo (VMC) | 变分蒙特卡洛 |
| Variational principle | 变分原理 |
| Trial wavefunction / ansatz | 试探波函数 |
| Local energy E_loc | 局域能量 |
| Configuration space | 构型空间 |
| Parameter space | 参数空间 |
| Stochastic Gradient Descent (SGD) | 随机梯度下降 |
| Stochastic Reconfiguration (SR) | 随机重整化 |
| Natural gradient | 自然梯度 |
| Quantum Geometric Tensor (QGT) | 量子几何张量 |
| Fisher information matrix | Fisher 信息矩阵 |
| Linear Method (LM) | 线性方法 |
| Exact Diagonalization (ED) | 精确对角化 |
| Jastrow factor | Jastrow 关联因子 |
| Slater determinant | Slater 行列式 |
| Neural Network Quantum State (NQS) | 神经网络量子态 |
| Signal-to-noise ratio (SNR) | 信噪比 |
| Statistical fluctuation | 统计涨落 |
| Covariance | 协方差 |
| Convergence | 收敛 |
| Comparison oracle | 比较预言机 |
| Normal direction estimation | 法线方向估计 |
| Line search | 线搜索 |
| Off-diagonal contribution | 非对角贡献 |

---

## Unit-by-Unit Q&A

### Unit 0: Big Picture
- **Q:** VMC 和经典 MC 的区别？
- **A:** VMC 有**双层结构**：
  - 内循环：固定 α，MC 采样构型 σ ~ |ψ(σ;α)|²，估计 E(α)——和经典 MC 完全一样
  - 外循环：比较不同 α 的能量，更新 α——这是 VMC 独有的优化层
  - 经典 MC 只有内循环（固定 T），没有外循环
- FFO-CO 的角色：替换外循环中的优化器

### Unit 1: Variational Principle & Simple Ansatz
- **变分原理**：E(α) = ⟨ψ(α)|H|ψ(α)⟩ / ⟨ψ(α)|ψ(α)⟩ ≥ E₀
  - 任意试探波函数能量不低于真实基态 → min E(α) 逼近 E₀
- **Ansatz 选择**（表达能力 vs 优化难度 trade-off）：
  - Jastrow: ψ(σ) = exp(Σ αᵢⱼ σᵢσⱼ)——两体经典关联，参数 ~N，物理含义清晰
    - 物理图像：类似经典 Ising Boltzmann 权重，αᵢⱼ 扮演"有效耦合"
  - Slater-Jastrow: det[ϕⱼ(rᵢ)] × exp(Σ αᵢⱼ σᵢσⱼ)——单粒子量子 + 两体关联
  - NQS: 神经网络(σ;θ)——万能逼近，参数 10⁴-10⁶，物理意义不明显
- **为什么 Jastrow 适合 FFO-CO 起步**：参数少（~N），landscape 不太复杂，有精确解可对照
- **Q:** Jastrow 是自由平面波吗？
- **A:** 不是。Jastrow = exp(Σ αᵢⱼ σᵢσⱼ) 是两体关联因子，物理上对应经典 Ising Boltzmann 权重。αᵢⱼ > 0 → 铁磁关联，αᵢⱼ < 0 → 反铁磁关联。

### Unit 2: MC Energy Estimation
- **为什么需要 MC？** 变分原理本身不需要 MC。MC 是因为 Hilbert 空间指数大（2ᴺ），求和 E(α) = Σ_σ |ψ(σ)|² E_loc(σ) / Z 有 2ᴺ 项不能穷举，但每项可逐点计算 → MC 采样
- **E_loc(σ;α)** = ⟨σ|H|ψ(α)⟩ / ψ(σ;α)
  - 对角部分：H|σ⟩ ∝ |σ⟩ → 经典能量
  - 非对角部分：H|σ⟩ 包含自旋翻转 → 需要 ψ(σ')/ψ(σ) 的比值
  - "局域"含义：E_loc 是在单个构型 σ 处评估的能量，包含波函数在 σ 附近邻域的信息
- **采样流程**：
  1. 提议 σ'（翻转一个自旋）
  2. 接受概率 A = min(1, |ψ(σ')|²/|ψ(σ)|²) = exp(2·log_ratio)
  3. log_ratio = Σ αᵢⱼ(σ'ᵢσ'ⱼ − σᵢσⱼ) 只涉及被翻转的自旋 → O(1) 计算
  4. Ė(α) = (1/Ns) Σ E_loc(σₖ)，统计误差 σ_E = std(E_loc)/√Ns
- **噪声来源**：有限样本 Ns → σ_E ~ σ(E_loc)/√Ns

### Unit 3: Optimization Methods
- **SGD**: Δα = −η∇E
  - 有限差分梯度：∂E/∂αᵢ ≈ (E(α+h·eᵢ)−E(α))/h → 噪声 2σ_E²/h²
  - 解析梯度：∂E/∂αᵢ = 2(⟨E_loc·Oᵢ⟩−⟨E_loc⟩⟨Oᵢ⟩)，Oᵢ = ∂ln ψ/∂αᵢ → 两个估计量协方差，噪声大
- **SR (Stochastic Reconfiguration)**: Δα = −η·S⁻¹·∇E
  - Sᵢⱼ = Cov(Oᵢ, Oⱼ) = Fisher information matrix / QGT
  - 物理含义：natural gradient，考虑参数空间曲率——平坦方向走大步，陡峭方向走小步
  - 问题：S⁻¹ 在高维不稳定（小特征值放大噪声），计算 O(p³)
  - 近似：对角近似 O(p)、CG 求解、正则化 S+λI
- **Linear Method (LM)**: 二阶展开 → 解线性方程组，精确但 O(p³)
- **Adam/RMSProp**: 自适应学习率，不利用量子结构
- **KFAC**: Kronecker 乘积近似 Fisher → 适合神经网络大参数量
- **方法对比**：

| 方法 | 信息量 | 计算成本 | 收敛质量 | 适合参数量 |
|------|--------|---------|---------|-----------|
| SGD | 一阶 | 低 | 差 | 任意 |
| Adam | 一阶+动量 | 低 | 中 | 任意 |
| SR | 二阶近似 | O(p²)−O(p³) | 好 | 中等(~10³) |
| LM | 完全二阶 | O(p³) | 很好 | 小 |
| KFAC | 近似二阶 | O(p) | 好 | 大(~10⁶) |

- **FFO-CO 定位**：不在上表中——不用 ∇E，只用 sign(ΔE)

### Unit 4: Noise Analysis — FFO-CO 的机会
- **噪声三层结构**：

| 量 | 噪声量级 | 信噪比 |
|---|---|---|
| Ė(α) | σ_E | — |
| ∇̂E (有限差分) | √2·σ_E/h | \|∂E/∂α\|·h/(√2·σ_E) |
| ∇̂E (解析) | ~10-100×σ_E | 取决于协方差质量 |
| sign(ΔE) | P(correct) = Φ(\|ΔE\|/(√2·σ_E)) | \|ΔE\|/σ_E 可通过大 h 增强 |

- **有限差分的噪声灾难**：Var(∇̂E) = 2σ_E²/h²，h↓ → 噪声↑↑
- **sign 的鲁棒性**：|ΔE|/σ_E > 3 → P(correct) > 97%
- **关键区别**：
  - sign 的信号 |ΔE| 可以通过增大比较半径 h 来增强（FFO-CO 自带）
  - 梯度的信号 |∂E/∂αᵢ| 是固定的，估计它又需要除以 h
  - 增大 Ns 对两者都有效，但 sign 可以"便宜地"增强信号，梯度不能
- **核心论点**："信息更多 ≠ 做得更好"——当梯度大小不可靠时，扔掉大小只用 sign，配合自适应 h，可能更稳健

### Unit 5: Implementation — 轻量 VMC + FFO-CO
- **为什么自写而不用 NetKet**：FFO-CO 不是梯度方法，optax 接口不自然；框架重 debug 不透明
- **代码架构**（~500 行 Python）：

```
vmc_ffoco/
├── ansatz.py      # Jastrow: log_psi, log_ratio, deriv
├── sampler.py     # Metropolis 采样（标准，复用 lecture 1）
├── optimizer.py   # FFO-CO / SGD / SR 优化器
└── main.py        # 主循环
```

- **关键函数**：
  - `log_ratio(σ, σ', i)` → O(1) 计算（只涉及被翻转自旋的 α 行）
  - `estimate_energy(ansatz, H, n_samples)` → MC 估计 E(α) ± σ_E
  - `compare(α₁, α₂)` → sign(ΔE) + 信噪比 |ΔE|/σ
  - `FFOCO.estimate_normal(α, energy_func)` → 正交补法线方向估计
- **FFO-CO 替换的是 optimizer**：其余模块（ansatz, sampler, E_loc）完全复用
- **测试体系**：

| 体系 | N | 参数 | 精确解 |
|------|---|------|--------|
| 1D Ising L=10 | 10 | 9 | Onsager |
| 1D Heisenberg L=10 | 10 | 9 | Bethe ansatz |
| 1D Heisenberg L=20 | 20 | 19 | Bethe ansatz |
| 2D Heisenberg 4×4 | 16 | 24 | ED |

- **实现路线图**：
  1. Jastrow + Metropolis + energy estimation (~100 行)
  2. SGD + 1D Heisenberg 测试 (~50 行)
  3. SR + 与 SGD 对比 (~80 行)
  4. FFO-CO 优化器 (~150 行)
  5. Benchmark 三者对比 + 可视化 (~100 行)

---

## Open-Source Packages Reference

| Package | Language | Custom Optimizer | Built-in SR | Domain | 推荐度 |
|---------|----------|-----------------|-------------|--------|--------|
| NetKet | Python/JAX | optax API | SR, QGT | 格点多体 | 可用但偏重 |
| FermiNet | Python/JAX | 需改内部代码 | KFAC | 量子化学 | 不适合 |
| DeepQMC | Python/JAX | Protocol 接口 | KFAC | 量子化学 | 不适合 |
| QuCumber | Python/PyTorch | PyTorch 优化器 | 无 | 态层析 | 不适合 |

**推荐方案**：自写轻量 VMC + FFO-CO，用 NumPy/JAX，不依赖 ML 框架。Benchmark 对照：SGD、SR 手写或借用 NetKet 的 SR 实现。

---

## FFO-CO in VMC Pipeline Summary

```
ansatz (Jastrow) → sampler (Metropolis) → estimate E(α) ± σ_E
                                              ↓
                              optimizer: SGD / SR / **FFO-CO**  ← 只换这里
                                              ↓
                                    update α → repeat
```

FFO-CO 输入：compare(α₁, α₂) → sign(ΔE)，不依赖 ∇E 或 S 矩阵。
其余模块（波函数、采样、E_loc 计算）完全复用。
