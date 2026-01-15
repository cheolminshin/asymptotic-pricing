# Asymptotic pricing (short-maturity, near-the-money) — Numerical experiments

This repository contains the c++ code used to generate the numerical simulation results reported in:

**“[Asymptotic pricing of short-maturity near-the-money options in stochastic volatility models](https://doi.org/10.1016/j.frl.2025.108922)”**
*Financial Research Letters, 2026*
(see **Table 1 (p.9)** and **Table 2 (p.10)** in the paper)

In particular, the experiments compare:
- **Baseline Monte Carlo (MC)** pricing, and
- A **proposed / new numerical method** introduced in the paper

for **Asian option** prices (and related quantities), focusing on the **short-maturity / near-the-money** regime under **stochastic volatility models(Heston / SABR)**.

> Goal: Provide a reproducible implementation that allows readers to regenerate Table 1 and Table 2, and to validate the agreement / differences between the baseline MC estimator and the proposed method.

---

## Usage

This repository provides a single entry-point program for reproducing the numerical results in the paper.

1. Download (or clone) this repository.
2. Compile and run `simulation_3.cpp`.
3. The program generates the following output files:
   - `output_table1(heston).csv` (Table 1: Heston model)
   - `output_table2(SABR).csv` (Table 2: SABR model)

### Compile & run (example)

```bash
g++ -std=c++17 simulation_3.cpp -o simulation_3
./simulation_3
