# Fock-State Generation via Repeated Post-Selected Measurements (MATLAB)

This repository contains a minimal MATLAB simulation of **Fock-state generation in a single-mode cavity** using **repeated projection (post-selected) measurements** on an auxiliary two-level system (qubit). The dynamics are modeled by a **Lindblad master equation** for an interacting qubit–cavity system, followed by a **measurement update rule** after each interaction round.

The code demonstrates how repeated *successful* measurement outcomes can progressively increase the overlap (fidelity) of the cavity state with a desired target Fock state (or a chosen Fock-state superposition).

---

## Files

- `one_mode.m`  
  Main script: sets parameters/operators, initializes the state, runs repeated evolution + projection cycles, and plots:
  - target overlap vs round number
  - final cavity photon-number distribution

- `Lindbladcav.m`  
  Right-hand-side function for `ode45`, implementing the Lindblad master equation:
  \[
  \dot{\rho} = -i[H,\rho] + \sum_j \mathcal{D}[L_j]\rho
  \]
  with dissipators for qubit relaxation, cavity loss, and qubit dephasing (configurable rates).

---

## Physical Model (as implemented)

### Hilbert space
- Qubit: dimension 2  
- Cavity: truncated Fock basis dimension `nn`  
- Total dimension: \(2 \times nn\)

### Hamiltonian (Jaynes–Cummings-type coupling)
In the script:
- Bare terms are defined (but currently weighted by `0*H0`), and the interaction term is enabled:
  - `Hi = g * kron([0,1;0,0], a)`
  - `H = Hi + Hi'`

### Dissipation channels (Lindblad form)
The ODE function `Lindbladcav.m` includes:
- Qubit relaxation with rate `gamma`
- Cavity loss with rate `kappa`
- Qubit dephasing with rate `gammap`

(See `Lindbladcav.m` for the exact superoperator expression.)

### Measurement / post-selection step (each round)
After evolving for an interaction time `tau`, the code applies a projector
\[
\rho \rightarrow \frac{(M_e \otimes I)\rho(M_e \otimes I)}{\mathrm{Tr}[(M_e \otimes I)\rho]}
\]
i.e. it **post-selects** the qubit being found in the state associated with `Me`.

The per-round success probability is accumulated in `Pg`.

---

## Requirements

- MATLAB (tested with MATLAB R2018a-style syntax)
- No external toolboxes required beyond standard ODE solvers (`ode45`)

---

## Quick Start

1. Clone/download the repository.
2. Open MATLAB in the repo folder.
3. Run:
   ```matlab
   main
