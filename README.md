# MPS-AVP

This repository builds on the MPS, DMRG, and TDVP directory developed by [Nishan Ranabhat](https://github.com/NishanRanabhat/Quench_dyn). It retains the history and source of [Quench_dyn](https://github.com/naraguho/Quench_dyn); its original README is preserved in [QuenchDyn_README.md](QuenchDyn_README.md). The new work here explores the adiabatic variational principle (AVP) as a route to **real-time dissipative dynamics of kinks** in the easy-axis XXZ chain. The AVP framework is introduced in:

> Yang Yang, Ali Rayat, Ho Jang, Yunhao Fan, Sankha Subhra Bakshi, and Gia-Wei Chern, **“Adiabatic Variational Principle: Born-Oppenheimer Theory for Emergent Electronic-Order Dynamics,”** [arXiv:2610.40117](https://arxiv.org/abs/2610.40117) (2026).

The real-space examples prepare localized one- and two-kink MPS states at zero longitudinal field. A separate constrained-DMRG example records an exploratory AVP target evolution, including stages where the kink profile may disappear. These codes document numerical procedures; they do not establish validated dissipative kink dynamics.

## Real-space kink examples

See [kinks-in-real-space/README.md](kinks-in-real-space/README.md) for the physical construction, the exact optimization steps, convergence checks, limitations, and output files.

The [constrained DMRG one-kink pipeline](constrained_DMRG/README.md) contains the local Lagrange-multiplier solver and target evolution used to study how an unrestricted kink can lose its profile as the field changes.

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --threads=1 --project=. kinks-in-real-space/one_kink.jl
julia --startup-file=no --threads=1 --project=. kinks-in-real-space/two_kinks.jl
```

Both examples are single Julia files. They load the Quench_dyn primitives from this repository's `src/` directory; neither requires a checkpoint, another example script, nor a site-specific filesystem path. They save their MPS, energy history, and spin profiles under `kinks-in-real-space/outputs/`. If Python `matplotlib` is available, they also save a one-panel PNG.

## Model and scope

We use an open spin-½ chain with $J=1$:

```math
H_{\mathrm{XXZ}}
=
\sum_{i=1}^{N-1}
\Bigl[
S_i^{x}\,S_{i+1}^{x}
\; + \;
S_i^{y}\,S_{i+1}^{y}
\; + \;
\Delta\,S_i^{z}\,S_{i+1}^{z}
\Bigr],
\qquad
N=100,
\qquad
\Delta=3.
```

The examples set all longitudinal fields to zero, including during preparation of the Néel background. A product Néel seed chooses one broken-symmetry branch. Local optimization around the defects then lowers the zero-field XXZ energy within a restricted MPS family.
