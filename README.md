# kotoba-lang/standard-model

[![CI](https://github.com/kotoba-lang/standard-model/actions/workflows/ci.yml/badge.svg)](https://github.com/kotoba-lang/standard-model/actions/workflows/ci.yml)

Classical Standard Model field algebra in zero-dep, portable `.cljc`: **gauge
fields**, **tensor fields**, **vector fields**, and **spinor fields**, and the
equations of motion that relate them. See
[`90-docs/adr/2607051500-kotoba-lang-standard-model-field-theory-library.md`](../../../90-docs/adr/2607051500-kotoba-lang-standard-model-field-theory-library.md)
in the superproject for the full design rationale and scope boundary.

**Scope: classical field content, evaluated numerically.** This is *not* a
quantum field theory engine -- no path-integral quantization, no
renormalization/running couplings, no scattering amplitudes/cross-sections, no
lattice simulation. Those are legitimate, much larger follow-up efforts and are
out of scope here.

## Namespaces

- **`kotoba.sm.complex`** -- complex scalar (`[re im]`) and complex-matrix
  algebra (add/mul/dagger/commutator/trace). The shared foundation `spinor`
  and `gauge` build their fixed-size matrices on top of.
- **`kotoba.sm.tensor`** -- Minkowski tensor algebra. Metric fixed to
  **mostly-minus** `eta = diag(1,-1,-1,-1)` (Bjorken-Drell convention),
  index raise/lower, rank-2 contraction, the Levi-Civita symbol.
- **`kotoba.sm.vector-field`** -- four-vectors, Lorentz boosts/rotations
  (the restricted Lorentz group SO(3,1)+), and a finite-difference
  four-gradient/four-divergence/d'Alembertian for numeric field functions.
- **`kotoba.sm.spinor`** -- Pauli matrices, Dirac gamma matrices (Dirac
  representation), chirality projectors, bilinear covariants, the
  free-particle plane-wave solution, and a Dirac-equation residual check.
- **`kotoba.sm.gauge`** -- generic compact-Lie-algebra layer: U(1)/SU(2)/SU(3)
  generators, structure constants derived generically (not hard-coded) via
  `f^abc = -2i Tr([T^a,T^b] T^c)`, the gauge covariant derivative, and the
  non-abelian Yang-Mills field-strength tensor.
- **`kotoba.sm.standard-model`** -- composition: the fermion generation table
  (checked against Gell-Mann-Nishijima `Q = T3 + Y/2`), electroweak symmetry
  breaking (Weinberg angle, photon/Z/W±, `M_W = M_Z cos(theta_W)`), the Higgs
  potential/vev/Yukawa mass generation, the CKM quark-mixing matrix, the full
  `SU(3)_c x SU(2)_L x U(1)_Y` covariant derivative, and Lagrangian-density
  term assembly (Yang-Mills + Dirac + Higgs + Yukawa).

## Develop

```sh
clojure -M:test
```
