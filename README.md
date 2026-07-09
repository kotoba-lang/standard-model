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
- **`kotoba.sm.gtg`** -- Gauge Theory Gravity (Lasenby-Doran-Gull 1998)
  **rotation-gauge sector (Phase 0a), position-gauge sector (Phase 0b), and a
  narrowly-scoped curvature quadratic invariant (Phase 0c) ONLY**. See "Gauge
  Theory Gravity scope" below before reading or extending this namespace --
  it is a narrow, literature-faithful port of two sectors plus one deliberately
  limited curvature invariant, not a gravity engine.

## Gauge Theory Gravity scope (`kotoba.sm.gtg`, Phase 0a + Phase 0b + Phase 0c)

Lasenby-Doran-Gull's 1998 observation is that the restricted Lorentz group
`SO(1,3)+` can be gauged with the exact same machinery `kotoba.sm.gauge`
already implements generically for the internal `U(1)/SU(2)/SU(3)` gauge
groups: the six bivector generators `T^{ab} = (i/4)[gamma^a,gamma^b]` of
`so(1,3)`, built from `kotoba.sm.spinor`'s Dirac gamma matrices, plug
straight into `kotoba.sm.gauge`'s generic structure-constant /
covariant-derivative / field-strength functions and reproduce the
rotation-gauge curvature `R_mu-nu` as an instance of that same generic
machinery -- no separate Yang-Mills-for-gravity implementation was written.

**Implemented (established literature content ported to this codebase, plus
new code that exercises it):**

1. the six `so(1,3)` bivector generators `T^{ab}`, a<b, as 4x4 complex
   matrices, reusing `kotoba.sm.spinor`'s existing gamma matrices and
   `sigma-munu`, and `kotoba.sm.complex`'s matrix algebra.
2. a numeric check that these generators close the Lorentz Lie algebra
   `[T^{ab},T^{cd}] = i(eta^{bc}T^{ad} - eta^{ac}T^{bd} - eta^{bd}T^{ac} +
   eta^{ad}T^{bc})`.
3. an **honesty check, not an assumption**: whether `kotoba.sm.gauge`'s
   (Phase-0a-era) compact-group trace normalization `Tr(T^aT^b)=1/2delta^ab`
   (documented there at the time as the basis for its `structure-constants`
   derivation, and verified in that namespace's tests only for the compact
   `SU(2)`/`SU(3)`) carries over to these noncompact `so(1,3)` generators. **It
   did not.** Numerically (`kotoba.sm.gtg/generator-trace-gram`,
   `kotoba.sm.gtg/compact-group-trace-normalization-holds?`, and
   `gtg_test.cljc`): the Gram matrix `Tr(T_A T_B)` is diagonal (same structure
   as the compact case) but its diagonal entries are **+/-1, not a uniform
   +1/2** -- exactly -1 for the 3 boost-type generators and +1 for the 3
   rotation-type generators, i.e. the indefinite Minkowski signature leaks
   directly into the generators' own trace pairing. Applying
   `kotoba.sm.gauge/structure-constants` (Phase-0a version) to this generator
   set unmodified therefore did **not** give the genuine `so(1,3)` structure
   constants: the raw output was off by a factor `2*Tr(T_C T_C)` from the true
   value, i.e. **the wrong magnitude (2x) for every triple, and additionally
   the wrong SIGN whenever the third generator was boost-type**. **This has
   since been fixed**: `kotoba.sm.gauge/structure-constants` was generalized
   to compute its own generator set's trace-Gram matrix `K[A][B]=Tr(T_AT_B)`
   instead of assuming a uniform `1/2delta^AB` (identical output for
   `SU(2)`/`SU(3)`, since their Gram matrix genuinely is `1/2delta^AB` -- see
   `gauge_test.cljc`'s `structure-constants-matches-legacy-*` regression
   tests), so applying it to this `so(1,3)` generator set now gives the
   genuine structure constants, matching `kotoba.sm.gtg/true-structure-constants`
   (`gtg_test.cljc`'s `structure-constants-now-match-genuine-values-after-gauge-fix`).
   A further, related noncompactness symptom recorded in the same test file
   (unaffected by the above fix): the 3 boost-type generators are
   anti-Hermitian (not Hermitian), while the 3 rotation-type generators are
   Hermitian -- the textbook reason (Peskin & Schroeder section 3.2) the Dirac
   spinor representation of the Lorentz group is not unitary. The
   `generator-trace-gram`/`true-structure-constants`/
   `compact-group-trace-normalization-holds?` functions are kept as an
   independent derivation and a live record of why the compact-group
   assumption fails here, even though `kotoba.sm.gauge/structure-constants`
   no longer needs a workaround.
4. the rotation gauge field `Omega_mu` (a 6-component real bivector-valued
   gauge potential, in the same data shape `kotoba.sm.gauge` already uses)
   and its field strength / curvature bivector `R_mu-nu`, via
   `kotoba.sm.gauge/field-strength` applied to this generator set unmodified.
   The pure curl term `d_mu Omega_nu - d_nu Omega_mu` carries no
   normalization caveat (it does not involve structure constants at all); the
   self-interaction term uses the now-fixed `structure-constants` from point 3.
5. the spin-connection covariant derivative on a Dirac spinor,
   `D_mu psi = d_mu psi - i g Omega_mu^k T_k psi`, via
   `kotoba.sm.gauge/covariant-derivative` applied to this generator set
   unmodified.
6. a flat-limit check: with `Omega_mu = 0` everywhere, the covariant Dirac
   residual built on (5) is exactly `kotoba.sm.spinor`'s existing free
   Dirac-equation residual, verified numerically identical (not merely
   close), not just structurally argued.

**Phase 0b (position-gauge sector) adds:**

7. the position gauge field `h_mu(x)`: at each spacetime point, an
   invertible linear map from gauge coordinates to physical spacetime,
   represented as a 4x4 real matrix `h[mu][nu]` = `h_mu^nu` (row `mu` IS the
   physical four-vector `h_mu`, in the same `[t x y z]` component convention
   `kotoba.sm.tensor`/`kotoba.sm.vector-field` use everywhere else).
   `kotoba.sm.gtg/position-gauge-identity` is the trivial choice
   `h_mu^nu = delta_mu^nu` ("gauge coordinates = physical coordinates"). As
   with `kotoba.sm.gauge`'s own documented omission of a general matrix
   inverse, this namespace does not implement a determinant/invertibility
   check for `h` -- invertibility is a precondition on `h`, not something
   computed here.
8. the derived metric `g_mu-nu(x) = h_mu(x) . h_nu(x) = eta_ab h_mu^a(x)
   h_nu^b(x)` (`kotoba.sm.gtg/derived-metric`, `derived-metric-field`),
   computed by treating each row of `h` as an upper-index four-vector and
   reusing `kotoba.sm.tensor/dot` (the existing Minkowski inner product) row
   by row -- no new inner-product machinery. **Verified EXACTLY** (bit-for-bit
   integer arithmetic, not merely close within a tolerance) to reduce to
   `kotoba.sm.tensor/metric` at `h = position-gauge-identity`
   (`gtg_test.cljc`'s `flat-limit-h-reproduces-minkowski-metric-exactly`).
9. a consistency check that `kotoba.sm.vector-field`'s pre-existing global
   `SO(3,1)+` representation (`boost`/`rotation-x`/`rotation-y`/`rotation-z`)
   is contained in this sector as the special case "`h_mu` = a constant
   Lorentz matrix, `Omega_mu = 0`": a constant `h` built from any Lorentz
   transformation reproduces the flat metric exactly (up to floating-point
   roundoff from `sqrt`/`cos`/`sin`), because `Lambda^T eta Lambda = eta`
   (`vector-field/lorentz?`'s defining property) algebraically implies
   `Lambda eta Lambda^T = eta`, which is exactly what `derived-metric`
   computes row-by-row for `h = Lambda` (derivation in `gtg.cljc`'s Phase-0b
   section-7 header comment; numeric check in `gtg_test.cljc`'s
   `global-lorentz-transformations-are-the-constant-flat-special-case`, for
   several concrete boosts and rotations).

**Phase 0c (curvature quadratic invariant, deliberately narrow) adds:**

10. `kotoba.sm.gtg/curvature-quadratic-invariant`: a scalar built from the
    rotation-gauge curvature bivector `R_mu-nu` **alone**, by lowering both
    spacetime indices with `kotoba.sm.tensor/lower2` and fully contracting
    with `kotoba.sm.tensor/full-contract` (reusing `kotoba.sm.tensor`'s
    **existing** rank-2 raise/lower/contract machinery -- no new tensor
    primitives added), per bivector-generator component `k`, summed over the
    6 components weighted by this namespace's own diagonal Killing-form-like
    trace pairing `generator-trace-gram` (point 3 above):
    `I = sum_k K[k][k] * (R_mu-nu^k R^{mu-nu}_k)` -- the
    `eta^mu-rho eta^nu-sigma R_mu-nu . R_rho-sigma` quadratic contraction, the
    so(1,3) analogue of `F_mu-nu^a F^mu-nu_a`.
    **This is honestly NOT the Lasenby-Doran-Gull curvature scalar
    `R = h^mu ^ h^nu . R(h-bar_nu, h-bar_mu)`.** That definition needs the
    position gauge field's *reciprocal frame* `h-bar`, i.e. a general 4x4
    matrix inverse -- not implemented anywhere in this codebase, and general-`h`
    curvature is explicitly out of scope (see below). The one case where
    `h-bar` is trivial (`h-bar = h` = the identity, no inverse needed) is
    `h = position-gauge-identity`, so this function is only *interpretable* in
    that special case, and it does not take an `h` argument at all (`R_mu-nu`
    is computed independently of `h`, see the still-deferred "full combined
    covariant derivative" note below). **Even at
    `h = position-gauge-identity`**, working through the true definition shows
    it is a **single LINEAR** contraction of `R_mu-nu` against the fixed
    bivector `h^mu^h^nu` -- the same relationship the (linear) Ricci scalar
    has to the full Riemann tensor -- while `curvature-quadratic-invariant`
    computes a **QUADRATIC** (`R` contracted against itself) invariant
    instead, the same relationship the Kretschmann scalar `R_abcd R^abcd` has
    to the Riemann tensor in ordinary GR. **Ricci scalar and Kretschmann
    scalar are different invariants of different degree in the curvature**:
    related only in that both vanish exactly when the full curvature vanishes
    (`gtg_test.cljc`'s `curvature-quadratic-invariant-vanishes-in-the-flat-limit`),
    not proportional to each other and not interchangeable in general -- so
    this function's output must not be read as "the GTG curvature scalar" or
    used as a stand-in for it. The true linear LDG contraction is deliberately
    **not** implemented either: pinning down which spacetime index pairs with
    which bivector-value index needs an LDG-convention detail this pass has no
    independent literature value to test against, so it is left as a named
    follow-up rather than guessed at. `curvature-quadratic-invariant` is also
    **not guaranteed non-negative** (`generator-trace-gram`'s `K` is
    indefinite, `+1`/`-1` -- a boost-dominated curvature can make it negative,
    another noncompactness symptom). Verified: exactly `0.0` in the flat limit
    (`Omega_mu = 0` everywhere), and a hand-checked nonzero value (`10.0`) for
    a concrete curvature bivector with a single rotation-type component,
    matching a by-hand `tensor/lower2` + `tensor/full-contract` computation
    (`gtg_test.cljc`'s `curvature-quadratic-invariant-nonzero-for-concrete-curvature`).

**Explicitly out of scope, not implemented here** (deliberately, deferred to
a later phase if ever pursued): the true Lasenby-Doran-Gull curvature scalar
`R` for **any** `h` (a general 4x4 matrix inverse / reciprocal frame `h-bar`
is not implemented anywhere in this codebase -- see point 10 above for
exactly how Phase 0c's quadratic invariant differs from it, even at
`h = position-gauge-identity`), the Einstein multivector, the GTG action
principle and field equations, any proof of equivalence to General
Relativity, and any dark-matter/dark-energy/de-Sitter extension. **Also out
of scope**: the FULL combined `h_mu` + `Omega_mu` GTG covariant derivative
(using `h` to relate the rotation-gauge connection to an actual
spacetime-vector derivative/torsion) -- Phase 0a's covariant derivative (5)
and Phase 0b's derived metric (8) are developed independently of each other
here, not yet combined, and Phase 0c's curvature invariant (10) is built from
Phase 0a's `R_mu-nu` alone, without reference to `h` at all; that combination
is deferred to a later phase (Phase 0d or beyond). **This namespace
implements the rotation-gauge sector's generator algebra/covariant
derivative, the position-gauge sector's derived metric, and one narrowly-
scoped curvature quadratic invariant only** -- it should not be described as
"GTG implemented" or as any kind of gravity engine.

## Develop

```sh
clojure -M:test
```
