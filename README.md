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
  **rotation-gauge sector (Phase 0a), position-gauge sector (Phase 0b), a
  narrowly-scoped curvature quadratic invariant (Phase 0c), a general-h
  reciprocal frame (Phase 0d), and a GA-native Minkowski-paired reciprocal
  frame (Phase 0e) ONLY**. See "Gauge Theory Gravity scope" below before
  reading or extending this namespace -- it is a narrow, literature-faithful
  port of two sectors, one deliberately limited curvature invariant, and two
  dual-basis constructions, not a gravity engine.

## Gauge Theory Gravity scope (`kotoba.sm.gtg`, Phase 0a + Phase 0b + Phase 0c + Phase 0d + Phase 0e)

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
   A second, independent bug was later found (and fixed) in
   `kotoba.sm.gauge/self-interaction-term` itself: it indexed the
   structure-constant array with the output/free generator index in the
   wrong array slot -- numerically invisible for the compact `U(1)`/`SU(2)`/
   `SU(3)` groups (whose structure constants are totally antisymmetric under
   any index permutation), but it broke the physically required
   `R_mu-nu = -R_nu-mu` antisymmetry for this `so(1,3)` generator set
   whenever the self-interaction term's two excited `Omega` components
   mixed a boost-type and a rotation-type index (this generator set's
   trace-Gram matrix is diagonal but **not** uniform, +1 for rotation-type
   vs -1 for boost-type -- see point 3). **This has since been fixed**
   (`kotoba.sm.gauge/self-interaction-term`'s docstring has the derivation);
   `R_mu-nu` is now antisymmetric for every excitation pattern, including
   mixed boost/rotation ones -- see `gtg_test.cljc`'s
   `rotation-field-strength-self-interaction-now-antisymmetric-*` tests for
   the before/after numeric record.
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

**Phase 0d (general-position-gauge reciprocal frame, a deliberately
two-staged addition) adds:**

11. **Stage 1 (implemented): a general `NxN` matrix inverse and a reciprocal
    (dual) frame for a general invertible `h`.** `kotoba.sm.tensor/mat-det`,
    `mat-minor`, `mat-adjugate`, `mat-inverse` -- ordinary Laplace-expansion/
    adjugate linear algebra (textbook, no physics convention involved),
    throwing `ex-info` on a `~0` determinant rather than dividing by it, the
    same spirit as `kotoba.sm.gauge/diagonal-gram-real`'s guard on a `~0`
    Gram-matrix diagonal entry -- and `kotoba.sm.gtg/reciprocal-frame`: the
    reciprocal (dual) frame `h-bar^mu` of a GENERAL invertible
    position-gauge-field VALUE `h` (not just `position-gauge-identity` or a
    constant Lorentz matrix, Phase 0b's scope), defined by the standard
    linear-algebra dual-basis biorthogonality relation
    `h-bar^mu . h_nu = delta^mu_nu`, where `.` is the PLAIN component
    (Kronecker) pairing, NOT `kotoba.sm.tensor/dot`'s Minkowski-metric
    pairing. This gives `Hbar = (H^-1)^T` (`mat-inverse` then
    `mat-transpose`) -- verified: `h * h^-1 = h^-1 * h = I` for a concrete
    non-identity invertible `h`; the biorthogonality relation holds
    numerically for that same `h`; a singular (`~0` determinant) `h` throws
    `ex-info`; and at `h = position-gauge-identity`, `h-bar = h` EXACTLY
    (bit-for-bit -- the identity matrix is its own inverse and its own
    transpose).

12. **Stage 2 (investigated, deliberately NOT implemented): the true
    Lasenby-Doran-Gull curvature scalar `R`.** Per the task instructions
    for this phase, stage 2 was to be implemented ONLY if the exact
    index/convention of the LDG formula
    `R = h^mu ^ h^nu . R(h-bar_nu, h-bar_mu)` could be confirmed from at
    least two independent literature sources that agree. The investigation
    found two independent, mutually-corroborating PRIMARY sources (by the
    original author team) that DO agree on the Ricci scalar's high-level
    structure:
    - Lasenby, Doran & Gull, "Gravity, Gauge Theories and Geometric
      Algebra", Phil. Trans. R. Soc. Lond. A (1998) 356 487-582, updated
      version [arXiv:gr-qc/0405033](https://arxiv.org/abs/gr-qc/0405033),
      section 4 "The field equations", eqns 4.9, 4.11, 4.12.
    - Lewis, Doran & Lasenby, "Quadratic Lagrangians and Topology in Gauge
      Theory Gravity", [arXiv:gr-qc/9910039](https://arxiv.org/abs/gr-qc/9910039),
      section 2, eqns 5, 6, 12, 14.

    Both sources give (translating into this namespace's notation)
    `R = sum_{a,b} gamma^a . (gamma^b . R(h_b ^ h_a))`, where `h_b = h(e_b)`
    is exactly Phase 0b's existing `h`, `gamma^a` is the reciprocal of the
    FIXED background orthonormal frame (trivial under the Minkowski metric,
    `kotoba.sm.tensor/raise` applied to the ordinary basis vectors -- **not**
    a reciprocal frame of `h`), and `R(x^y)` is the bilinear extension of
    Phase 0a's `rotation-field-strength` to an arbitrary bivector argument
    (rigorously justified from `R_mu-nu`'s own defining bilinearity in its
    two vector slots, needing no further literature confirmation). This is
    a genuine, useful finding on its own: it shows Stage 1's
    `reciprocal-frame` (a reciprocal frame *of* `h`) is **not** what this
    formula's outer contraction needs at all, correcting the phase's
    original working assumption. **However**, translating
    `gamma^b . R(h_b^h_a)` (a vector-dot-bivector Clifford-algebra
    contraction) into this codebase's six-generator-component
    `R[mu][nu][k]` array representation requires an ADDITIONAL, nontrivial
    step -- which spacetime/generator-index slot the contracting vector's
    index pairs with, and with what sign/factor -- that neither source
    spells out in a form directly transcribable to this array
    representation. A hand re-derivation via the standard GA identity
    `v.(a^b)=(v.a)b-(v.b)a` produced a candidate formula, but with
    unresolved factor-of-2/sign risk once expressed over the six generator
    components, and no independent worked numeric example was found to test
    a candidate implementation against. This is exactly the class of
    silent-index-bug risk `kotoba.sm.gauge/self-interaction-term`'s
    since-fixed slot-order bug (see `rotation-field-strength`'s docstring)
    and `curvature-quadratic-invariant`'s own honesty note (point 10 above)
    already flag for this precise problem. Per the task's own instructions,
    **declining to implement here is the correct outcome, not a failure**:
    `curvature-scalar` remains unimplemented, and this writeup (plus the
    matching docstring on `kotoba.sm.gtg/reciprocal-frame`) records the
    investigation so a future pass does not have to repeat the literature
    search from scratch.

**Phase 0e (GA-native Minkowski-paired reciprocal frame, plus a bounded,
declined stage-2 curvature-scalar investigation) adds:**

13. **`kotoba.sm.gtg/reciprocal-frame-minkowski`: the GA-native reciprocal
    frame, using the Minkowski pairing instead of Stage 1's plain pairing.**
    Stage 1's `reciprocal-frame` (point 11 above) uses the plain
    component/Kronecker pairing; its own HONESTY NOTE already derives, but
    stops short of implementing, the alternative built from
    `kotoba.sm.tensor/dot`'s Minkowski invariant pairing instead -- the
    pairing Geometric Algebra / GTG actually use for a scalar product.
    `reciprocal-frame-minkowski` implements and independently re-verifies
    exactly that: `Hbar_eta = (H^-1)^T . eta` (re-derived from the
    biorthogonality condition `sum_a eta[a][a] Hbar[mu][a] H[nu][a] =
    delta[mu][nu]`, i.e. the matrix equation `Hbar . eta . H^T = I`).
    **Verified** (`gtg_test.cljc`): (a) the Minkowski-paired biorthogonality
    relation holds for a concrete non-identity invertible `h`, and the
    result differs numerically from `reciprocal-frame`'s plain-pairing
    output for the *same* `h` -- concretely, row 0 is `[0.5 0 0.5 0]`
    (Minkowski-paired) versus `[0.5 0 -0.5 0]` (plain-paired), agreeing on
    the timelike column and flipping sign on every spatial column exactly as
    predicted by `Hbar_eta = Hbar_flat . eta` (right-multiplying by the
    diagonal `eta` scales column `a` by `eta[a][a]`). (b) at
    `h = position-gauge-identity`, `h-bar_eta = eta` **exactly**
    (bit-for-bit): `h-bar^0 = h_0` but `h-bar^i = -h_i` for the 3 spatial
    rows -- genuinely different from `reciprocal-frame`'s `h-bar = h` at the
    same point, confirming the divergence `reciprocal-frame`'s HONESTY NOTE
    predicted without numerically checking it. (c) for a constant Lorentz
    transformation `h = Lambda` (`kotoba.sm.vector-field`'s existing
    boost/rotation matrices), `h-bar_eta = eta . Lambda` **exactly** (up to
    floating-point roundoff) -- a closed form derived from `vf/lorentz?`'s
    defining property `Lambda^T eta Lambda = eta`, checked for several
    concrete boosts/rotations.

14. **A bounded, declined attempt at independently verifying the true LDG
    curvature scalar `R` via linearized weak-field General Relativity
    (per this phase's task instructions: independent numerical
    verification, *not* a further literature read).** Plan: (i) pick a small
    metric perturbation `h_mu-nu(x)` with the known linearized-GR Ricci
    scalar `R^(1) = d^mu d^nu h_mu-nu - []h` (a standard textbook result);
    (ii) build a position-gauge field `h_mu(x) = delta_mu^nu +
    (1/2)h_mu^nu(x)` reproducing it via `derived-metric` to linear order
    (worked out and confirmed straightforward: `derived-metric` of that `h`
    equals `eta_mu-nu + h_mu-nu(x) + O(h^2)` exactly as expected, no new
    machinery needed); (iii) obtain the *corresponding* rotation-gauge field
    `Omega_mu(x)` from `h_mu(x)`'s derivatives; (iv) run
    `rotation-field-strength`; (v) test a candidate `R[mu][nu][k]->scalar`
    formula against `R^(1)`. **This attempt stopped at step (iii).** This
    codebase has no `h_mu(x) -> Omega_mu(x)` correspondence anywhere
    (grep-confirmed), and unlike the ordinary general-relativity
    vierbein-postulate torsion-free spin connection (which relates a
    curved-spacetime tetrad `e^a_mu` -- one flat index, one curved-coordinate
    index -- to Christoffel symbols), GTG's `h_mu^nu` has *both* indices as
    components of the *same* flat background vector space: `h` maps a flat
    gauge-coordinate direction to a flat physical-spacetime direction, with
    no curved-coordinate side to import the ordinary tetrad-postulate
    formula onto. The genuine GTG-specific "intrinsic" Omega-from-h relation
    (part of LDG's own gauge-invariance/field-equation construction) is
    exactly the kind of GTG-paper-specific, index/sign-convention-bearing
    formula this pass was told not to re-derive from the literature this
    round, and it cannot be safely reconstructed from generic (non-GTG)
    general-relativity textbook knowledge -- ordinary tetrad calculus does
    not carry over to GTG's flat, non-curved-coordinate formalism without
    exactly the kind of paper-specific translation step this codebase's own
    established practice already treats as too risky to guess at without an
    independent reference value (same class of risk as point 12's blocker,
    `rotation-field-strength`'s and `curvature-quadratic-invariant`'s
    docstrings). So this investigation stopped at the *same class* of
    blocker Phase 0d stage 2 hit -- an unverifiable index/sign convention,
    no independent numeric reference value -- one step *earlier* in the
    pipeline (`h -> Omega`, rather than `R -> scalar`). **Per the task's own
    instructions, declining to implement here is the correct outcome, not a
    failure**: `curvature-scalar` remains unimplemented.

**Explicitly out of scope, not implemented here** (deliberately, deferred to
a later phase if ever pursued): the true Lasenby-Doran-Gull curvature scalar
`R` for **any** `h` (see point 12 above for the Phase 0d investigation and
exactly which remaining step blocked implementation, point 10 above for how
Phase 0c's quadratic invariant differs from it even at
`h = position-gauge-identity`, and point 14 above for Phase 0e's further,
declined attempt via independent numerical verification and exactly where
*that* attempt stopped), the Einstein multivector, the GTG action principle
and field equations, any proof of equivalence to General Relativity, and any
dark-matter/dark-energy/de-Sitter extension. **Also out of scope**: the FULL
combined `h_mu` + `Omega_mu` GTG covariant derivative (using `h` to relate
the rotation-gauge connection to an actual spacetime-vector
derivative/torsion) -- Phase 0a's covariant derivative (5) and Phase 0b's
derived metric (8) are developed independently of each other here, not yet
combined, and Phase 0c's curvature invariant (10) is built from Phase 0a's
`R_mu-nu` alone, without reference to `h` at all (neither Phase 0d's
`reciprocal-frame` nor Phase 0e's `reciprocal-frame-minkowski` changes this --
neither is wired into either, and point 14 above records exactly why Phase 0e
could not bridge this gap either); that combination is deferred to a later
phase (Phase 0f or beyond). **This namespace implements the rotation-gauge
sector's generator algebra/covariant derivative, the position-gauge sector's
derived metric, one narrowly-scoped curvature quadratic invariant, and two
general-h dual-basis reciprocal frames (plain-pairing and Minkowski-pairing)
only** -- it should not be described as "GTG implemented" or as any kind of
gravity engine.

## Develop

```sh
clojure -M:test
```
