(ns kotoba.sm.gtg
  "Gauge Theory Gravity (Lasenby-Doran-Gull 1998, 'Gravity, gauge theories and
  geometric algebra') -- ROTATION-GAUGE SECTOR (Phase 0a), POSITION-GAUGE
  SECTOR (Phase 0b), A NARROWLY-SCOPED CURVATURE QUADRATIC INVARIANT
  (Phase 0c, h = position-gauge-identity ONLY -- NOT the LDG curvature
  scalar), AND A GENERAL-h RECIPROCAL FRAME (Phase 0d) ONLY. See the 'Scope'
  section of the repo README for the ADR reference.

  GTG's central observation (well known in the literature, not a new result
  of this namespace) is that the restricted Lorentz group SO(1,3)+ can be
  gauged exactly like the internal Yang-Mills groups already implemented in
  `kotoba.sm.gauge`: the six bivector generators T^{ab} of so(1,3), built
  from `kotoba.sm.spinor`'s existing Dirac gamma matrices, plug directly into
  `kotoba.sm.gauge`'s *generic* (not hard-coded to U(1)/SU(2)/SU(3))
  structure-constant / covariant-derivative / field-strength functions and
  reproduce the rotation-gauge-field curvature R_mu-nu -- the GTG analogue of
  the Yang-Mills field strength -- as a direct instance of that same generic
  machinery.

  SCOPE -- this namespace implements ONLY:
    1. the six so(1,3) bivector generators T^{ab} = (i/4)[gamma^a,gamma^b],
       a<b, as 4x4 complex matrices (`generators`, `generator`).
    2. a numeric check (in gtg_test.cljc) that these six generators close the
       Lorentz Lie algebra: [T^{ab},T^{cd}] = i(eta^{bc}T^{ad} - eta^{ac}T^{bd}
       - eta^{bd}T^{ac} + eta^{ad}T^{bc}) (`lorentz-algebra-rhs`).
    3. an HONEST numeric check of whether `kotoba.sm.gauge`'s (as of Phase 0a)
       compact-group trace normalization Tr(T^aT^b) = 1/2 delta^ab carries
       over to this noncompact generator set. IT DID NOT -- see
       `generator-trace-gram`, `true-structure-constants`,
       `compact-group-trace-normalization-holds?` and gtg_test.cljc for the
       actual numbers and the failure mode this Phase 0a check exposed in
       `kotoba.sm.gauge/structure-constants`'s output (a factor of 2*K[C][C]
       off, wrong SIGN for boost-type indices). THIS HAS SINCE BEEN FIXED:
       `kotoba.sm.gauge/structure-constants` was generalized (post-Phase-0a)
       to compute its own generator set's trace-Gram matrix K[A][B]=Tr(T_AT_B)
       instead of assuming it is 1/2delta^AB, so it is no longer restricted to
       the compact-group normalization and its output for this so(1,3)
       generator set is now the GENUINE structure constants -- see
       `rotation-raw-structure-constants` and gtg_test.cljc's
       `structure-constants-now-match-genuine-values-after-gauge-fix` (the
       renamed/inverted former `raw-structure-constants-diverge-from-genuine-ones`).
       This section's functions (`generator-trace-gram`, `true-structure-constants`,
       `compact-group-trace-normalization-holds?`, `basis-projection`) are KEPT as
       an independent derivation of the genuine structure constants and a live,
       re-checkable record of WHY Tr(T^aT^b)=1/2delta^ab fails to hold for a
       noncompact generator set -- useful documentation and cross-check value
       even though `kotoba.sm.gauge/structure-constants` no longer needs the
       workaround.
       A related, independently-checked noncompactness symptom (gtg_test.cljc,
       'generators-are-4x4-complex-matrices'): the 3 boost-type generators
       T^{0i} are ANTI-Hermitian, not Hermitian, while the 3 rotation-type
       generators T^{ij} are Hermitian -- the textbook reason (Peskin &
       Schroeder section 3.2, whose S^{mu-nu}=(i/4)[gamma^mu,gamma^nu] is
       literally this namespace's T^{ab}) the Dirac spinor representation of
       the Lorentz group is not unitary.
    4. the rotation gauge field Omega_mu (a 6-component real bivector-valued
       gauge potential, same data shape `kotoba.sm.gauge` already uses for
       A_mu^a) and its field strength / curvature bivector R_mu-nu, obtained
       by handing this generator set to `kotoba.sm.gauge/field-strength`
       unmodified (`rotation-field-strength`). R_mu-nu^k = -R_nu-mu^k
       (the physically required antisymmetry) now holds for EVERY excited
       Omega component, including mixed boost/rotation self-interaction --
       `kotoba.sm.gauge/self-interaction-term` had an independent,
       since-fixed structure-constant slot-order bug (invisible for compact
       U(1)/SU(2)/SU(3) groups, but visible here because this generator
       set's trace-Gram matrix is diagonal yet non-uniform) that used to
       break this antisymmetry for mixed boost/rotation excitations; see
       `rotation-field-strength`'s docstring and gtg_test.cljc's
       `rotation-field-strength-self-interaction-now-antisymmetric-*` tests
       for the before/after record.
    5. the spin-connection covariant derivative on a Dirac spinor,
       D_mu psi = d_mu psi - i g Omega_mu^{ab} T_{ab} psi, obtained by
       handing this generator set to `kotoba.sm.gauge/covariant-derivative`
       unmodified (`spinor-covariant-derivative`).
    6. a flat-limit check (in gtg_test.cljc) that with Omega_mu = 0
       everywhere, the covariant Dirac residual built on top of (5) is
       EXACTLY (not just approximately) `kotoba.sm.spinor`'s existing free
       Dirac-equation residual (`dirac-residual-with-rotation-gauge`).

  PHASE 0b (position-gauge sector) adds:
    7. the position gauge field h_mu(x): at each spacetime point, an
       invertible linear map from gauge coordinates to physical spacetime,
       represented as a 4x4 real matrix `h[mu][nu]` = h_mu^nu (row `mu` IS
       the physical four-vector h_mu, in the same [t x y z] component
       convention `kotoba.sm.tensor`/`kotoba.sm.vector-field` use everywhere
       else). `position-gauge-identity` is the trivial choice h_mu^nu =
       delta_mu^nu ('gauge coordinates = physical coordinates'). AT PHASE 0b
       TIME this namespace did not implement a general matrix inverse to
       check h's invertibility numerically or to build a reciprocal frame for
       general h (the same deliberate omission `kotoba.sm.gauge` documents
       for its own structure-constant derivation) -- invertibility was a
       precondition on `h`, not something computed. PHASE 0d (item 11 below)
       adds a general NxN matrix inverse (`kotoba.sm.tensor/mat-inverse`) and
       uses it to build a reciprocal frame `reciprocal-frame` for ANY
       invertible h, not just `position-gauge-identity` -- invertibility is
       now CHECKED (throws `ex-info` on a ~0 determinant) rather than merely
       assumed.
    8. the derived metric g_mu-nu(x) = h_mu(x) . h_nu(x) = eta_ab h_mu^a(x)
       h_nu^b(x) (`derived-metric`, `derived-metric-field`), computed by
       treating each row of `h` as an upper-index four-vector and reusing
       `kotoba.sm.tensor/dot` (the existing Minkowski inner product) row by
       row -- no new inner-product machinery. A flat-limit check (in
       gtg_test.cljc) that at h = `position-gauge-identity` this is EXACTLY
       (bit-for-bit, integer arithmetic throughout, not merely close)
       `kotoba.sm.tensor/metric`.
    9. a consistency check that `kotoba.sm.vector-field`'s PRE-EXISTING
       global SO(3,1)+ representation (`boost`/`rotation-x`/`rotation-y`/
       `rotation-z`) is CONTAINED in this sector as the special case
       'h_mu = a constant Lorentz matrix, Omega_mu = 0' -- a constant h built
       from any Lorentz transformation reproduces the flat metric exactly
       (up to floating-point roundoff from sqrt/cos/sin), because
       Lambda^T eta Lambda = eta (`vector-field/lorentz?`'s defining
       property) algebraically implies Lambda eta Lambda^T = eta, which is
       exactly what `derived-metric` computes row-by-row for h = Lambda --
       see the section-7 header comment below for the derivation, and
       gtg_test.cljc's
       `global-lorentz-transformations-are-the-constant-flat-special-case`
       for the numeric check (several concrete boosts and rotations).
       `derived-metric-matches-flat?` packages this equality check (used for
       both (8)'s exact identity case and (9)'s Lorentz-matrix case).

  PHASE 0c (curvature quadratic invariant, DELIBERATELY NARROW SCOPE) adds:
    10. `curvature-quadratic-invariant`: a scalar built from the rotation-gauge
       curvature bivector R_mu-nu ALONE (`rotation-field-strength`'s output),
       by lowering both spacetime indices with `kotoba.sm.tensor/lower2` and
       fully contracting with `kotoba.sm.tensor/full-contract` (reusing
       `kotoba.sm.tensor`'s EXISTING rank-2 raise/lower/contract machinery,
       no new tensor primitives), per bivector-generator component k, then
       summing the 6 components weighted by this namespace's own diagonal
       Killing-form-like trace pairing `generator-trace-gram` (section 3
       above) -- i.e. I = sum_k K[k][k] (R_mu-nu^k R^{mu-nu}_k), the
       'eta^mu-rho eta^nu-sigma R_mu-nu . R_rho-sigma' quadratic contraction
       collapsed via the diagonal metric.
       THIS IS DELIBERATELY NOT the Lasenby-Doran-Gull curvature scalar
       R = h^mu ^ h^nu . R(h-bar_nu, h-bar_mu). That definition needs the
       position gauge field's RECIPROCAL FRAME h-bar, i.e. a general 4x4
       matrix inverse, which is not implemented ANYWHERE in this codebase
       (see item 7's note and `kotoba.sm.gauge`'s own analogous documented
       omission) -- computing R for a general h is explicitly OUT OF SCOPE
       and deferred (see below). The one case where h-bar IS trivial to get
       right is h = `position-gauge-identity` (h-bar = h = the identity), so
       `curvature-quadratic-invariant` is scoped to being INTERPRETABLE only
       in that special case; it does not take an `h` argument at all because
       R_mu-nu itself (Phase 0a) is computed independently of h (Phase 0b) --
       see the FULL COMBINED covariant derivative note below, still deferred.
       EVEN AT h = `position-gauge-identity`, `curvature-quadratic-invariant`
       is honestly NOT a port of the LDG R: working through the definition,
       R = h^mu ^ h^nu . R(h-bar_nu, h-bar_mu) is a single LINEAR contraction
       of R_mu-nu (dot the FIXED bivector h^mu^h^nu, which does not depend on
       R at all, against R -- bilinear in two DIFFERENT bivectors, one of
       which is a constant), the same relationship the ordinary linear Ricci
       scalar R has to the full Riemann tensor. `curvature-quadratic-invariant`
       is instead QUADRATIC in R_mu-nu (R contracted against itself) -- the
       same relationship the (quadratic) Kretschmann scalar R_abcd R^abcd has
       to the Riemann tensor in ordinary GR. Ricci scalar and Kretschmann
       scalar are DIFFERENT invariants of DIFFERENT degree in the curvature:
       related only in that both vanish exactly when the full curvature
       vanishes (gtg_test.cljc's flat-limit check for this one), not
       proportional to each other and not interchangeable in general. This
       namespace does not implement the true linear LDG contraction either --
       doing so correctly requires pinning down, from the LDG geometric-
       algebra convention (not guessed at here), exactly which of R_mu-nu's
       two spacetime indices pairs with which of the curvature bivector's own
       two indices, and there is no independent literature reference value in
       this codebase to test that pairing against -- so it is left as a
       named, explicit follow-up rather than risked as a plausible-looking
       but silently wrong implementation. `curvature-quadratic-invariant` is
       also NOT guaranteed non-negative: `generator-trace-gram`'s K is
       indefinite (+1 rotation-type, -1 boost-type), so a boost-dominated
       curvature can make it negative -- a further symptom, like
       `compact-group-trace-normalization-holds?`'s finding, of this
       generator set's noncompactness.

  PHASE 0d (general-position-gauge reciprocal frame, DELIBERATELY STAGED --
  see the README's Phase 0d section for the two-stage split) adds:
    11. a GENERAL NxN matrix inverse (`kotoba.sm.tensor/mat-det`,
       `mat-minor`, `mat-adjugate`, `mat-inverse` -- ordinary Laplace-
       expansion/adjugate linear algebra, no physics convention involved,
       throws `ex-info` on a ~0 determinant rather than dividing by it, same
       spirit as `kotoba.sm.gauge/diagonal-gram-real`'s guard) and
       `reciprocal-frame`: the GTG reciprocal (dual) frame h-bar^mu of a
       GENERAL invertible position-gauge-field VALUE h (not just
       `position-gauge-identity` or a constant Lorentz matrix, Phase 0b's
       scope), defined by the STANDARD linear-algebra dual-basis
       biorthogonality relation h-bar^mu . h_nu = delta^mu_nu, where '.' here
       is the PLAIN component (Kronecker) pairing sum_a h-bar^mu[a] h_nu[a]
       -- NOT `kotoba.sm.tensor/dot`'s Minkowski-metric pairing. Given this
       pairing, biorthogonality is a matrix equation Hbar . H^T = I, so
       Hbar = (H^-1)^T -- `reciprocal-frame` computes exactly this, via
       `mat-inverse` then `kotoba.sm.tensor/mat-transpose`; gtg_test.cljc
       checks the biorthogonality relation holds numerically for a concrete
       non-identity invertible h, that h*h^-1=I for that same h, that a
       singular (determinant-0) h throws `ex-info`, and that at
       h = `position-gauge-identity`, h-bar = h EXACTLY (the identity matrix
       is its own inverse and its own transpose).

       *** HONESTY NOTE ***: this is the ORDINARY linear-algebra dual basis
       of the ROWS of h, treating them as an abstract basis of R^4 under the
       PLAIN component pairing -- it is deliberately NOT asserted to be an
       ingredient of the Lasenby-Doran-Gull curvature scalar. A literature
       follow-up (stage 2 of this phase, see `reciprocal-frame`'s own
       docstring and the README's Phase 0d section for the full writeup and
       both source URLs) found TWO independent, mutually-corroborating
       primary sources (Lasenby/Doran/Gull 1998/2004, arXiv gr-qc/0405033;
       Lewis/Doran/Lasenby 1999, arXiv gr-qc/9910039) pinning down the
       ACTUAL Ricci-scalar formula's high-level structure -- and it turns
       out that formula's outer contraction uses the FIXED BACKGROUND
       frame's OWN reciprocal (trivial under the metric, no matrix inverse
       of h needed at all), NOT `reciprocal-frame` (a reciprocal frame OF h)
       as originally guessed. So `reciprocal-frame` is confirmed NOT to be
       what R needs, rather than left merely unconfirmed either way. R
       itself remains unimplemented because a DIFFERENT step -- translating
       the confirmed formula's 'vector . bivector' Clifford-algebra
       contraction into this codebase's six-generator-component array
       representation -- has its own unresolved index/sign risk that this
       pass could not independently verify (see `reciprocal-frame`'s
       docstring for the details). `reciprocal-frame` remains useful in its
       own right as h's standard dual basis (e.g. for a future combined
       covariant derivative, still out of scope), just not as a confirmed
       LDG-R ingredient.

  EXPLICITLY OUT OF SCOPE, NOT IMPLEMENTED HERE (Phase 0a+0b+0c+0d): the true
  Lasenby-Doran-Gull curvature scalar R (see item 10's note for exactly how
  Phase 0c's quadratic invariant differs from it, and item 11's honesty note
  for the Phase 0d literature investigation's findings and exactly which
  remaining step blocked implementation);
  the Einstein multivector, the GTG action principle and field equations, any
  proof of equivalence to General Relativity, and any dark-matter/dark-
  energy/de-Sitter extension. ALSO OUT OF SCOPE: the FULL combined
  h_mu+Omega_mu GTG covariant derivative (i.e. using h to relate the
  rotation-gauge connection to an actual spacetime-vector derivative/torsion)
  -- this namespace's Phase-0a covariant derivative (5) and Phase-0b derived
  metric (8) are developed independently of each other, not yet combined, and
  Phase 0c's curvature invariant (10) is built from Phase 0a's R_mu-nu alone,
  without reference to h at all (Phase 0d's `reciprocal-frame` does not
  change this -- it is not yet wired into either (5) or (10)). Those are
  legitimate, much larger follow-up efforts, deliberately deferred to a later
  phase (Phase 0e or beyond) -- this namespace is a limited,
  literature-faithful port of the rotation-gauge and position-gauge SECTORS,
  one narrowly-scoped curvature invariant, and a general-h dual-basis
  reciprocal frame ONLY, and must not be read as 'GTG is implemented here' or
  as any kind of gravity engine."
  (:require [kotoba.sm.complex :as c]
            [kotoba.sm.tensor :as tensor]
            [kotoba.sm.spinor :as spinor]
            [kotoba.sm.gauge :as gauge]))

(defn- abs-num [x] (if (neg? x) (- x) x))

;; ---------------------------------------------------------------------------
;; 1. so(1,3) bivector generators T^{ab} = (i/4)[gamma^a,gamma^b], a<b
;; ---------------------------------------------------------------------------

(def generator-index-pairs
  "The 6 independent (a,b), a<b index pairs, in the fixed enumeration order
  used throughout this namespace to index the 6-component rotation-gauge
  field / generator vector: 0=(0 1) 'boost-x', 1=(0 2) 'boost-y',
  2=(0 3) 'boost-z', 3=(1 2) 'rotation-z', 4=(1 3) 'rotation-y',
  5=(2 3) 'rotation-x' (the usual identification of the 3 boost + 3 rotation
  generators of the restricted Lorentz group)."
  [[0 1] [0 2] [0 3] [1 2] [1 3] [2 3]])

(defn bivector-generator
  "T^{ab} = (i/4)[gamma^a,gamma^b] for a<b, a 4x4 complex matrix. Built from
  `kotoba.sm.spinor`'s existing Dirac gamma matrices and
  `kotoba.sm.complex`'s matrix ops; equal to (1/2) `kotoba.sm.spinor/sigma-munu`
  (sigma^{ab} = (i/2)[gamma^a,gamma^b] is already defined there as 'the
  antisymmetric tensor bilinear generator' -- reused here rather than
  reimplementing the same commutator)."
  [a b]
  (c/m-rscale 0.5 (spinor/sigma-munu a b)))

(def generators
  "The 6 so(1,3) generators, in `generator-index-pairs` order. This is the
  generator SET this whole namespace hands to `kotoba.sm.gauge`'s generic
  functions in place of the SU(2)/SU(3) generator sets they were originally
  exercised against."
  (mapv (fn [[a b]] (bivector-generator a b)) generator-index-pairs))

(defn generator
  "T^{ab} for ANY a,b in 0..3 (not just a<b), via antisymmetry T^{ba}=-T^{ab}
  and T^{aa}=0 (the zero 4x4 matrix). Needed because the Lorentz-algebra
  commutation relation (`lorentz-algebra-rhs`) references T^{ad} etc. for
  index orderings outside the a<b enumeration of `generators`."
  [a b]
  (cond
    (= a b) (c/m-zero 4 4)
    (< a b) (bivector-generator a b)
    :else (c/m-rscale -1 (bivector-generator b a))))

;; ---------------------------------------------------------------------------
;; 2. Lorentz algebra closure: [T^{ab},T^{cd}] = i(eta^{bc}T^{ad} - eta^{ac}T^{bd}
;;    - eta^{bd}T^{ac} + eta^{ad}T^{bc})
;; ---------------------------------------------------------------------------

(defn- eta [i j] (get-in tensor/metric [i j]))

(defn lorentz-algebra-rhs
  "The right-hand side of the defining so(1,3) commutation relation
  [T^{ab},T^{cd}] = i(eta^{bc}T^{ad} - eta^{ac}T^{bd} - eta^{bd}T^{ac} +
  eta^{ad}T^{bc}), evaluated with this namespace's generators and
  `kotoba.sm.tensor`'s mostly-minus metric. gtg_test.cljc checks this against
  the actual commutator `kotoba.sm.complex/m-commutator` of the corresponding
  generators, for every a<b, c<d pair."
  [a b c d]
  (let [term1 (c/m-rscale (eta b c) (generator a d))
        term2 (c/m-rscale (- (eta a c)) (generator b d))
        term3 (c/m-rscale (- (eta b d)) (generator a c))
        term4 (c/m-rscale (eta a d) (generator b c))
        total (reduce c/m-add [term1 term2 term3 term4])]
    (c/m-scale c/i total)))

;; ---------------------------------------------------------------------------
;; 3. Trace-normalization honesty check -- Tr(T^aT^b)=1/2delta^ab (the
;;    compact-group normalization `kotoba.sm.gauge/structure-constants` used
;;    to hard-code) does NOT hold for these noncompact so(1,3) generators.
;;    Phase 0a recorded this as a live, re-checkable predicate/value (not
;;    just asserted in the test file), and it exposed a real bug in
;;    `kotoba.sm.gauge/structure-constants` (see `rotation-raw-structure-constants`).
;;    That bug has SINCE BEEN FIXED: `kotoba.sm.gauge/structure-constants` now
;;    computes its OWN generator set's trace-Gram matrix instead of assuming
;;    it is 1/2delta^AB, so it is correct for this so(1,3) generator set too.
;;    The functions below are KEPT (not deleted) as an independent derivation
;;    of the genuine structure constants and documentation of WHY the old
;;    compact-group assumption fails here -- useful on their own merits, and
;;    now doubling as the reference oracle that proves the gauge.cljc fix is
;;    correct (gtg_test.cljc's
;;    `structure-constants-now-match-genuine-values-after-gauge-fix`).
;; ---------------------------------------------------------------------------

(defn generator-trace-gram
  "K[A][B] = Tr(T_A T_B) for the 6 basis generators (`generator-index-pairs`
  order) -- the so(1,3) analogue of the compact-group normalization
  Tr(T^aT^b)=1/2delta^ab that `kotoba.sm.gauge/structure-constants` used to
  hard-code (it now computes this same kind of Gram matrix generically from
  whatever generator set it is given, `kotoba.sm.gauge/generator-gram`).
  Returned as a plain 6x6 real matrix; gtg_test.cljc additionally checks the
  imaginary part of every trace is ~0 before this real part is taken at face
  value."
  []
  (vec (for [A (range 6)]
         (vec (for [B (range 6)]
                (c/re (c/m-trace (c/m-mul (nth generators A) (nth generators B)))))))))

(defn compact-group-trace-normalization-holds?
  "Does Tr(T^A T^B) = 1/2 delta^AB hold for these so(1,3) generators?
  Numerically: NO. `generator-trace-gram` is diagonal (off-diagonal terms
  vanish, as for the compact groups) but the diagonal entries are +/-1, not a
  uniform +1/2: -1 for the 3 boost-type generators (index 0,1,2 = pairs
  containing the timelike index 0) and +1 for the 3 rotation-type generators
  (index 3,4,5 = purely spatial pairs) -- the indefinite Minkowski signature
  shows up directly in the generators' own Killing-form-like trace pairing.
  This predicate is exposed live (not just asserted once in the test file) as
  a permanent record of why this generator set needed
  `kotoba.sm.gauge/structure-constants`'s Gram-matrix generalization (it is
  diagonal, so that generalization applies -- see the `kotoba.sm.gauge`
  namespace docstring's SCOPE note on non-diagonal Gram matrices)."
  ([] (compact-group-trace-normalization-holds? 1e-9))
  ([eps]
   (let [K (generator-trace-gram)]
     (every? true?
             (for [A (range 6) B (range 6)]
               (< (abs-num (- (get-in K [A B]) (if (= A B) 0.5 0.0))) eps))))))

(defn- basis-projection
  "Coefficient of generator C in a bivector-space matrix X = sum_D coeff_D T_D
  -- valid because `generator-trace-gram` is DIAGONAL in this basis (checked
  by gtg_test.cljc), so coeff_C = Tr(X T_C) / K[C][C] with no cross terms."
  [X C K]
  (/ (c/re (c/m-trace (c/m-mul X (nth generators C))))
     (get-in K [C C])))

(defn true-structure-constants
  "The GENUINE so(1,3) structure constants in this basis, f-true[A][B][C],
  defined by [T_A,T_B] = i f-true[A][B][C] T_C (sum over C) -- obtained by
  taking the actual commutator of the basis generators and projecting it back
  onto the basis via `basis-projection` (using the REAL, non-1/2 trace
  pairing `generator-trace-gram`), independently of
  `kotoba.sm.gauge/structure-constants`'s implementation. Now that
  `kotoba.sm.gauge/structure-constants` has been generalized to compute its
  own Gram matrix, this independently-derived value and
  `rotation-raw-structure-constants` MATCH (gtg_test.cljc's
  `structure-constants-now-match-genuine-values-after-gauge-fix`) -- this
  function is kept anyway as an independent cross-check/derivation, not
  merely a historical artifact."
  []
  (let [K (generator-trace-gram)]
    (vec (for [A (range 6)]
           (vec (for [B (range 6)]
                  (let [[a b] (nth generator-index-pairs A)
                        [c d] (nth generator-index-pairs B)
                        comm (c/m-commutator (generator a b) (generator c d))
                        ;; comm = i * (sum_C f[C] T_C)  =>  (sum_C f[C] T_C) = -i * comm
                        y (c/m-scale (c/c 0 -1) comm)]
                    (vec (for [C (range 6)]
                           (basis-projection y C K))))))))))

(defn rotation-raw-structure-constants
  "`kotoba.sm.gauge/structure-constants` applied to this namespace's 6 so(1,3)
  generators. Named so `rotation-field-strength` and the test suite share one
  computed value instead of recomputing it ad hoc; kept its Phase-0a name
  ('raw') for continuity even though it is no longer 'raw' in the sense of
  being uncorrected -- see FIXED note below.

  FIXED (previously a CAVEAT -- see `compact-group-trace-normalization-holds?`
  / `true-structure-constants` / gtg_test.cljc): Phase 0a found that
  `kotoba.sm.gauge/structure-constants` assumed Tr(T^aT^b)=1/2delta^ab, which
  does NOT hold for this generator set, making its output for this call off
  from the genuine so(1,3) structure constants by a factor of 2*K[C][C] --
  i.e. +2 for a rotation-type third index and -2 (wrong SIGN, not just wrong
  magnitude) for a boost-type third index. `kotoba.sm.gauge/structure-constants`
  has SINCE BEEN GENERALIZED to compute its own generator set's trace-Gram
  matrix K[A][B]=Tr(T_AT_B) (which for THIS generator set is diagonal, so the
  generalization applies cleanly -- see `kotoba.sm.gauge`'s namespace
  docstring) instead of assuming it is 1/2delta^AB, so this call now returns
  the GENUINE so(1,3) structure constants, numerically matching
  `true-structure-constants` (verified in gtg_test.cljc)."
  []
  (gauge/structure-constants generators))

;; ---------------------------------------------------------------------------
;; 4. rotation gauge field Omega_mu and its field strength (curvature bivector)
;;    R_mu-nu, via kotoba.sm.gauge/field-strength applied to `generators`.
;; ---------------------------------------------------------------------------

(defn rotation-gauge-field-gradient
  "d Omega_nu^k / d x^mu via central finite differences, for a rotation-gauge
  field function Omega-field: R^4 -> (vector-of-4 vector-of-6 reals),
  Omega-field(x)[nu][k] = Omega_nu^k(x) (k indexes `generator-index-pairs`).
  A thin, GTG-named pass-through to `kotoba.sm.gauge/gauge-field-gradient` --
  this namespace does not reimplement finite differencing."
  ([Omega-field x] (gauge/gauge-field-gradient Omega-field x))
  ([Omega-field x h] (gauge/gauge-field-gradient Omega-field x h)))

(defn rotation-field-strength
  "R_mu-nu = d_mu Omega_nu - d_nu Omega_mu + g f^abc Omega_mu^b Omega_nu^c,
  the GTG curvature bivector, via `kotoba.sm.gauge/field-strength` applied to
  this namespace's 6 so(1,3) generators -- Lasenby-Doran-Gull 1998's
  observation that this is exactly the non-abelian Yang-Mills field-strength
  formula, unmodified, with the Lorentz-algebra bivector generators standing
  in for an internal gauge group.

  `d-Omega[mu][nu][k]` = d Omega_nu^k / d x^mu (see
  `rotation-gauge-field-gradient`); `Omega[mu][k]` = Omega_mu^k(x) at the same
  spacetime point. The CURL term d_mu Omega_nu - d_nu Omega_mu is pure
  antisymmetric differencing and never depended on the structure-constant
  normalization; the SELF-INTERACTION term g f^abc Omega_mu^b Omega_nu^c uses
  `rotation-raw-structure-constants`, which (see the FIXED note there) now
  returns the genuine so(1,3) structure constants.

  R_mu-nu^k = -R_nu-mu^k (the physically required curvature-bivector
  antisymmetry) NOW HOLDS EVEN WHEN THE EXCITED Omega COMPONENTS MIX
  BOOST-TYPE AND ROTATION-TYPE INDICES. This used to fail for mixed
  boost/rotation self-interaction: `kotoba.sm.gauge/self-interaction-term`
  read the structure-constant array with the output/free index `a` in the
  wrong array slot (a bug independent of, and found after, the Phase-0a
  trace-Gram-normalization fix noted above), which is numerically invisible
  for compact groups (U(1)/SU(2)/SU(3), whose structure constants are
  totally antisymmetric under any index permutation) but broke R's
  mu<->nu antisymmetry here whenever the self-interaction term's two
  summed-over Omega components had DIFFERENT `generator-trace-gram` signs
  (one boost K=-1, one rotation K=+1) -- see
  `kotoba.sm.gauge`'s `self-interaction-term` docstring for the fix and
  `gtg_test.cljc`'s
  `rotation-field-strength-self-interaction-now-antisymmetric-for-mixed-trace-gram-sign-components-after-index-order-fix`
  for the before/after numeric record."
  [d-Omega Omega g]
  (gauge/field-strength (rotation-raw-structure-constants) g d-Omega Omega))

;; ---------------------------------------------------------------------------
;; 5. spin-connection covariant derivative on a Dirac spinor
;; ---------------------------------------------------------------------------

(defn spinor-covariant-derivative
  "The spin-connection covariant derivative on a Dirac spinor,
  D_mu psi = d_mu psi - i g Omega_mu^{k} T_k psi (sum over the 6 generators,
  k indexing `generator-index-pairs`), via
  `kotoba.sm.gauge/covariant-derivative` applied to this namespace's 6 so(1,3)
  generators -- unmodified, the same generic function `kotoba.sm.standard-model`
  already reuses for the 3 internal SM gauge groups.

  `d-mu-psi` is the ordinary partial derivative of psi along this direction;
  `Omega-mu` is the 6-vector of rotation-gauge components Omega_mu^k(x) at
  this point; `g` is a coupling constant. At Omega-mu = [0 0 0 0 0 0] this
  reduces EXACTLY to `d-mu-psi` (checked in gtg_test.cljc against
  `kotoba.sm.spinor`'s free Dirac-equation residual, not just numerically
  close)."
  [d-mu-psi Omega-mu g psi]
  (gauge/covariant-derivative d-mu-psi generators Omega-mu g psi))

;; ---------------------------------------------------------------------------
;; 6. flat-limit check support: covariant Dirac-equation residual under the
;;    rotation gauge connection
;; ---------------------------------------------------------------------------

(defn dirac-residual-with-rotation-gauge
  "(i gamma^mu D_mu - m) psi(x), the covariant Dirac-equation residual under
  the rotation-gauge (spin) connection Omega_mu, evaluated numerically:
  the ordinary derivative via `kotoba.sm.spinor/spinor-four-gradient` (reused,
  not reimplemented), then `spinor-covariant-derivative`'s gauge correction on
  top. `Omega-mu-field` is a function R^4 -> (vector-of-4 vector-of-6 reals),
  same shape as `rotation-gauge-field-gradient`'s `Omega-field` argument.

  At Omega-mu-field = (constantly [[0 0 0 0 0 0] [0 0 0 0 0 0] [0 0 0 0 0 0]
  [0 0 0 0 0 0]]), this is EXACTLY (bit-identical, since
  `kotoba.sm.gauge/gauge-correction` returns the literal zero vector when
  every A_mu^a is 0 -- not merely small) `kotoba.sm.spinor/dirac-residual` on
  the same `psi-field`/`m`/`x` -- the flat-limit check in gtg_test.cljc."
  [psi-field Omega-mu-field m g x]
  (let [dpsi (spinor/spinor-four-gradient psi-field x)
        psi-x (psi-field x)
        Omega-x (Omega-mu-field x)
        D-psi (vec (for [mu (range 4)]
                      (spinor-covariant-derivative (nth dpsi mu) (nth Omega-x mu) g psi-x)))
        terms (for [mu (range 4)] (c/m-vec (nth spinor/gammas mu) (nth D-psi mu)))
        slash-D (reduce c/v-add (vec (repeat 4 c/zero)) terms)
        i-slash-D (c/v-scale c/i slash-D)
        m-psi (c/v-scale (c/c m 0) psi-x)]
    (c/v-sub i-slash-D m-psi)))

;; ---------------------------------------------------------------------------
;; 7+8+9. position gauge field h_mu(x) and the derived metric g_mu-nu(x),
;;    Phase 0b -- see the namespace docstring's Phase-0b SCOPE note (items
;;    7-9) for the full description.
;;
;;    COMPATIBILITY WITH THE EXISTING GLOBAL SO(3,1)+ REPRESENTATION (item 9):
;;    a CONSTANT position-gauge field h = Lambda, built from ANY
;;    `kotoba.sm.vector-field/boost`/`rotation-x`/`rotation-y`/`rotation-z`
;;    matrix (with Omega_mu = 0, i.e. no rotation-gauge connection -- Phase
;;    0a's sector switched off), is included in this sector as a special case
;;    that reproduces the flat metric EXACTLY, not approximately. Derivation:
;;    `vector-field/lorentz?` checks Lambda^T eta Lambda = eta. Right-multiply
;;    both sides by Lambda^-1: Lambda^T eta = eta Lambda^-1, then left-
;;    multiply by eta^-1 = eta (eta is its own inverse in this basis):
;;    Lambda^-1 = eta Lambda^T eta. Substitute into Lambda Lambda^-1 = I:
;;    Lambda eta Lambda^T eta = I, then right-multiply by eta^-1 = eta:
;;    Lambda eta Lambda^T = eta. So Lambda^T eta Lambda = eta ALSO implies the
;;    'other' contraction Lambda eta Lambda^T = eta -- and `derived-metric`
;;    applied to h = Lambda computes exactly
;;    (derived-metric Lambda)[mu][nu] = sum_a eta[a][a] Lambda[mu][a] Lambda[nu][a]
;;    which IS (Lambda eta Lambda^T)[mu][nu] (eta diagonal collapses the
;;    double sum to the single shared index a). Hence `derived-metric` of any
;;    Lorentz matrix equals `kotoba.sm.tensor/metric` exactly (up to
;;    floating-point roundoff from the sqrt/cos/sin inside `boost`/
;;    `rotation-*`), numerically checked in gtg_test.cljc's
;;    `global-lorentz-transformations-are-the-constant-flat-special-case`.
;; ---------------------------------------------------------------------------

(def position-gauge-identity
  "The trivial position gauge field h_mu^nu = delta_mu^nu (the 4x4 identity
  matrix): 'gauge coordinates = physical coordinates'. `h[mu][nu]` is
  h_mu^nu -- the physical-spacetime component `nu` of the image, under the
  position-gauge map h, of the gauge-coordinate-direction-`mu` basis vector,
  so row `mu` of this matrix IS the physical four-vector h_mu (same
  [t x y z]-component convention `kotoba.sm.tensor`/`kotoba.sm.vector-field`
  use everywhere else). Kept as plain integers (not floats) so that
  `derived-metric` applied to it stays exact integer arithmetic throughout,
  for a bit-for-bit (not merely approximate) flat-limit check against
  `kotoba.sm.tensor/metric` -- see gtg_test.cljc's
  `flat-limit-h-reproduces-minkowski-metric-exactly`."
  [[1 0 0 0]
   [0 1 0 0]
   [0 0 1 0]
   [0 0 0 1]])

(defn derived-metric
  "The GTG derived metric g_mu-nu = h_mu . h_nu = eta_ab h_mu^a h_nu^b, for a
  position-gauge-field VALUE `h` (a 4x4 real matrix, row `mu` = the physical
  four-vector h_mu -- see `position-gauge-identity`) AT A SINGLE spacetime
  point. Computed by treating each row of `h` as an upper-index four-vector
  and taking `kotoba.sm.tensor/dot` (the EXISTING Minkowski inner product,
  eta = diag(1,-1,-1,-1)) of every pair of rows -- this function supplies
  only the h_mu . h_nu bookkeeping on top of that existing inner product, no
  new linear-algebra machinery.

  Returned as a plain 4x4 real matrix, the same shape as
  `kotoba.sm.tensor/metric` itself -- which is exactly what this reduces to
  at `h` = `position-gauge-identity` (gtg_test.cljc's
  `flat-limit-h-reproduces-minkowski-metric-exactly`), and, more generally,
  at `h` = any constant Lorentz transformation matrix from
  `kotoba.sm.vector-field` (this namespace's Phase-0b section-7 header
  comment above, and gtg_test.cljc's
  `global-lorentz-transformations-are-the-constant-flat-special-case`)."
  [h]
  (vec (for [mu (range 4)]
         (vec (for [nu (range 4)]
                (tensor/dot (nth h mu) (nth h nu)))))))

(defn derived-metric-field
  "g_mu-nu(x) for a position-gauge FIELD function h-field: R^4 -> 4x4 matrix
  (h-field(x)[mu][nu] = h_mu^nu(x)), at spacetime point x -- `derived-metric`
  applied to (h-field x). The field-valued analogue of `derived-metric`, the
  same relationship `rotation-gauge-field-gradient`/`rotation-field-strength`
  already have to a single `Omega` VALUE via their `Omega-field` argument."
  [h-field x]
  (derived-metric (h-field x)))

(defn derived-metric-matches-flat?
  "Does the derived metric g_mu-nu built from position-gauge-field VALUE `h`
  equal `kotoba.sm.tensor/metric`, component by component, within `eps`?
  Used both for the exact flat-limit check (`h` = `position-gauge-identity`)
  and the global-Lorentz-as-special-case check (`h` = a constant Lorentz
  matrix from `kotoba.sm.vector-field`) -- see this namespace's Phase-0b
  section-7 header comment above for why both are expected to hold exactly,
  up to floating-point roundoff."
  ([h] (derived-metric-matches-flat? h 1e-9))
  ([h eps]
   (let [g (derived-metric h)]
     (every? true?
             (for [mu (range 4) nu (range 4)]
               (< (abs-num (- (get-in g [mu nu]) (get-in tensor/metric [mu nu]))) eps))))))

;; ---------------------------------------------------------------------------
;; 10. Phase 0c -- curvature quadratic invariant, narrowly scoped to
;;    h = position-gauge-identity. See the namespace docstring's Phase-0c
;;    SCOPE note (item 10) for the full derivation and, crucially, for why
;;    this is honestly NOT the Lasenby-Doran-Gull curvature scalar R.
;; ---------------------------------------------------------------------------

(defn curvature-bivector-component
  "R_mu-nu^k (`rotation-field-strength`'s output `R[mu][nu][k]`) with the
  bivector-generator index `k` held FIXED, reshaped into a plain 4x4 real
  rank-2 tensor R^{mu-nu} (upper spacetime indices) -- exactly the shape
  `kotoba.sm.tensor`'s raise/lower/contract functions (`lower2`,
  `full-contract`) already expect. Internal-use helper for
  `curvature-quadratic-invariant`, exposed because it is independently
  useful/checkable (e.g. against `rotation-field-strength-curl-term-matches-
  yang-mills-form`'s hand-built `R`)."
  [R k]
  (vec (for [mu (range 4)]
         (vec (for [nu (range 4)]
                (get-in R [mu nu k]))))))

(defn curvature-quadratic-invariant
  "A curvature INVARIANT built from the rotation-gauge curvature bivector
  R_mu-nu (`rotation-field-strength`'s output) ALONE:

    I = sum_k K[k][k] * (R_mu-nu^k R^{mu-nu}_k)

  where, for each of the 6 bivector-generator components `k`
  (`generator-index-pairs`), `R_mu-nu^k R^{mu-nu}_k` is
  `kotoba.sm.tensor/full-contract` of `curvature-bivector-component`'s
  fixed-k rank-2 slice against its own `kotoba.sm.tensor/lower2` (BOTH
  spacetime indices lowered via the Minkowski metric eta -- `tensor.cljc`'s
  EXISTING raise/lower/contract machinery, reused unmodified, no new
  tensor-algebra primitives added here), and `K` is `generator-trace-gram`
  (this namespace's own diagonal so(1,3) Killing-form-like trace pairing from
  section 3 above: +1 for the 3 rotation-type generators, -1 for the 3
  boost-type ones). This is exactly the
  'eta^mu-rho eta^nu-sigma R_mu-nu . R_rho-sigma' quadratic contraction (a
  bivector-valued Yang-Mills-field-strength-squared-style invariant, the
  so(1,3) analogue of F_mu-nu^a F^mu-nu_a), collapsed via the diagonal metric
  so only rho=mu, sigma=nu survive.

  *** HONESTY NOTE, READ BEFORE USING THIS AS 'the curvature scalar': *** this
  is NOT the Lasenby-Doran-Gull (1998) curvature scalar
  R = h^mu ^ h^nu . R(h-bar_nu, h-bar_mu). That true definition, even
  restricted to h = `position-gauge-identity` (the one case where the
  reciprocal frame h-bar = h = the identity trivially, needing no matrix
  inverse), is a SINGLE LINEAR contraction of R_mu-nu against the FIXED
  bivector h^mu^h^nu (which does not itself depend on R) -- the same
  relationship the (linear) Ricci scalar has to the full Riemann tensor. This
  function computes a QUADRATIC (R contracted against itself) invariant
  instead -- the same relationship the Kretschmann scalar R_abcd R^abcd has
  to the Riemann tensor in ordinary GR. Ricci scalar and Kretschmann scalar
  are DIFFERENT invariants of DIFFERENT degree in the curvature: related only
  in that both vanish exactly when the full curvature vanishes (see
  gtg_test.cljc's flat-limit test), NOT proportional to each other and NOT
  interchangeable in general -- so this function's numeric output should not
  be read as 'the GTG curvature scalar' or used as a stand-in for it. The
  true linear LDG contraction is deliberately NOT implemented here either
  (see the namespace docstring's Phase-0c SCOPE note, item 10, for why:
  pinning down which spacetime index pairs with which bivector-value index
  needs an LDG-convention detail this pass has no independent literature
  value to test against, so it is left as a named follow-up rather than
  guessed at).

  Also NOT guaranteed non-negative: `generator-trace-gram`'s K is indefinite
  (+1/-1), so a boost-dominated curvature can make this negative -- a further
  symptom, like `compact-group-trace-normalization-holds?`'s finding, of this
  generator set's noncompactness."
  [R]
  (let [K (generator-trace-gram)]
    (reduce +
            (for [k (range 6)]
              (let [Rk (curvature-bivector-component R k)
                    Rk-lower (tensor/lower2 Rk)]
                (* (get-in K [k k]) (tensor/full-contract Rk Rk-lower)))))))

;; ---------------------------------------------------------------------------
;; 11. Phase 0d -- reciprocal (dual) frame h-bar^mu for a GENERAL invertible
;;    position-gauge field h, via `kotoba.sm.tensor/mat-inverse` (a general
;;    NxN matrix inverse, ordinary linear algebra, no physics convention).
;;    See the namespace docstring's Phase-0d SCOPE note (item 11) for the
;;    full derivation and, crucially, for the honesty note on why this is a
;;    standard dual-basis construction and NOT asserted to be a literature-
;;    exact port of the Lasenby-Doran-Gull geometric-algebra reciprocal frame
;;    for a general h.
;; ---------------------------------------------------------------------------

(defn reciprocal-frame
  "The reciprocal (dual) frame h-bar^mu of a GENERAL invertible position-gauge
  field VALUE `h` (a 4x4 real matrix, row `mu` = the physical four-vector
  h_mu -- same layout as `position-gauge-identity`/`derived-metric`, see item
  7's docstring), defined by the STANDARD linear-algebra dual-basis
  biorthogonality relation

    h-bar^mu . h_nu = delta^mu_nu

  where '.' here is the PLAIN component (Kronecker) pairing
  sum_a h-bar^mu[a] h_nu[a] -- the ordinary dual-basis construction for ANY
  basis {h_mu} of R^4, independent of any metric and independent of any
  GTG-paper-specific convention (see this function's HONESTY NOTE below for
  exactly how this differs from, and is NOT verified against, the
  Lasenby-Doran-Gull geometric-algebra reciprocal frame, which uses the
  Minkowski-metric pairing `kotoba.sm.tensor/dot` instead of this plain
  component pairing).

  DERIVATION: write `h` as a 4x4 matrix H (row `mu` = h_mu) and the sought
  reciprocal frame as a 4x4 matrix Hbar (row `mu` = h-bar^mu). The
  biorthogonality condition, component by component, is
  sum_a Hbar[mu][a] * H[nu][a] = delta[mu][nu]; the left side is exactly
  (Hbar . H^T)[mu][nu] (ordinary matrix product), so the condition is the
  matrix equation Hbar . H^T = I, i.e. Hbar = (H^T)^-1 = (H^-1)^T. This
  function computes exactly that: `kotoba.sm.tensor/mat-inverse` (throws if
  `h` is singular) then `kotoba.sm.tensor/mat-transpose`.

  VERIFICATION (gtg_test.cljc): (a) for a concrete invertible, non-identity
  4x4 matrix, `kotoba.sm.tensor/mat-inverse` genuinely satisfies
  h * h^-1 = h^-1 * h = the 4x4 identity. (b) the biorthogonality relation
  h-bar^mu . h_nu = delta^mu_nu (plain component pairing) holds numerically
  for that same non-identity h. (c) at h = `position-gauge-identity`,
  h-bar = h EXACTLY (the identity matrix is its own inverse and its own
  transpose, so Hbar = (I^-1)^T = I^T = I = H -- bit-for-bit, not merely
  close). (d) a singular (determinant ~0) h throws `ex-info` (propagated
  from `mat-inverse`) rather than silently dividing by ~0 or returning a
  garbage matrix.

  HONESTY NOTE, READ BEFORE USING THIS AS 'the LDG reciprocal frame' for a
  general h: *** this function's biorthogonality pairing is the PLAIN
  component/Kronecker pairing, NOT the Minkowski-metric GA scalar product
  `kotoba.sm.tensor/dot` that a GA-native reciprocal frame would generally be
  defined with. *** Working the Minkowski-pairing version through the same
  derivation: sum_a eta[a][a] Hbar[mu][a] H[nu][a] = delta[mu][nu] is the
  matrix equation Hbar . eta . H^T = I, giving Hbar_eta = (H^-1)^T . eta (eta
  is its own inverse in this basis) -- a DIFFERENT matrix from this
  function's Hbar in general, and in particular at h =
  `position-gauge-identity` it gives Hbar_eta = eta = diag(1,-1,-1,-1), i.e.
  h-bar^0 = h_0 but h-bar^i = -h_i for the 3 spatial rows -- NOT h-bar = h.

  *** FOLLOW-UP LITERATURE CHECK (Phase 0d, stage 2 investigation): this
  function's h-of-h reciprocal frame turns out NOT to be what the actual
  Lasenby-Doran-Gull curvature scalar needs at all. *** Reading the primary
  source directly (Lasenby, Doran & Gull, 'Gravity, Gauge Theories and
  Geometric Algebra', Phil. Trans. R. Soc. Lond. A (1998) 356 487-582,
  updated arXiv version gr-qc/0405033, section 4 'The field equations',
  eqns 4.9, 4.11, 4.12) and a second, independent, mutually-corroborating
  source by two of the same three authors (Lewis, Doran & Lasenby,
  'Quadratic Lagrangians and Topology in Gauge Theory Gravity',
  arXiv gr-qc/9910039, section 2, eqns 5, 6, 12, 14), the Ricci SCALAR is:

    R = sum_{a,b} gamma^a . (gamma^b . R(h_b ^ h_a))

  where h_b = h(gamma_b) (source 1's h_mu, EXACTLY this namespace's Phase 0b
  `h`/`position-gauge-identity` convention -- gr-qc/0405033 eq 3.35's
  g_mu = h^-1(e_mu) is a DIFFERENT, h^-1-built quantity, not this one),
  gamma^a is the reciprocal of the FIXED BACKGROUND orthonormal frame
  gamma_a (TRIVIAL under the Minkowski metric -- gamma^0=gamma_0,
  gamma^i=-gamma_i, i.e. `kotoba.sm.tensor/raise` applied to the ordinary
  basis vectors, NOT a reciprocal frame of h at all), and R(x^y) is the
  BILINEAR extension of Phase 0a's `rotation-field-strength` to an arbitrary
  bivector argument (justified rigorously, not guessed: R_ab is DEFINED as
  D_a Omega_b - D_b Omega_a + Omega_a x Omega_b, manifestly bilinear in a,b
  since Omega(a) is itself linear in a by construction -- so
  R(h_b^h_a) = sum_{rho,sigma} h[b][rho] h[a][sigma] R[rho][sigma] follows
  from `rotation-field-strength`'s own defining bilinearity, needing no new
  literature confirmation).

  Both sources independently agree on this formula's STRUCTURE, and in
  particular BOTH show the outer contraction uses the FIXED frame's OWN
  reciprocal (gamma^a, no inverse of h needed), NOT `reciprocal-frame`
  (h-bar of h) as originally guessed when this function was written -- a
  genuine, literature-grounded correction to that initial assumption.
  HOWEVER, translating 'gamma^b . R(h_b^h_a)' (a vector . bivector -> vector
  Clifford-algebra contraction) into this codebase's R[mu][nu][k]
  six-generator-component array representation (rather than a native GA
  multivector) requires an ADDITIONAL, nontrivial contraction step -- which
  of R's spacetime-index/generator-index slots the contracting vector's
  index pairs with, and with what sign/factor -- that NEITHER source spells
  out in a form directly transcribable to this array representation, and
  that this pass could not independently re-derive with confidence within a
  bounded research effort (a hand re-derivation via the standard GA identity
  v.(a^b)=(v.a)b-(v.b)a produced a candidate formula but with unresolved
  factor-of-2/sign risk once expressed over the SIX generator components
  rather than a full antisymmetric rank-2 slice, and no independent
  numeric-example reference value was found to test a candidate against).
  This is exactly the class of error this codebase's own established
  practice already treats as too risky to guess at without independent
  verification -- see `kotoba.sm.gauge/self-interaction-term`'s
  since-fixed slot-order bug (gtg.cljc's rotation-field-strength docstring)
  and `curvature-quadratic-invariant`'s own honesty note (item 10) for two
  prior instances of the SAME 'which index pairs with which' risk in this
  exact curvature-scalar problem. `curvature-scalar` is therefore
  DELIBERATELY NOT implemented in this namespace -- see the README's Phase
  0d section for the fuller writeup of this investigation, including both
  source URLs. `reciprocal-frame` is offered only as a standard,
  independently-correct dual-basis construction of h's own reciprocal frame
  (useful in its own right, e.g. for a future combined covariant derivative,
  still out of scope -- see the namespace docstring's out-of-scope list),
  not asserted to be an ingredient of the LDG curvature scalar for a general
  h."
  ([h] (reciprocal-frame h 1e-9))
  ([h eps]
   (tensor/mat-transpose (tensor/mat-inverse h eps))))
