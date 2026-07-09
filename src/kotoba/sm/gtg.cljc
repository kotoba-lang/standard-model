(ns kotoba.sm.gtg
  "Gauge Theory Gravity (Lasenby-Doran-Gull 1998, 'Gravity, gauge theories and
  geometric algebra') -- ROTATION-GAUGE SECTOR ONLY, Phase 0a. See the
  'Scope' section of the repo README for the ADR reference.

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
       unmodified (`rotation-field-strength`).
    5. the spin-connection covariant derivative on a Dirac spinor,
       D_mu psi = d_mu psi - i g Omega_mu^{ab} T_{ab} psi, obtained by
       handing this generator set to `kotoba.sm.gauge/covariant-derivative`
       unmodified (`spinor-covariant-derivative`).
    6. a flat-limit check (in gtg_test.cljc) that with Omega_mu = 0
       everywhere, the covariant Dirac residual built on top of (5) is
       EXACTLY (not just approximately) `kotoba.sm.spinor`'s existing free
       Dirac-equation residual (`dirac-residual-with-rotation-gauge`).

  EXPLICITLY OUT OF SCOPE, NOT IMPLEMENTED HERE: the position gauge field
  h_mu, the derived spacetime metric g_mu-nu = h_mu . h_nu, the Riemann/Ricci
  curvature scalars built from R_mu-nu, the Einstein multivector, the GTG
  action principle and field equations, any proof of equivalence to General
  Relativity, and any dark-matter/dark-energy/de-Sitter extension. Those are
  legitimate, much larger follow-up efforts, deliberately deferred to a later
  phase -- this namespace is a limited, literature-faithful port of the
  rotation-gauge SECTOR ONLY, and must not be read as 'GTG is implemented
  here' or as any kind of gravity engine."
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
  returns the genuine so(1,3) structure constants."
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
