(ns kotoba.sm.gtg-test
  "Tests for the GTG (Lasenby-Doran-Gull 1998) rotation-gauge sector (Phase
  0a) and position-gauge sector (Phase 0b) -- see kotoba.sm.gtg's namespace
  docstring for the exact scope. The `trace-normalization-*` tests below are
  the load-bearing ones: they record, numerically, that kotoba.sm.gauge's
  (Phase-0a-era, compact-group-only) trace normalization
  Tr(T^aT^b)=1/2delta^ab does NOT hold for the noncompact so(1,3) bivector
  generators.

  `structure-constants-now-match-genuine-values-after-gauge-fix` (formerly
  `raw-structure-constants-diverge-from-genuine-ones`, see that deftest's own
  docstring for the before/after) records the CONSEQUENCE of that fact as it
  stood at Phase 0a time: `kotoba.sm.gauge/structure-constants`'s raw output
  used to diverge from the genuine so(1,3) structure constants. That has
  SINCE BEEN FIXED by generalizing `kotoba.sm.gauge/structure-constants` to
  compute its own generator set's trace-Gram matrix instead of assuming it is
  1/2delta^AB -- so raw and genuine now MATCH. These numbers are recorded
  honestly either way, not massaged to force agreement.

  `flat-limit-h-reproduces-minkowski-metric-exactly` and
  `global-lorentz-transformations-are-the-constant-flat-special-case`
  (Phase 0b, position-gauge sector) are the load-bearing tests for the
  position gauge field h_mu and its derived metric g_mu-nu -- see
  kotoba.sm.gtg's Phase-0b SCOPE note (items 7-9) for the exact claims being
  checked."
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.complex :as c]
            [kotoba.sm.spinor :as spinor]
            [kotoba.sm.tensor :as tensor]
            [kotoba.sm.vector-field :as vf]
            [kotoba.sm.gtg :as gtg]))

(defn- close? [a b] (< (Math/abs (double (- a b))) 1e-6))

;; ---------------------------------------------------------------------------
;; 1+2. generators are built from the existing gamma matrices, and close the
;;      so(1,3) Lie algebra: [T^{ab},T^{cd}] = i(eta^{bc}T^{ad} - eta^{ac}T^{bd}
;;      - eta^{bd}T^{ac} + eta^{ad}T^{bc})
;; ---------------------------------------------------------------------------

(deftest generators-are-4x4-complex-matrices
  (testing "there are 6 generators, each a 4x4 complex matrix"
    (is (= 6 (count gtg/generators)))
    (is (every? #(= [4 4] (c/m-dims %)) gtg/generators)))
  (testing "boost-type generators T^{0i} (index 0,1,2) are ANTI-Hermitian, not Hermitian -- another
            noncompactness symptom, and the textbook reason the Dirac spinor rep of the Lorentz
            group is not unitary (Peskin & Schroeder section 3.2's S^{mu-nu}=(i/4)[gamma^mu,gamma^nu]
            is literally this namespace's T^{ab}, and they note the same anti-Hermitian boost /
            Hermitian rotation split)"
    (doseq [A [0 1 2]]
      (let [T (nth gtg/generators A)]
        (is (c/m-approx= T (c/m-rscale -1 (c/m-dagger T)) 1e-9) (str "T index " A " anti-Hermitian")))))
  (testing "rotation-type generators T^{ij} (index 3,4,5, purely spatial pairs) ARE Hermitian,
            same as the compact-group generators kotoba.sm.gauge was originally tested against"
    (doseq [A [3 4 5]]
      (is (c/m-hermitian? (nth gtg/generators A) 1e-9) (str "T index " A " Hermitian")))))

(deftest lorentz-algebra-closure
  (testing "a single concrete pair: [T^{01},T^{02}] against the closed-form RHS"
    (is (c/m-approx= (c/m-commutator (gtg/generator 0 1) (gtg/generator 0 2))
                      (gtg/lorentz-algebra-rhs 0 1 0 2)
                      1e-9)))
  (testing "all 36 (a<b, c<d) generator pairs satisfy the so(1,3) commutation relation"
    (doseq [[a b] gtg/generator-index-pairs
            [cc d] gtg/generator-index-pairs]
      (is (c/m-approx= (c/m-commutator (gtg/generator a b) (gtg/generator cc d))
                        (gtg/lorentz-algebra-rhs a b cc d)
                        1e-9)
          (str "[T^" a b ",T^" cc d "]")))))

;; ---------------------------------------------------------------------------
;; 3. THE load-bearing honesty check: Tr(T^aT^b)=1/2delta^ab does NOT hold for
;;    these noncompact so(1,3) generators.
;; ---------------------------------------------------------------------------

(deftest trace-normalization-does-not-match-compact-group-convention
  (testing "kotoba.sm.gauge's docstring convention Tr(T^aT^b)=1/2delta^ab does NOT hold here"
    (is (false? (gtg/compact-group-trace-normalization-holds?))))
  (testing "the actual Gram matrix K[A][B]=Tr(T_A T_B) is diagonal (off-diagonal ~0, same
            structure as the compact-group case) but +/-1 on the diagonal, not a uniform +1/2"
    (let [K (gtg/generator-trace-gram)]
      (doseq [A (range 6) B (range 6)]
        (when (not= A B)
          (is (close? (get-in K [A B]) 0.0) (str "K[" A "][" B "] off-diagonal"))))
      (testing "boost-type generators (index 0,1,2 = pairs containing timelike index 0): Tr(T T) = -1"
        (doseq [A [0 1 2]]
          (is (close? (get-in K [A A]) -1.0) (str "K[" A "][" A "]"))))
      (testing "rotation-type generators (index 3,4,5 = purely spatial pairs): Tr(T T) = +1"
        (doseq [A [3 4 5]]
          (is (close? (get-in K [A A]) 1.0) (str "K[" A "][" A "]"))))))
  (testing "in particular the magnitude is 1, not 1/2, for every generator (a further, independent
            way the compact-group assumption fails beyond the boost-generator sign flip)"
    (let [K (gtg/generator-trace-gram)]
      (doseq [A (range 6)]
        (is (close? (Math/abs (double (get-in K [A A]))) 1.0))))))

(deftest structure-constants-now-match-genuine-values-after-gauge-fix
  "BEFORE/AFTER: this deftest used to be named
  `raw-structure-constants-diverge-from-genuine-ones` and asserted the
  OPPOSITE of what it asserts now. At Phase 0a time,
  `kotoba.sm.gauge/structure-constants` hard-coded the compact-group
  assumption Tr(T^aT^b)=1/2delta^ab, so its raw output for this so(1,3)
  generator set (`rotation-raw-structure-constants`) was f-raw[A][B][C] =
  2*K[C][C]*f-true[A][B][C] -- a uniform +2 factor when C is a rotation-type
  generator, and a WRONG-SIGNED -2 factor when C is boost-type -- and never
  literally equaled the genuine values `true-structure-constants` for any
  nonzero triple. `kotoba.sm.gauge/structure-constants` has SINCE BEEN FIXED
  to compute its own generator set's trace-Gram matrix K[A][B]=Tr(T_AT_B)
  instead of assuming it is 1/2delta^AB (see the `kotoba.sm.gauge` namespace
  docstring), so raw and true now MATCH -- the assertions below are the same
  shape as the old ones, with `not (close? ...)` flipped to `close? ...`, so
  the diff against the old test makes the before/after explicit."
  (testing "kotoba.sm.gauge/structure-constants applied to these generators now gives
            f-raw[A][B][C] = f-true[A][B][C] exactly (within floating tolerance) for
            every triple -- the fix eliminated the old 2*K[C][C] divergence entirely"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)]
      (doseq [A (range 6) B (range 6) C (range 6)]
        (is (close? (get-in f-raw [A B C]) (get-in f-true [A B C]))
            (str "f-raw[" A "][" B "][" C "] vs f-true")))))
  (testing "and, as a check that the fix isn't accidentally satisfying the OLD (broken)
            2*K[C][C] relationship too (which would only be possible if 2*K[C][C]=1, false for
            every K[C][C]=+/-1 in this generator set), raw no longer equals the historically
            broken value the old hard-coded -2i formula used to produce, for any nonzero triple"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)
          K (gtg/generator-trace-gram)
          nonzero-triples (for [A (range 6) B (range 6) C (range 6)
                                 :when (> (Math/abs (double (get-in f-true [A B C]))) 1e-9)]
                             [A B C])]
      (is (= 24 (count nonzero-triples)))
      (doseq [[A B C] nonzero-triples]
        (let [old-broken-value (* 2 (get-in K [C C]) (get-in f-true [A B C]))]
          (is (not (close? (get-in f-raw [A B C]) old-broken-value))
              (str "[" A " " B " " C "] raw should NOT equal the pre-fix 2*K[C][C]*f-true value"))))))
  (testing "concretely: [T^{01},T^{12}] has genuine coefficient -1 on T^{02} (basis index 1); the
            raw output for that same triple is now -1 too (Phase 0a: it used to be the
            wrong-signed +2, since T^{02} is boost-type, K=-1)"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)]
      (is (close? (get-in f-true [0 3 1]) -1.0))
      (is (close? (get-in f-raw [0 3 1]) -1.0))
      (is (close? (get-in f-raw [0 3 1]) (get-in f-true [0 3 1])))))
  (testing "and for a rotation-type third index: [T^{12},T^{13}] on T^{23} (basis index 5) -- raw
            is now +1, matching true (Phase 0a: it used to be +2, right sign but 2x magnitude)"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)]
      (is (close? (get-in f-true [3 4 5]) 1.0))
      (is (close? (get-in f-raw [3 4 5]) 1.0)))))

;; ---------------------------------------------------------------------------
;; 4. rotation gauge field Omega_mu and its field strength R_mu-nu, via
;;    kotoba.sm.gauge/field-strength applied to `generators`.
;; ---------------------------------------------------------------------------

(deftest rotation-field-strength-curl-term-matches-yang-mills-form
  (testing "with the gauge-field VALUE zero at the point (so the self-interaction term
            g f^abc Omega_mu^b Omega_nu^c vanishes regardless of the structure-constant
            normalization caveat), R_mu-nu = d_mu Omega_nu - d_nu Omega_mu exactly -- the same
            functional form as kotoba.sm.gauge's abelian-U(1)/non-abelian curl term
            (mirrors gauge_test.cljc's u1-field-strength-is-ordinary-em-tensor)"
    (let [zero6 (vec (repeat 6 0.0))
          zero-nu (vec (repeat 4 zero6))
          d-Omega (-> (vec (repeat 4 zero-nu))
                      (assoc-in [0 0] (assoc zero6 3 1.0))   ;; d/dt Omega_t^{12} = 1
                      (assoc-in [1 0] (assoc zero6 3 2.0))   ;; d/dx Omega_t^{12} = 2
                      (assoc-in [2 1] (assoc zero6 3 3.0)))  ;; d/dy Omega_x^{12} = 3
          Omega (vec (repeat 4 zero6))
          R (gtg/rotation-field-strength d-Omega Omega 1.0)]
      (is (close? (get-in R [1 0 3]) (- (get-in d-Omega [1 0 3]) (get-in d-Omega [0 1 3]))))
      (is (close? (get-in R [1 0 3]) 2.0))
      (is (close? (get-in R [0 1 3]) -2.0))
      (is (close? (get-in R [2 1 3]) 3.0)))))

(deftest rotation-field-strength-vanishes-for-constant-single-component-field
  (testing "a CONSTANT field with only one bivector component turned on has zero curl (constant)
            AND zero self-interaction (f^{akk}=0 by total antisymmetry, true regardless of the
            raw-vs-true normalization question) -- a normalization-agnostic sanity check that
            rotation-gauge-field-gradient + rotation-field-strength compose correctly end to end"
    (let [k 2
          val 0.7
          Omega-field (fn [_x] (vec (repeat 4 (assoc (vec (repeat 6 0.0)) k val))))
          x [0.1 0.2 0.3 0.4]
          d-Omega (gtg/rotation-gauge-field-gradient Omega-field x)
          Omega (Omega-field x)
          R (gtg/rotation-field-strength d-Omega Omega 1.0)]
      (is (every? #(close? % 0.0) (flatten R))))))

;; ---------------------------------------------------------------------------
;; 5+6. spin-connection covariant derivative on a Dirac spinor, and the flat
;;      limit against kotoba.sm.spinor's free Dirac equation
;; ---------------------------------------------------------------------------

(deftest spinor-covariant-derivative-reduces-to-partial-when-omega-off
  (testing "D_mu psi = d_mu psi when Omega_mu^k = 0 for all 6 k"
    (let [psi [(c/c 1 0) (c/c 0 1) (c/c 0.5 0.5) (c/c -1 0)]
          d-psi [(c/c 0.1 0) (c/c 0 0.2) (c/c 0.3 0) (c/c 0 -0.1)]
          D (gtg/spinor-covariant-derivative d-psi (vec (repeat 6 0.0)) 1.0 psi)]
      (is (c/v-approx= D d-psi 1e-12)))))

(deftest flat-limit-matches-free-dirac-equation-exactly
  (testing "with Omega_mu = 0 everywhere, the covariant Dirac residual built on
            spinor-covariant-derivative is EXACTLY kotoba.sm.spinor/dirac-residual on the same
            plane-wave field, mass, and point -- not merely close within finite-difference error,
            but numerically identical (gauge-correction is the literal zero vector when every
            gauge-field component is 0, so the two finite-difference computations coincide term
            for term)"
    (binding [spinor/*h* 1e-5]
      (let [m 1.0 p [0.3 0.1 0.2]
            psi-field (spinor/plane-wave :up p m)
            Omega-mu-field (constantly (vec (repeat 4 (vec (repeat 6 0.0)))))
            x [0.4 0.1 -0.2 0.3]
            r-plain (spinor/dirac-residual psi-field m x)
            r-gauge (gtg/dirac-residual-with-rotation-gauge psi-field Omega-mu-field m 1.0 x)]
        (is (c/v-approx= r-plain r-gauge 1e-15))
        (testing "and (as an independent existing-namespace cross-check) both residuals are
                  themselves ~0, since plane-wave is an exact free-Dirac-equation solution"
          (is (every? (fn [z] (< (c/abs2 z) 1e-6)) r-plain))
          (is (every? (fn [z] (< (c/abs2 z) 1e-6)) r-gauge)))))))

;; ---------------------------------------------------------------------------
;; 7+8+9. position gauge field h_mu(x) and the derived metric g_mu-nu(x),
;;    Phase 0b
;; ---------------------------------------------------------------------------

(deftest flat-limit-h-reproduces-minkowski-metric-exactly
  (testing "h = position-gauge-identity (h_mu^nu = delta_mu^nu) gives g_mu-nu
            EXACTLY kotoba.sm.tensor/metric -- bit-for-bit (integer arithmetic
            throughout), not merely close within a tolerance"
    (let [g (gtg/derived-metric gtg/position-gauge-identity)]
      (is (= g tensor/metric))
      (doseq [mu (range 4) nu (range 4)]
        (is (= (get-in g [mu nu]) (get-in tensor/metric [mu nu]))
            (str "g[" mu "][" nu "] exact match")))))
  (testing "derived-metric-matches-flat? agrees, at the default epsilon"
    (is (true? (gtg/derived-metric-matches-flat? gtg/position-gauge-identity))))
  (testing "derived-metric-field composes h-field + point evaluation correctly for a
            constant identity field, at several distinct spacetime points"
    (let [h-field (constantly gtg/position-gauge-identity)]
      (doseq [x [[0 0 0 0] [1.0 2.0 -3.0 0.5] [-4.0 0.0 1.0 7.0]]]
        (is (= tensor/metric (gtg/derived-metric-field h-field x)))))))

(deftest global-lorentz-transformations-are-the-constant-flat-special-case
  (testing "a constant position-gauge field built from ANY existing kotoba.sm.vector-field
            Lorentz boost/rotation matrix reproduces the flat metric exactly (within
            floating-point roundoff from sqrt/cos/sin, not merely a loose 'close') -- the
            mathematical content is Lambda eta Lambda^T = eta, which follows algebraically
            from Lambda^T eta Lambda = eta (vf/lorentz?'s defining property), see this
            namespace's Phase-0b section-7 header comment for the derivation"
    (doseq [lambda [(vf/boost-x 0.6)
                     (vf/boost-y 0.3)
                     (vf/boost-z -0.45)
                     (vf/boost [0.3 0.2 0.1])
                     (vf/rotation-x 0.9)
                     (vf/rotation-y 1.234)
                     (vf/rotation-z -2.1)]]
      (is (vf/lorentz? lambda) "sanity: lambda actually IS a Lorentz transformation")
      (is (gtg/derived-metric-matches-flat? lambda 1e-9)
          (str "derived metric from constant h=" lambda " should equal the flat metric"))))
  (testing "concretely, for a single boost, every derived-metric component matches
            kotoba.sm.tensor/metric numerically"
    (let [lambda (vf/boost-x 0.6)
          g (gtg/derived-metric lambda)]
      (doseq [mu (range 4) nu (range 4)]
        (is (close? (get-in g [mu nu]) (get-in tensor/metric [mu nu]))
            (str "g[" mu "][" nu "]")))))
  (testing "the identity Lorentz transformation (zero boost, vf/boost's own literal-integer
            identity-matrix branch) reproduces g_mu-nu EXACTLY, the same bit-for-bit
            equality as position-gauge-identity itself -- not a coincidence, since
            vf/boost's zero-velocity branch literally returns position-gauge-identity's
            same 4x4 integer identity matrix"
    (let [g (gtg/derived-metric (vf/boost [0 0 0]))]
      (is (= g tensor/metric)))))
