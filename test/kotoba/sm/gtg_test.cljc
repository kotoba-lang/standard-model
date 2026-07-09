(ns kotoba.sm.gtg-test
  "Tests for the GTG (Lasenby-Doran-Gull 1998) rotation-gauge sector, Phase
  0a -- see kotoba.sm.gtg's namespace docstring for the exact scope. The
  `trace-normalization-*` tests below are the load-bearing ones: they record,
  numerically, that kotoba.sm.gauge's compact-group trace normalization
  Tr(T^aT^b)=1/2delta^ab does NOT hold for the noncompact so(1,3) bivector
  generators, and exactly how kotoba.sm.gauge/structure-constants's raw
  output diverges from the genuine so(1,3) structure constants as a result.
  These numbers are recorded honestly, not massaged to force agreement."
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.complex :as c]
            [kotoba.sm.spinor :as spinor]
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

(deftest raw-structure-constants-diverge-from-genuine-ones
  (testing "kotoba.sm.gauge/structure-constants applied as-is to these generators gives
            f-raw[A][B][C] = 2*K[C][C]*f-true[A][B][C] -- i.e. a uniform +2 factor when C is a
            rotation-type generator, and a WRONG-SIGNED -2 factor when C is boost-type"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)
          K (gtg/generator-trace-gram)]
      (doseq [A (range 6) B (range 6) C (range 6)]
        (is (close? (get-in f-raw [A B C])
                     (* 2 (get-in K [C C]) (get-in f-true [A B C])))
            (str "f-raw[" A "][" B "][" C "] vs 2*K[C][C]*f-true")))))
  (testing "concretely: [T^{01},T^{12}] has genuine coefficient -1 on T^{02} (basis index 1), but
            kotoba.sm.gauge/structure-constants's raw output for that same triple is +2 -- both
            the magnitude (2x) and the SIGN are wrong, because T^{02} is boost-type (K=-1)"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)]
      (is (close? (get-in f-true [0 3 1]) -1.0))
      (is (close? (get-in f-raw [0 3 1]) 2.0))
      (is (not (close? (get-in f-raw [0 3 1]) (get-in f-true [0 3 1]))))))
  (testing "by contrast, for a rotation-type third index the raw output has the right SIGN, just
            2x the genuine magnitude: [T^{12},T^{13}] on T^{23} (basis index 5)"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)]
      (is (close? (get-in f-true [3 4 5]) 1.0))
      (is (close? (get-in f-raw [3 4 5]) 2.0))))
  (testing "raw structure constants never literally equal the genuine ones for any nonzero triple
            (the naive Tr=1/2delta reuse is not merely off by a global constant that could be
            divided out -- the boost/rotation sign split makes it structurally wrong)"
    (let [f-raw (gtg/rotation-raw-structure-constants)
          f-true (gtg/true-structure-constants)
          nonzero-triples (for [A (range 6) B (range 6) C (range 6)
                                 :when (> (Math/abs (double (get-in f-true [A B C]))) 1e-9)]
                             [A B C])]
      (is (= 24 (count nonzero-triples)))
      (doseq [[A B C] nonzero-triples]
        (is (not (close? (get-in f-raw [A B C]) (get-in f-true [A B C])))
            (str "[" A " " B " " C "] raw should NOT equal true"))))))

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
