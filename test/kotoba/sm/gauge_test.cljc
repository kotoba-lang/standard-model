(ns kotoba.sm.gauge-test
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.complex :as c]
            [kotoba.sm.gauge :as g]))

(defn- close? [a b] (< (Math/abs (double (- a b))) 1e-6))

;; ---------------------------------------------------------------------------
;; backward-compatibility regression: `structure-constants` was generalized
;; (kotoba.sm.gauge namespace docstring) from a hard-coded
;; f^abc = -2i Tr([T^a,T^b]T^c) formula -- which silently assumed
;; Tr(T^aT^b)=1/2delta^ab -- to a generic f^{ABD} = -i Tr([T_A,T_B]T_D) / K[D][D]
;; formula that computes the generators' own trace-Gram matrix K instead of
;; assuming it. `su2-structure-constants`/`su3-structure-constants` below
;; (pre-existing) already regression-test against textbook values; the tests
;; in this section additionally regression-test the new generic formula
;; against an independent re-implementation of the OLD hard-coded formula
;; itself, proving the generalization changed nothing for the
;; Tr(T^aT^b)=1/2delta^ab generator sets it always used to run against.
;; ---------------------------------------------------------------------------

(defn- legacy-half-normalized-structure-constants
  "Re-implementation, kept ONLY in this test file, of the formula
  `kotoba.sm.gauge/structure-constants` used before it was generalized to
  compute the generators' own trace-Gram matrix: f^abc = -2i Tr([T^a,T^b]T^c),
  valid only under the Tr(T^aT^b)=1/2delta^ab convention (Pauli/2,
  Gell-Mann/2). An independent oracle for the regression tests below -- NOT
  itself exercised by any production code path."
  [generators]
  (let [n (count generators)]
    (vec (for [a (range n)]
           (vec (for [b (range n)]
                  (vec (for [cc (range n)]
                         (let [comm (c/m-commutator (nth generators a) (nth generators b))
                               prod (c/m-mul comm (nth generators cc))
                               tr (c/m-trace prod)]
                           (c/re (c/c* (c/c 0 -2) tr)))))))))))

(deftest structure-constants-matches-legacy-half-normalized-formula
  (testing "SU(2): the new generic Gram-matrix formula reproduces the OLD hard-coded
            f^{abc}=-2i Tr([T^a,T^b]T^c) formula for every triple (within floating
            tolerance), since Tr(T^aT^b)=1/2delta^ab holds for the Pauli/2 generators"
    (let [new-f (g/structure-constants g/su2-generators)
          old-f (legacy-half-normalized-structure-constants g/su2-generators)]
      (doseq [a (range 3) b (range 3) cc (range 3)]
        (is (close? (get-in new-f [a b cc]) (get-in old-f [a b cc]))
            (str "f[" a "][" b "][" cc "]")))))
  (testing "SU(3): same equivalence for the Gell-Mann/2 generators"
    (let [new-f (g/structure-constants g/su3-generators)
          old-f (legacy-half-normalized-structure-constants g/su3-generators)]
      (doseq [a (range 8) b (range 8) cc (range 8)]
        (is (close? (get-in new-f [a b cc]) (get-in old-f [a b cc]))
            (str "f[" a "][" b "][" cc "]"))))))

(deftest gram-matrix-is-uniform-half-for-su2-and-su3
  (testing "the generalization's premise, checked directly: K[A][B]=Tr(T_AT_B) is
            1/2delta^AB for both SU(2) Pauli/2 and SU(3) Gell-Mann/2 generators
            (this is WHY the new generic formula collapses to the old hard-coded one)"
    (doseq [gens [g/su2-generators g/su3-generators]]
      (let [K (g/generator-gram gens)
            n (count gens)]
        (doseq [A (range n) B (range n)]
          (is (close? (c/re (get-in K [A B])) (if (= A B) 0.5 0.0))
              (str "K[" A "][" B "] re"))
          (testing "imaginary part of the (complex-valued, in general) trace is ~0"
            (is (close? (c/im (get-in K [A B])) 0.0)
                (str "K[" A "][" B "] im"))))))))

(deftest non-diagonal-gram-throws
  (testing "structure-constants explicitly scopes out non-orthogonal generator bases
            (a non-diagonal trace-Gram matrix K[A][B]=Tr(T_AT_B)) by throwing rather
            than silently returning wrong values -- exercised here with two identical
            (hence non-orthogonal, K[0][1]=Tr(I*I)=2≠0) 2x2 'generators'"
    (let [I (c/m-identity 2)]
      (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                    (g/structure-constants [I I]))))))

(deftest su2-structure-constants
  (testing "f^123 = epsilon^123 = 1 for SU(2) (0-indexed: f[0][1][2])"
    (let [f (g/structure-constants g/su2-generators)]
      (is (close? (get-in f [0 1 2]) 1.0))
      (is (close? (get-in f [1 0 2]) -1.0))
      (is (close? (get-in f [0 0 2]) 0.0)))))

(deftest su3-structure-constants
  (testing "f^123 = 1 (0-indexed [0 1 2])"
    (let [f (g/structure-constants g/su3-generators)]
      (is (close? (get-in f [0 1 2]) 1.0))))
  (testing "f^458 = sqrt(3)/2 (0-indexed [3 4 7])"
    (let [f (g/structure-constants g/su3-generators)]
      (is (close? (get-in f [3 4 7]) (/ (Math/sqrt 3.0) 2.0)))))
  (testing "f^147 = 1/2 (0-indexed [0 3 6])"
    (let [f (g/structure-constants g/su3-generators)]
      (is (close? (get-in f [0 3 6]) 0.5)))))

(deftest u1-abelian
  (testing "U(1) is abelian: the single generator commutes with itself, so f=0"
    (let [gens (g/u1-generators 1.0 1)
          f (g/structure-constants gens)]
      (is (close? (get-in f [0 0 0]) 0.0)))))

(deftest covariant-derivative-reduces-to-partial-when-field-off
  (testing "D_mu psi = d_mu psi when A_mu^a = 0"
    (let [psi [(c/c 1 0) (c/c 0 1)]
          d-psi [(c/c 0.1 0) (c/c 0 0.2)]
          D (g/covariant-derivative d-psi g/su2-generators [0.0 0.0 0.0] 1.0 psi)]
      (is (c/v-approx= D d-psi)))))

(deftest u1-field-strength-is-ordinary-em-tensor
  (testing "for abelian U(1), F_mu-nu = d_mu A_nu - d_nu A_mu (no self-interaction term)"
    ;; A_nu(x) = [t+2x, 3y, 0, 0] (single U(1) component); pick d-A directly
    ;; (bypassing finite differences) to check the algebra in field-strength itself.
    (let [f-abc (g/structure-constants (g/u1-generators 1.0 1))
          ;; d-A[mu][nu][a]: only nonzero entries are d/dt(A_t)=1 (mu=0,nu=0) and
          ;; d/dx(A_t)=2 (mu=1,nu=0), d/dy(A_x)=3 (mu=2,nu=1)
          zero3 (vec (repeat 4 (vec (repeat 4 [0.0]))))
          d-A (-> zero3
                  (assoc-in [0 0] [1.0])
                  (assoc-in [1 0] [2.0])
                  (assoc-in [2 1] [3.0]))
          A [[0.0] [0.0] [0.0] [0.0]] ;; field value itself irrelevant for abelian F (g f^abc term vanishes)
          F (g/field-strength f-abc 1.0 d-A A)]
      (is (close? (get-in F [1 0 0]) (- (get-in d-A [1 0 0]) (get-in d-A [0 1 0]))))
      (is (close? (get-in F [1 0 0]) 2.0))
      (is (close? (get-in F [0 1 0]) -2.0))
      (is (close? (get-in F [2 1 0]) 3.0)))))
