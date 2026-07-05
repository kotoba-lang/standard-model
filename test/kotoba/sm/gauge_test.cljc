(ns kotoba.sm.gauge-test
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.complex :as c]
            [kotoba.sm.gauge :as g]))

(defn- close? [a b] (< (Math/abs (double (- a b))) 1e-6))

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
