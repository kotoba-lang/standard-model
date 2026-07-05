(ns kotoba.sm.tensor-test
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.tensor :as t]))

(defn- close? [a b] (< (Math/abs (double (- a b))) 1e-9))

(deftest metric-and-invariants
  (testing "a timelike unit vector has norm2 = +1"
    (is (close? (t/norm2 [1 0 0 0]) 1.0)))
  (testing "a spacelike unit vector has norm2 = -1"
    (is (close? (t/norm2 [0 1 0 0]) -1.0)))
  (testing "lowering then raising a vector is the identity (eta is its own inverse)"
    (is (= (t/raise (t/lower [3 1 4 1])) [3 1 4 1])))
  (testing "dot product of a four-momentum at rest gives m^2"
    (is (close? (t/dot [5 0 0 0] [5 0 0 0]) 25.0))))

(deftest rank2-tensors
  (testing "outer product then lowering both indices matches direct component computation"
    (let [v [1 2 3 4] w [5 6 7 8]
          T (t/outer v w)
          T-lower (t/lower2 T)]
      ;; T_mu-nu = eta_mu-mu eta_nu-nu T^mu-nu (diagonal metric, no sum)
      (is (close? (get-in T-lower [0 0]) (* 1 1 (get-in T [0 0]))))
      (is (close? (get-in T-lower [1 1]) (* -1 -1 (get-in T [1 1]))))
      (is (close? (get-in T-lower [0 1]) (* 1 -1 (get-in T [0 1]))))))
  (testing "antisymmetrize kills a symmetric tensor"
    (let [v [1 2 3 4]
          sym (t/outer v v)]
      (is (every? #(close? % 0.0) (flatten (t/antisymmetrize sym))))))
  (testing "full-contract of eta^mu-nu with eta_mu-nu is the spacetime dimension (4)"
    (is (close? (t/full-contract t/metric t/metric) 4.0))))

(deftest levi-civita
  (testing "epsilon^0123 = +1 (identity permutation)"
    (is (= (t/levi-civita4 0 1 2 3) 1)))
  (testing "a single transposition flips the sign"
    (is (= (t/levi-civita4 1 0 2 3) -1)))
  (testing "a repeated index gives 0"
    (is (= (t/levi-civita4 0 0 1 2) 0)))
  (testing "a 4-cycle (0123)->(1230) is an odd permutation, gives -1"
    (is (= (t/levi-civita4 1 2 3 0) -1))))
