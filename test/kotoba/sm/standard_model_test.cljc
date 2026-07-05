(ns kotoba.sm.standard-model-test
  (:require [clojure.test :refer [deftest is testing]]
            [kotoba.sm.complex :as c]
            [kotoba.sm.tensor :as tensor]
            [kotoba.sm.gauge :as gauge]
            [kotoba.sm.standard-model :as sm]))

(defn- close? [a b] (< (Math/abs (double (- a b))) 1e-6))

(deftest fermion-table
  (testing "the fermion table has 3 generations x 7 entries (2 doublet + 2 up/down-singlet quark, 2 doublet + 1 singlet lepton)"
    (is (= 21 (count (sm/fermion-content)))))
  (testing "every entry satisfies Gell-Mann-Nishijima Q = T3 + Y/2"
    (is (sm/gell-mann-nishijima-ok?)))
  (testing "the up quark has charge +2/3"
    (is (close? (:Q (first (filter #(and (= (:name %) "u") (= (:chirality %) :L)) (sm/fermion-content)))) (/ 2.0 3.0))))
  (testing "the electron has charge -1"
    (is (close? (:Q (first (filter #(and (= (:name %) "e") (= (:chirality %) :L)) (sm/fermion-content)))) -1.0)))
  (testing "the left-handed neutrino has charge 0"
    (is (close? (:Q (first (filter #(and (= (:name %) "nu_e") (= (:chirality %) :L)) (sm/fermion-content)))) 0.0))))

(deftest electroweak-mixing
  (testing "the mixing matrix is orthogonal (rotation): R^T R = I"
    (let [theta-w (sm/weinberg-angle 0.65 0.36)
          R (sm/mixing-matrix theta-w)
          RT-R (tensor/mat-mat (tensor/mat-transpose R) R)]
      (is (close? (get-in RT-R [0 0]) 1.0))
      (is (close? (get-in RT-R [1 1]) 1.0))
      (is (close? (get-in RT-R [0 1]) 0.0))))
  (testing "M_W = M_Z cos(theta_W) (tree-level relation, derived from couplings + vev, not assumed)"
    (let [g 0.65 g-prime 0.36 v 246.0
          theta-w (sm/weinberg-angle g g-prime)
          mw (sm/w-mass g v)
          mz (sm/z-mass g g-prime v)]
      (is (close? mw (* mz (Math/cos theta-w))))))
  (testing "photon and Z are an orthonormal recombination of B and W3 (energy/norm preserved)"
    (let [theta-w (sm/weinberg-angle 0.65 0.36)
          [A Z] (sm/photon-Z 1.0 2.0 theta-w)]
      (is (close? (+ (* A A) (* Z Z)) (+ (* 1.0 1.0) (* 2.0 2.0)))))))

(deftest higgs-mechanism
  (testing "higgs-vev and higgs-mass are consistent: sqrt(2 lambda) v = sqrt(-2 mu2)"
    (let [mu2 -100.0 lambda 0.13
          v (sm/higgs-vev mu2 lambda)
          mh (sm/higgs-mass mu2 lambda)]
      (is (close? mh (* (Math/sqrt (* 2 lambda)) v)))))
  (testing "the potential V(x)=mu2 x + lambda x^2 is minimized at x_min = -mu2/(2 lambda)
            (derivative-free check: V(x_min) < V(0.9 x_min) and < V(1.1 x_min))"
    (let [mu2 -100.0 lambda 0.13
          x-min (/ (- mu2) (* 2 lambda))]
      (is (< (sm/higgs-potential mu2 lambda x-min) (sm/higgs-potential mu2 lambda (* 0.9 x-min))))
      (is (< (sm/higgs-potential mu2 lambda x-min) (sm/higgs-potential mu2 lambda (* 1.1 x-min))))))
  (testing "Yukawa mass generation m_f = y_f v / sqrt(2)"
    (is (close? (sm/yukawa-mass 1.0 246.0) (/ 246.0 (Math/sqrt 2.0))))))

(deftest ckm-unitarity
  (testing "the CKM matrix (PDG parameterization) is exactly unitary: V-dagger V = I"
    (let [V (sm/ckm-matrix 0.227 0.0035 0.0413 1.2)
          Vd (c/m-dagger V)
          prod (c/m-mul Vd V)]
      (is (c/m-approx= prod (c/m-identity 3) 1e-9)))))

(deftest covariant-derivative-sm-reduces-to-partial-when-fields-off
  (testing "for a color-singlet, weak-doublet lepton with all gauge fields off, D_mu = d_mu"
    (let [fermion {:color :singlet :weak :doublet :Y -1.0}
          psi [(c/c 1 0) (c/c 0 1)]
          d-psi [(c/c 0.1 0) (c/c 0 0.2)]
          D (sm/covariant-derivative-sm d-psi fermion [] [0.0 0.0 0.0] 0.0 1.0 1.0 0.6 psi)]
      (is (c/v-approx= D d-psi))))
  (testing "for a color-triplet, weak-singlet quark with all gauge fields off, D_mu = d_mu"
    (let [fermion {:color :triplet :weak :singlet :Y (/ 4.0 3.0)}
          psi [(c/c 1 0) (c/c 0 1) (c/c 1 1)]
          d-psi [(c/c 0.1 0) (c/c 0 0.2) (c/c 0 0)]
          D (sm/covariant-derivative-sm d-psi fermion (vec (repeat 8 0.0)) [] 0.0 1.2 0.65 0.36 psi)]
      (is (c/v-approx= D d-psi)))))

(deftest lagrangian-density-pieces
  (testing "Yang-Mills density of a zero field strength is zero"
    (let [F (vec (repeat 4 (vec (repeat 4 (vec (repeat 8 0.0))))))]
      (is (close? (sm/yang-mills-density F) 0.0))))
  (testing "Yang-Mills density for U(1) (no generators) is zero by convention"
    (let [F (vec (repeat 4 (vec (repeat 4 []))))]
      (is (close? (sm/yang-mills-density F) 0.0))))
  (testing "total-lagrangian-density sums all pieces"
    (is (c/approx= (sm/total-lagrangian-density [1.0 2.0] [(c/c 3 0)] (c/c 4 0) [(c/c 5 0)])
                    (c/c 15 0)))))
