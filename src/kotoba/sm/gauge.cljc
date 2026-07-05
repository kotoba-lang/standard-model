(ns kotoba.sm.gauge
  "Generic compact-Lie-algebra layer: generator sets for U(1)/SU(2)/SU(3),
  structure constants derived generically (not hard-coded) from any Hermitian
  generator set normalized Tr(T^a T^b) = 1/2 delta^ab (the standard convention
  for the Pauli/2 and Gell-Mann/2 generators used here), the gauge covariant
  derivative, and the non-abelian Yang-Mills field-strength tensor."
  (:require [kotoba.sm.complex :as c]
            [kotoba.sm.spinor :as spinor]))

;; ---------------------------------------------------------------------------
;; structure constants — generic: f^abc = -2i Tr([T^a,T^b] T^c)
;; (derivation: [T^a,T^b] = i f^abc T^c; trace both sides against T^d using
;;  Tr(T^cT^d) = 1/2 delta^cd  =>  Tr([T^a,T^b]T^d) = (i/2) f^abd)
;; ---------------------------------------------------------------------------

(defn structure-constants
  "f[a][b][cc] for a Hermitian generator set normalized Tr(T^aT^b)=1/2 delta^ab."
  [generators]
  (let [n (count generators)]
    (vec (for [a (range n)]
           (vec (for [b (range n)]
                  (vec (for [cc (range n)]
                         (let [comm (c/m-commutator (nth generators a) (nth generators b))
                               prod (c/m-mul comm (nth generators cc))
                               tr (c/m-trace prod)]
                           (c/re (c/c* (c/c 0 -2) tr)))))))))))

;; ---------------------------------------------------------------------------
;; U(1) — abelian hypercharge/electric-charge generator
;; ---------------------------------------------------------------------------

(defn u1-generators
  "A single U(1) generator = charge * identity, acting on a `dim`-component
  field of that charge. Structure constants are trivially all zero (abelian)."
  [charge dim]
  [(c/m-rscale charge (c/m-identity dim))])

;; ---------------------------------------------------------------------------
;; SU(2) — weak isospin, generators = Pauli/2 (fundamental/doublet rep)
;; ---------------------------------------------------------------------------

(def su2-generators (mapv #(c/m-rscale 0.5 %) spinor/sigmas))

;; ---------------------------------------------------------------------------
;; SU(3) — color, generators = Gell-Mann/2 (fundamental/triplet rep)
;; ---------------------------------------------------------------------------

#?(:clj (defn- sqrt3 [] (Math/sqrt 3.0))
   :cljs (defn- sqrt3 [] (.sqrt js/Math 3.0)))

(def gell-mann
  (let [ii c/i neg-i (c/c 0 -1) z c/zero o c/one neg-o (c/c -1 0)
        s (/ 1.0 (sqrt3))]
    [[[z o z] [o z z] [z z z]]                                   ; lambda-1
     [[z neg-i z] [ii z z] [z z z]]                               ; lambda-2
     [[o z z] [z neg-o z] [z z z]]                                ; lambda-3
     [[z z o] [z z z] [o z z]]                                    ; lambda-4
     [[z z neg-i] [z z z] [ii z z]]                               ; lambda-5
     [[z z z] [z z o] [z o z]]                                    ; lambda-6
     [[z z z] [z z neg-i] [z ii z]]                                ; lambda-7
     [[(c/c s) z z] [z (c/c s) z] [z z (c/c (* -2 s))]]]))         ; lambda-8

(def su3-generators (mapv #(c/m-rscale 0.5 %) gell-mann))

;; ---------------------------------------------------------------------------
;; covariant derivative
;; ---------------------------------------------------------------------------

(defn gauge-correction
  "-i g sum_a A_mu^a (T^a psi) -- the pure gauge-interaction term, i.e. the part
  of the covariant derivative beyond the ordinary partial derivative. Kept
  separate from `covariant-derivative` so `kotoba.sm.standard-model` can sum
  the corrections from all three SM factor groups against a single partial
  derivative instead of nesting three additions of it."
  [Ts A-mu g psi]
  (if (empty? Ts)
    (vec (repeat (count psi) c/zero))
    (let [zero-vec (vec (repeat (count psi) c/zero))]
      (reduce c/v-add zero-vec
              (map (fn [Ta Aa]
                     (c/v-scale (c/c 0 (- (* g Aa))) (c/m-vec Ta psi)))
                   Ts A-mu)))))

(defn covariant-derivative
  "D_mu psi = d_mu psi - i g sum_a A_mu^a (T^a psi). `d-mu-psi` is the ordinary
  partial derivative of psi along this direction (numeric or supplied
  analytically); `Ts` are the generator matrices for psi's representation;
  `A-mu` is the vector of gauge-field components A_mu^a at this spacetime
  point, one per generator; `g` is the coupling constant."
  [d-mu-psi Ts A-mu g psi]
  (c/v-add d-mu-psi (gauge-correction Ts A-mu g psi)))

;; ---------------------------------------------------------------------------
;; field strength
;; ---------------------------------------------------------------------------

(defn gauge-field-gradient
  "Given a gauge-field function A-field: R^4 -> (vector-of-4 vector-of-n reals,
  A[nu][a] = A_nu^a(x)), compute d-A[mu][nu][a] = d/dx^mu A_nu^a(x) via central
  finite differences."
  ([A-field x] (gauge-field-gradient A-field x 1e-5))
  ([A-field x h]
   (vec (for [mu (range 4)]
          (let [xp (update x mu + h)
                xm (update x mu - h)
                Ap (A-field xp)
                Am (A-field xm)]
            (vec (for [nu (range 4)]
                   (mapv (fn [ap am] (/ (- ap am) (* 2 h))) (nth Ap nu) (nth Am nu)))))))))

(defn- self-interaction-term
  "g sum_b,c f^abc A_mu^b A_nu^c, for one generator index a."
  [f-abc g A mu nu a n]
  (let [A-mu (nth A mu)
        A-nu (nth A nu)]
    (* g (reduce +
                 (for [b (range n) cc (range n)]
                   (* (get-in f-abc [a b cc]) (nth A-mu b) (nth A-nu cc)))))))

(defn field-strength
  "F_mu-nu^a = d_mu A_nu^a - d_nu A_mu^a + g f^abc A_mu^b A_nu^c.
  `f-abc` are structure constants (see `structure-constants`, all-zero for
  abelian U(1)); `d-A[mu][nu][a]` = d_mu A_nu^a; `A[mu][a]` = A_mu^a."
  [f-abc g d-A A]
  (let [n (count (first A))]
    (vec
     (for [mu (range 4)]
       (vec
        (for [nu (range 4)]
          (vec
           (for [a (range n)]
             (let [curl (- (get-in d-A [mu nu a]) (get-in d-A [nu mu a]))
                   self (self-interaction-term f-abc g A mu nu a n)]
               (+ curl self))))))))))
