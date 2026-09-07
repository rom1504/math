# Wave 51 main audit: a hypercontractive completion boundary layer

Status: independently derived candidate theorem, pending cross-audit against
the agent reports and the ledger definitions.

Fix the local state in (10.1243).  As a function of the outside Rademacher
spins, the centered completion increment `Z` is a multilinear polynomial of
degree at most two, with `E Z=0` and `E Z^2=V`.  If `V>0`, Bonami--Beckner at
`2 -> 4` gives

    ||Z||_4 <= 3 ||Z||_2,

and interpolation (`||Z||_2 <= ||Z||_1^(1/3)||Z||_4^(2/3)`) gives

    E|Z| >= (E Z^2)^(3/2)/(E Z^4)^(1/2) >= sqrt(V)/9.

Writing `Z_- = max(-Z,0)`, centering yields

    E Z_- = E|Z|/2 >= sqrt(V)/18.

For every `r>=0`, split at `r` and apply Cauchy--Schwarz:

    E Z_- <= r + E[Z_- 1_{Z_->r}]
            <= r + sqrt(V P{Z<-r}).

Therefore the following is verified provided the standard degree-two
Bonami--Beckner normalization above is the one used:

    P{Z<=-r} >= (1/18-r/sqrt(V))_+^2.                 (A)

In particular `P{Z<=0}>=1/324`.  If `V=0`, then `Z=0` almost surely, so the
nonnegative-threshold CDF is one and every strictly negative-threshold CDF is
zero.

In the notation `g=(1-p_2)e-h` of (10.1244), define

    B_t(S,a) = 1/324,                                 g>=0, V>0;
               (1/18+g/(p_2 sqrt(V)))_+^2,            g<0, V>0;
               1,                                     g>=0, V=0;
               0,                                     g<0, V=0.

Then direct substitution into (A) proves

    Z_t >= E_{S,a} max{C_t(S,a), B_t(S,a)},            (B)

where `C_t` is the two-support/two-moment envelope from (10.1245)--(10.1246).
The maximum is valid pointwise because both expressions lower-bound the same
conditional CDF.

Consequences and scope:

- Every positive-margin local state contributes a fixed constant, even when
  `g^2/V` is arbitrarily small.  Thus the ratio loss in the positive Cantelli
  branch is unnecessary.
- A new negative-margin mechanism is available whenever
  `-g < p_2 sqrt(V)/18`; it does not require small Parseval defect `D`.
- It remains open to prove saved mass of the enlarged boundary-layer event
  `g >= -p_2 sqrt(V)/36` (the latter would contribute at least `1/1296` per
  state).  This abundance statement, uniformly at the target exponent, would
  prove (10.795).
- The lemma alone does not control negative margins much larger than one
  conditional standard deviation and proves no asymptotic abundance.

Potential audit points: confirm that the fixed-state completion increment is
indeed degree at most two after orientation is frozen; ensure `p_2>0`; decide
whether to use `<` or `<=` only in the harmless direction shown above.
