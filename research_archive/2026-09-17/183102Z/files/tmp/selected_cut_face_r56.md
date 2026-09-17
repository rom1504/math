# Wave 56: selected nonuniform cut laws and the exact collision obligation

Status: exact finite-law statements below are verified algebraically.  They
do **not** prove the missing active-face incidence bound.  No uniform local
state or uniform outside-completion law is used.

## 1. Bare incidence in parent-deficit coordinates

Fix an exact minimizer `A` of order `n`, put `q=Q(A)=q_n`, and write

```
p2=(m)_2/(n)_2,                 B=(p^(3/2)-p2)q.
```

For an oriented global cut `d=(sigma,x)`, let

```
E(d)=sigma x^T A x,             Delta(d)=q-E(d)>=0,
R(d)=x^T A^2 x,
delta_S(d)=Q(A[S])-sigma x_S^T A[S]x_S>=0.
```

Direct substitution in (10.792) gives the exact identity

```
hat ell(S,d)=delta_S(d)-p2 Delta(d)-B.                 (R56.1)
```

Thus, for a tolerance `t`, define

```
F_t={(S,d): delta_S(d)<=B+t+p2 Delta(d)},
u_t(d)=U_m{S:(S,d) in F_t}.                            (R56.2)
```

This is exactly the bare favorable incidence in (10.795), not the stronger
complement incidence.  The negative parent-deficit term is retained.

## 2. Captured extraction for an arbitrarily selected global-cut prior

Let `nu` be **any** selector-independent probability law on global cuts.  It
may be singular, supported on the parent ground face, supported on a
near-ground layer, row tilted, or chosen by a deterministic global
optimization.  Put

```
Z_nu = E_{d~nu} u_t(d).
```

Assume `Z_nu>0` and let `mu` be the captured output law

```
mu(d)=nu(d)u_t(d)/Z_nu.
```

Write `Rbar=E_mu R` and `Dbar=E_mu Delta`.  For arbitrary `a,b>1` satisfying

```
theta=1-1/a-1/b>0,
```

there is a deterministic cut `d` such that

```
u_t(d)>=theta Z_nu,      R(d)<=a Rbar,      Delta(d)<=b Dbar.  (R56.3)
```

Proof: under `mu`, Markov and the union bound give

```
mu{R<=a Rbar, Delta<=b Dbar}>=theta.
```

Consequently the `nu`-weighted favorable mass of this good-cost set is at
least `theta Z_nu`.  Since its total `nu` mass is at most one, one member has
`u_t(d)>=theta Z_nu`.  It has the two displayed cost bounds by membership.

If `nu` is supported on the exact parent ground face, `Delta=0` identically
and the deficit clause can be omitted.  For every `a>1`, one then obtains

```
u_t(d)>=(1-1/a)Z_nu,             R(d)<=a Rbar.          (R56.4)
```

In particular `a=2` loses only a factor two.  With both costs, `a=b=3`
gives the convenient factor-three statement

```
u_t(d)>=Z_nu/3,    R(d)<=3Rbar,    Delta(d)<=3Dbar.     (R56.5)
```

The constants trade continuously according to (R56.3); there is no need to
bury the tradeoff in big-O notation.

The entropy cost is exact.  Condition the product law `U_m tensor nu` on
`F_t` and call the result `P`.  The chain rule gives

```
-log Z_nu
 = I_P(S;D)+D(P_S||U_m)+D(P_D||nu).                    (R56.6)
```

Thus a selected-prior package

```
-log Z_nu=O(H),             Rbar=O(n^(9/4-c))           (R56.7)
```

already produces the fixed cut required by (10.795).  If desired, a
captured near-ground condition `Dbar=O(n^(3/2-c))` survives extraction too,
but the bare lemma itself does not require it.

The important research obligation is the first clause of (R56.7): one must
construct a selector-independent active-face/near-ground law whose *captured
incidence* is saved.  (R56.3) does not establish this clause.

## 3. Exact arbitrary-channel collapse

There is a reference-free form which makes the collision obligation
explicit.  Let `P` be any joint selector-cut law supported on `F_t`, and put

```
K(P)=I_P(S;D)+D(P_S||U_m).
```

For every output `d`, the posterior is supported on a set of uniform mass
`u_t(d)`.  KL projection onto a set therefore gives

```
D(P_{S|d}||U_m)
 =D(P_{S|d}||U_m(.|F_t^d))+log(1/u_t(d))
 >=log(1/u_t(d)).
```

Averaging and using the chain rule proves

```
K(P)>=E_D log(1/u_t(D)).                              (R56.8)
```

More precisely, for every `lambda,eta>=0`,

```
inf_{P:F_t a.s.} {K(P)+lambda E R(D)+eta E Delta(D)}
 = min_{d:u_t(d)>0}
   {log(1/u_t(d))+lambda R(d)+eta Delta(d)}.           (R56.9)
```

The lower bound follows from (R56.8) and averaging.  Equality takes `D=d`
constant and `S` distributed as `U_m` conditioned on `F_t^d`.

As a quantitative corollary, if

```
K(P)<=K,       E R(D)<=R0,       E Delta(D)<=D0,
```

then, when all three positive budgets are retained, choose
`lambda=K/R0`, `eta=K/D0` in (R56.9).  Some deterministic output satisfies

```
u_t(d)>=exp(-3K),       R(d)<=3R0,       Delta(d)<=3D0. (R56.10)
```

If the law is ground-face-supported, omit the deficit term and improve 3 to
2.  Formula (R56.3) is much sharper whenever the channel actually comes from
conditioning a selected product prior: it preserves `Z_nu` up to a constant
rather than replacing it by a power.

## 4. Why this does not trigger the Wave 55 localization theorem

The proof of (10.1315)--(10.1318) uses a uniform oriented local state and a
uniform independent outside completion.  Under that law, the completed word
is the uniform full Rademacher word and its quadratic completion increment
has the Hanson--Wright law in (10.1315).

Here `nu` is an arbitrary globally selected law and can be supported on an
exponentially singular parent face or even on one cut.  There is no local
state/outside-completion factorization and no reason for the completed word
to be Rademacher.  Therefore (10.1315) cannot be invoked.  This is exactly
the surviving distinction noted after (10.1318): the conclusion is a
selected fixed column.

## 5. Sharp scoped obstruction: selector-wise active witnesses are useless

The following abstract incidence model shows that pointwise existence,
exact parent activity, and perfect row control alone cannot prove (R56.7).
Let `X` be the `N=binom(n,m)` selectors and let the cut labels be
`D={d_S:S in X}`.  Declare

```
F={(S,d_S):S in X},       Delta(d_S)=0,       R(d_S)=R0.
```

Thus every selector has a favorable exact-active, row-`R0` witness.  Yet

```
u(d_S)=1/N=exp{-Theta(n)}
```

at fixed density, so no column has entropy `O(H)` for `H=o(n)`.  More
sharply, every joint law supported on this matching obeys

```
I(S;D)+D(P_S||U_m)=log N.                              (R56.11)
```

Indeed `D` is a bijective deterministic image of `S`, so the two terms are
`H(P_S)` and `log N-H(P_S)`.  Equality holds in (R56.8).  This model is an
incidence-level no-go, not a claimed realization by a complete signing.  It
proves that any actual signing theorem must add a multiplicity/overlap or
low-information conclusion; the words "active", "near ground", and "low
row" do not supply one by themselves.

## 6. Audit against Wave 55 migration

The migration theorem (10.1326) does not currently bound the new collision
quantity.

1. For each high-field child incidence it produces *some* parent state
   reversing a fixed fraction of an extracted edge block, but it gives no
   bound on the local deficit `delta_S(omega)`.  Hence the witness is not
   known to lie in `F_t`; reversal of child-positive internal edges can in
   fact increase that local deficit unless other edges compensate.
2. Even if favorability were added, (10.1326) gives neither a bound on the
   entropy of the witness map nor a lower bound on the number of child
   fibres sharing one witness.  The matching model shows that existential
   witnesses can have `K=Theta(n)`.
3. The finite failure of pointwise common-witness compatibility recorded
   after (10.1326) is consistent with both gaps.

Even granting a stronger minimax upgrade of migration--namely, one
selector-independent law `mu` on `Delta=o(T_n)` states for which every
selector's extracted excess block is escaped with constant `mu`-probability--
only the entropy of the **escape** incidence would be `O(1)`.  To invoke
(R56.3) one still needs either

```
escape(S,d) => (S,d) in F_t
```

or a saved lower bound on their intersection, plus the captured estimate
`E[R(D)|F_t]=O(n^(9/4-c))`.  Neither follows from (10.1326): escape controls
one signed internal edge sum, whereas `F_t` controls the complete local
deficit in (R56.1), and a generic near-ground state has only the much weaker
row bound `R=O(n^(5/2))`.  Thus the hypothesized minimax law would solve the
output-overlap problem for its own block event, but not yet the two required
favorability and project-row clauses for the bare fixed-cut tail.

Accordingly the exact bridge from migration to this route would have to say:
on a selector set of mass `exp{-O(H)}`, choose migrated witnesses which are
bare-favorable and whose joint law has `K=O(H)` (or, more directly, construct
`nu` satisfying (R56.7)).  Neither clause follows from Wave 55.

## 7. Resulting scoped judgment

(R56.3)--(R56.9) give a clean selected nonuniform-law interface and preserve
the full parent-deficit term.  They rigorously avoid the retired uniform
completion mechanism.  They also show exactly why this is not yet a proof:
the genuinely new mathematical content must be a saved active-face or
near-ground incidence `Z_nu`, or an equivalent `K=O(H)` overlap theorem.
Pointwise box/migration witnesses without multiplicity cannot provide it.
