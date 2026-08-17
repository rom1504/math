# Independent audit: bounded-cap optimizer fibres and witness covers

**Verdict:** PASS.  The orbit count, probability estimate, cap accounting,
packing constants, child-cap averaging argument, and repaired projective
denominator are all valid.  Several constants are deliberately conservative
but none is used in the wrong direction.  I ran
`experiments/verify_bounded_cap_optimizer_switching.py`; all exhaustive and
randomized finite checks pass.

## 1. WS.1: switching orbits and common witnesses

The switching group has size `2^(k-1)`.  Its action on complete signings is
free: if `s_i s_j A_ij=A_ij` on every edge, then all `s_i` are equal, which
is the identity modulo global sign.  Hence the exact orbit count is

```math
2^{\binom{k}{2}}/2^{k-1}=2^{\binom{k}{2}-k+1}.
```

Every orbit has a representative maximized at `1`: choose a top Boolean
state and switch by it.  A further switch by the prescribed `u` gives a
representative maximized at both `u` and `-u`.

For a fixed spin, `H_A(x)` is a sum of `E=binom(k,2)` independent signs.
At `t=2k^(3/2)`, Hoeffding gives

```math
2\exp[-t^2/(2E)]
=2\exp[-4k^2/(k-1)]
\le2e^{-4k}.
```

The union bound over `2^k` spins is therefore exactly bounded by the
quantity in WS.7, which is below one half for `k>=3`.  Since cap is switching
invariant and all orbits have equal size, more than half the orbits are good.
Selecting one prescribed-optimizer representative from each gives at least
`2^(E-k)` distinct children.

Both `H_A` and `|u dot x|` attain their maxima at `+-u`, so the common-witness
claim follows.  The appended graph is complete on `k+1` vertices, and

```math
|H_A(x)+y u\cdot x|\le Q(A)+k
```

proves the parent-cap claim.  There is no missing factor from the quadratic
normalization.

## 2. WS.2: one scalar coordinate

Inside a common-witness fibre,

```math
R_g(A)=H_A(u)+g(u).
```

The completed-parent cap places this scalar in an interval of length
`2C Lambda^(3/2)k^(3/2)`.  A one-dimensional packing with spacing
`epsilon k^(3/2)` has at most

```math
1+floor(2C Lambda^(3/2)/epsilon)
```

points.  This checks WS.10 exactly.

## 3. WS.3: approximate dictionaries

For a common dictionary `U`, max nonexpansiveness gives

```math
|R_A^U(g)-R_(A')^U(g)|
\le\max_{u\in U}|H_A(u)-H_(A')(u)|.
```

Each truncation error lies in `[0,tau k^(3/2)]`.  The draft pays two such
errors by the triangle inequality, yielding the valid, though nonsharp,
separation `(epsilon-2tau)k^(3/2)`.  (Because both errors are one-sided,
`epsilon-tau` can in fact be obtained, but no conclusion relies on this
improvement.)

Every evaluation coordinate lies in an interval of length
`2C_0 k^(3/2)`.  Dividing it into

```math
1+ceil(2C_0/(epsilon-2tau))
```

half-open bins makes each bin strictly narrower than the packing distance.
The product count gives WS.16 and logarithms give WS.17.  Thus an
`exp(alpha k)` response packing really forces a linear-size reusable witness
dictionary.

The context--switching corollary is a direct application: the union of all
query witness sets is an exact cover, its size is at most their total size,
and pigeonhole gives the two stated alternatives.

## 4. Child cap from parent cap

For a hollow quadratic parent,

```math
P_A(x,y)=H_A(x)+x^TBy+H_C(y),
```

uniform averaging over Boolean `y` kills both the bilinear term and every
hollow quadratic monomial in `H_C`.  Hence

```math
H_A(x)=E_y P_A(x,y).
```

Pointwise boundedness by `Q(P_A)` then gives
`|H_A(x)|<=Q(P_A)` and therefore `Q(A)<=Q(P_A)`.  The substitution
`C_0=C Lambda^(3/2)` is correct.  This argument would need adjustment in the
presence of a nonzero constant calibration or diagonal term; neither is
present in the stated parent class.

## 5. Projective variant and WS.23

Let `e_A(g)=R_A(g)-R_A^U(g)`, so `e_A in [0,tau k^(3/2)]`.  The difference
`e_A-e_(A')` has oscillation at most `2tau k^(3/2)`.  Therefore projective
distance loses at most `tau k^(3/2)` under truncation.  If

```math
d_proj(A,A')>=epsilon k^(3/2),
```

the truncated response difference has projective distance at least
`(epsilon-tau)k^(3/2)`.

Writing `d_u=H_A(u)-H_(A')(u)`, every truncated response difference lies
between `min_U d_u` and `max_U d_u`.  Hence its oscillation is at most the
range of `(d_u)`.  After anchoring one fixed `u_0`, at least one coordinate
`d_u-d_(u_0)` has magnitude at least half that range.  This actually gives
an anchored separation of at least `(epsilon-tau)k^(3/2)`.

The repaired draft uses the weaker denominator `epsilon-2tau`, so WS.23 is
certainly valid under `2tau<epsilon`.  Anchored coordinates lie in
`[-2C_0k^(3/2),2C_0k^(3/2)]`, an interval of length `4C_0k^(3/2)`, and there
are `|U|-1` of them.  This verifies both the numerator `4C_0`, the exponent,
and the repaired denominator.

For an absolute outer response, the same proof applies after treating the
sign channel as part of the declared future and using signed evaluations
`sH_A(u)`.  This qualification is implicit in the draft's “signed witness”
language.

## 6. Scope

The results do not bound the raw number of children in an optimizer fibre;
WS.1 proves that such a bound is false.  They bound only response-distinct
children per reusable approximate witness dictionary.  A linear dictionary
is fully compatible with the linear hidden rate of the flat Gram family, so
WS.3 is a strict narrowing rather than an exact-sign compiler impossibility.
It also does not decide whether witness identity is child-dependent: an
extensive collection of common query pins remains an allowed realization.

These qualifications are explicit in the draft.  Subject to them, the note
is rigorous and ready for canonical use.
