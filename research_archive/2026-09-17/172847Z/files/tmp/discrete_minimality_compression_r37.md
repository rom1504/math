# Wave 37: simultaneous edge flips reduce to a row-good fractional congestion statistic

## Outcome

Full simultaneous edge-flip minimality does give a favorable-selector cover
in a high-ratio window, but it does not compress that cover.  The exact
missing statistic is a row-good **fractional covering number** defined below.
A project-scale upper bound on this number proves the restriction recurrence
immediately through (10.1022).  Conversely, a persistent normalized gap
forces this number to have logarithm `Omega(n^(3/4))`.

There is also a sharp obstruction to compressing the original full edge-cube
cover: even a fractional radius-`R` cover has total weight `exp(Omega(n))`.
Thus one cannot obtain the desired `exp(O(n^(3/4-c)))` family by pruning the
full deep-hole cover.  Any positive theorem must be selector-specific and
must prove weighted congestion after imposing the row cap.  Bare
minimality, the old nonnegative Gibbs aggregations, and counting all
low-row cuts do not provide that statistic.

No convergence proof or actual-minimizer counterexample is obtained.

## 1. Complement flips give an exact favorable-selector cover

Let `A` be an exact order-`n` minimizer, `Q(A)=q_n`, and let
`d=(sigma,x)` be an oriented state.  Write

```math
E_d=\langle A,d\rangle=q_n-\Delta_d,
\qquad
c_S(d)=\sigma x_S^TA[S]x_S.
```

For an `m`-selector `S`, flip every edge **outside** the internal edge set
`E(S)`:

```math
F_S=E(K_n)\setminus E(S).
```

The oriented energy after this simultaneous perturbation is exactly

```math
\langle A^{F_S},d\rangle
=E_d-4\sum_{e\in F_S}s_e(d)
=2c_S(d)-E_d.
\tag{D37.1}
```

Since `q_n` is the minimum cap over all signings, `Q(A^(F_S))>=q_n`.
Thus for every selector there is a state satisfying

```math
\boxed{
2c_S(d)-E_d\ge q_n
\quad\Longleftrightarrow\quad
c_S(d)\ge q_n-\frac{\Delta_d}{2}.
}
\tag{D37.2}
```

This is precisely the restriction of the full deep-hole cover (10.733) to
the structured perturbations `F_S`.

It has the useful sign missing from the internal-edge perturbation.  Recall

```math
h_d(S)=\widehat\ell(S,d)
=Q(A[S])-c_S(d)-p_2\Delta_d-B_{n,m}.
```

Every principal cap is at most `q_n`: extend a child ground by independent
uniform outside spins and average its oriented full energy.  Consequently,
if

```math
p_2=\frac{m(m-1)}{n(n-1)}\ge\frac12,
```

then every certificate in (D37.2) obeys

```math
\boxed{
h_d(S)
\le Q(A[S])-q_n+(1/2-p_2)\Delta_d-B_{n,m}
\le-B_{n,m}\le0.
}
\tag{D37.3}
```

Thus in any fixed ratio window above `1/sqrt(2)` (with the harmless finite
correction in `p_2`), simultaneous minimality covers **every** selector by
an actually favorable arbitrary-cut incidence.  The witness may depend on
`S`, and (D37.2) gives no upper bound on its row square.

The restriction to `p_2>=1/2` is explicit.  For smaller ratios the term
`(1/2-p_2)Delta_d` has the wrong sign, so complement flips alone do not give
(D37.3).  A convergence argument may work in a high-ratio active window,
but this identity is not claimed uniformly at all fixed densities.

## 2. The exact missing statistic

For a row cap `R_*`, let

```math
\mathcal D_{R_*}=\{d:R_2(d)\le R_*\}
```

and define the certificate incidence

```math
I_d=\{S:|S|=m,\ 2c_S(d)-E_d\ge q_n\}.
```

Define the row-good fractional covering number

```math
\boxed{
\tau_{\rm flip}(A,m;R_*)
=\min\left\{
\sum_{d\in\mathcal D_{R_*}}w_d:
w_d\ge0,\quad
\sum_{d:S\in I_d}w_d\ge1\ \text{for every }S
\right\},
}
\tag{D37.4}
```

with value `+infinity` if a selector has no row-good certificate.  Its exact
dual is

```math
\boxed{
\tau_{\rm flip}
=\max\left\{
\sum_S y_S:
y_S\ge0,\quad
\sum_{S\in I_d}y_S\le1
\ \text{for every }d\in\mathcal D_{R_*}
\right\}.
}
\tag{D37.5}
```

This is the sought weighted-congestion statistic.  It is stronger than
nonemptiness of every selector row and weaker than asking for a small
integral family.

There is an exact sufficient theorem.

> **Row-good flip-congestion theorem (open).**  Fix `0<c<1/4` and a
> high-ratio active window with `p_2>=1/2`.  Uniformly for every relevant
> exact-minimizer target pair,
> 
> ```math
> R_*=O(n^{9/4-c}),
> \qquad
> \log\tau_{\rm flip}(A,m;R_*)=O(n^{3/4-c}).
> \tag{D37.6}
> ```

If (D37.6) holds, average the covering constraints in (D37.4) over `U_m`:

```math
1\le\sum_dw_dU_m(I_d).
```

Hence one row-good cut has `U_m(I_d)>=1/tau_flip`.  By (D37.3), the same cut
has

```math
U_m\{h_d\le0\}\ge\tau_{\rm flip}^{-1}.
\tag{D37.7}
```

Putting `L=log tau_flip`, `t=0`, and `R=R_*` in (10.1023) proves

```math
q_m\le(m/n)^{3/2}q_n+O(n^{3/2-c}).
\tag{D37.8}
```

With the established high-ratio landing scheme, this proves convergence.
Thus (D37.6) is a precise theorem, not merely a request to count low-row
cuts.

The converse is equally exact.  If a fixed normalized gap
`G_(n,m)>=delta n^(3/2)` persists and
`R_*=O(n^(9/4-c))`, then (10.1022) bounds every row-good incidence by

```math
U_m(I_d)\le\exp\{-\Omega(n^{3/4})\}.
```

Averaging (D37.4) therefore forces

```math
\boxed{
\log\tau_{\rm flip}(A,m;R_*)=\Omega(n^{3/4}).
}
\tag{D37.9}
```

This is the row-sensitive analogue of the codebook converse (10.744).  It
pinpoints the only useful contradiction: exact minimality would have to
upper-bound (D37.4) at the strictly smaller exponent in (D37.6).

## 3. Why the full deep-hole cover cannot supply (D37.6)

Let `N=binom(n,2)` and `R=N/2-q_n/4`.  The center attached to a lifted state
is `C_d={e:s_e(d)=-1}`.  Equation (10.733) says that the radius-`R` balls
about all such centers cover the full `N`-dimensional edge cube.

Suppose even fractional weights `w_d` retain that full cover:

```math
\sum_{d:\,|F\mathbin\triangle C_d|\le R}w_d\ge1
\qquad\text{for every }F\subseteq E(K_n).
```

Averaging over a uniform edge set gives the sphere-covering lower bound

```math
\sum_dw_d
\ge\frac{2^N}{\sum_{j\le R}\binom Nj}.
\tag{D37.10}
```

Since `R=N/2-q_n/4`, Hoeffding's binomial tail yields

```math
\frac1{2^N}\sum_{j\le R}\binom Nj
\le\exp\left\{-\frac{q_n^2}{8N}\right\}.
```

Therefore

```math
\boxed{
\sum_dw_d\ge
\exp\left\{\frac{q_n^2}{8N}\right\}
=\exp\{\Omega(n)\}.
}
\tag{D37.11}
```

The last step uses the established `q_n=Omega(n^(3/2))` lower bound.
This applies a fortiori after imposing a row cap.  Hence no subfamily and no
fractional weighting of total mass `exp(O(n^(3/4-c)))` can preserve the
**full** simultaneous-flip cover.  The full edge cube contains far more
constraints than the selector family relevant to restriction.

This rules out the most direct proposed compression.  A positive argument
must throw away almost all edge perturbations and prove the selector-specific
bound (D37.6).  The wrong-way transforms (10.735)--(10.736) and
(10.754)--(10.757) do not do this: they give lower bounds on nonnegative
Gibbs averages or affinities and permit witness migration.  They supply no
upper bound on the primal mass in (D37.4), and their exact weighted-cover
duals do not control the row-restricted column congestion in (D37.5).

## 4. Finite exact diagnostics

The checker enumerates every lifted state and every selector in the
high-ratio choices `A_6,m=5`, `A_8,m=6`, and `A_9,m=7`.  It verifies
(D37.1)--(D37.3) exactly and solves (D37.4) by linear programming.

| minimizer | witness degree per selector | unrestricted `tau_flip` | selectors with no witness at cap `R_2<=n(n-1)` |
|:--|:--:|--:|--:|
| `A_6,m=5` | `12` | `2` | `0` |
| `A_8,m=6` | `8--10` | `10` | `4` |
| `A_9,m=7` | `5--15` | `55/4` | `5` |

These are finite diagnostics, not asymptotic evidence for (D37.6).  They do
show two points exactly:

1. complement-flip certificates really are favorable at `p_2>=1/2`;
2. the trace/Markov fact that many cuts have row at most the mean does not
   allow pointwise pruning of the witness cover—some selectors in `A_8` and
   `A_9` lose every certificate at the mean row cap.

The project cap `n^(9/4-c)` is asymptotically larger than `n^2`, so the
second observation is only a mechanism warning.  It is not a finite
falsifier of (D37.6).

## 5. Updated assessment of this route

Simultaneous edge minimality contains one useful new sign choice:
complement-edge flips yield the exact favorable cover (D37.3) at high ratio.
But the cover has up to `2^n` migrating cuts, and the full deep-hole cover
provably cannot be compressed below `exp(Omega(n))`, even fractionally.

The only viable continuation is to prove the selector-specific,
row-restricted congestion estimate (D37.6), probably from an additional
exact-minimizer statistic coupling the witness multiplicity to `R_2`.
Neither minimality nonemptiness, generic union bounds, nor the existing Gibbs
moments contain that statistic.  An unbounded exact-minimizer family with
`tau_flip=exp(omega(n^(3/4-c)))` for every admissible fixed `c` would falsify
this implementation, but not the arbitrary-cut route itself.

The independent checker is
`tmp/discrete_minimality_compression_r37_check.py`.
