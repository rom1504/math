# Independent audit: deep-hole projective packing frontier

**Verdict: PASS after repairs.**

This audit checks the frozen version of
`deep_hole_projective_packing_frontier.md` with SHA-256
`01e1080a6b5efaceb73ceabea776ee96378947d736c781a5dde60fb7ce189532`.
The first draft had three scope errors and one false finite-order clause in
PP.4.  All theorem-level errors were repaired before this verdict.  PP.1--PP.4
now prove exactly the advertised limited conclusions; they do not prove
Cut-DH(3), a fixed-scale packing, or convergence.

## 1. Normalization and PP.1

For the augmented cut code `C_n` of length `N=binom(n,2)`, if `y` encodes
the signing, then

```math
\langle a,z_c\rangle=N-2d(y,c),
\qquad
Q(a)=N-2d(y,C_n).
```

Thus the covering radius is `(N-M_n)/2`, and code radius excess `s` is
positive-energy deficit `2s` whenever the shell is thin enough to remain
positive.  The repaired projective weight formula is

```math
\min\{k(n-k),N-k(n-k)\};
```

its nonzero minimum is `n-1` for `n>=5`.

PP.1 is correct.  For every coordinate `e`, a nearest codeword to `y+e`
must disagree with `y` at `e`; otherwise it would lie strictly closer than
the deep-hole radius to `y`.  Its error support has size at most `rho+1`.
These supports cover all `N` coordinates, so at least
`ceil(N/(rho+1))` shell words are needed.  If `1 in C` and
`rho+1<N/2`, the shell contains at most one word from each antipodal pair.
For `C_n`, `M_n>2` therefore gives three projective words and AO separation
parameter `2(n-1)/N=4/n`.

The conditional three-subcover inequality is also exact.  If
`N_1 union N_2 union N_3=[N]`, then for every pair `{i,j}`, with `k` the
remaining index,

```math
N\le |N_i|+|N_i\mathbin\triangle N_j|+|N_k|,
```

while both `d(c_i,c_j)` and its complement are at least
`N-2(rho+1)`.  Hence the projective distance is at least `M_n-2` for
the augmented cut code.  This is only `Theta(n^(-1/2))` in normalized AO
scale, so the claim that AO.2's `O(n^(5/4))` simultaneous error overwhelms
its `Theta(n)` signal is correct.

## 2. PP.2: the generic antipodal countermodel

The replacement of the earlier Hamming-factor sketch by

```math
C_m=\operatorname{Rep}_{m^4-m^3}\oplus\mathbb F_2^{m^3}
```

is valid and cleaner.  Direct-sum covering radii add.  A balanced received
word on the repetition block is a deep hole of radius

```math
(m^4-m^3)/2=N_m/2-N_m^{3/4}/2,
```

and the dimension is `m^3+1=Theta(N_m^(3/4))`.  The shell identity

```math
\mathcal L_s(y)=
\{(r,b):r\in\{0,1\},\ d(b,y_B)\le s\}
```

is exact.  Same-repetition pairs have distance at most `2s`; opposite pairs
have complementary distance at most `D_m=m^3`.  Consequently the whole
projective diameter is at most `max(D_m,2s)=o(N_m)` for every
`s=o(N_m^(3/4))`.

The stronger flip-certificate paragraph is also correct.  For an `r`-set
`F`, take `T=F cap H`, `h=|T|`, and whichever of `U union T` and `V union T`
contains more of `F`.  Its intersection has size at least `(r+h)/2`, exactly
the amount required for

```math
d(y+F,c)\le rho.
```

Thus this is a scalable counterexample to any argument using only antipodal
linearity, deepest-coset stability, the `N^(3/4)` deficit, and all finite
flip-cover constraints.  It is not a counterexample with the augmented cut
code's `Theta(sqrt N)` dimension or cut-weight algebra.

## 3. PP.3: sphere-covering barrier

With the repaired assumptions `k>=1` and `t>=0`, PP.3 is correct.  Covering
the `D`-cube gives

```math
1\le 2^k\Pr\{\operatorname{Bin}(D,1/2)\le D/2-t\},
```

and Hoeffding gives the upper bound `exp(-2t^2/D)`.  Therefore

```math
D\ge {2t^2\over k\log 2}.
```

For `t=Theta(n^(3/2))` and `k=O(n)`, this forces `D=Omega(n^2)`.
The repaired text correctly draws only the scoped conclusion: PP.2's
particular `D=o(N)` diameter certificate is unavailable.  The inequality
does not rule out a different localization mechanism inside a linear-size
`B` block or a nonseparable augmented-cut construction.

## 4. PP.4: cut-specific agreement core

The algebraic part, PP.15--PP.16, is correct without an extra cut-shape
hypothesis.  If the actual distance is `N-D`, then the two error supports
have symmetric difference `N-D` and agreement size `D`, so

```math
|B|={w_0+w_1-N+D\over2},
\quad
|Z|={N-w_0-w_1+D\over2},
\quad
|Z|-|B|=N-w_0-w_1\ge M_n-2s.
```

In particular `D>=M_n-2s`.

The original PP.4 incorrectly inferred from only `D<N/2` that an actual
large lift is the complement of a cut.  A balanced cut can itself have
weight above `N/2`; exact active-shell counterexamples already occur for
`n=5,6`.  The repaired hypothesis

```math
D<N-\lfloor n^2/4\rfloor
```

fixes this: `N-D` then exceeds the maximum cut size, so the difference must
be a complemented cut.  The agreement coordinates are exactly
`delta(S)`, and after switching the two signed words are `1` and
`-c(v_S)`.  Adding their energy lower bounds yields

```math
\sum_{e\in\delta(S)}a_e\ge M_n-2s.
```

The hypothesis is automatic when `D=o(N)`.  Then
`D=|S|(n-|S|)=o(n^2)` allows the smaller shore convention
`|S|=o(n)`, while `D>=M_n-2s` gives
`|S|>= (M_n-2s)/n=Omega(sqrt n)` for `s=o(n^(3/2))`.
This is a genuine cut-specific description of the surviving collapsed
two-cap case, not a contradiction to it.

## 5. Independent finite checks

I independently enumerated linear codes generated by up to three vectors at
lengths `N<=6`: PP.1 and its support-cover proof passed `37,974` deep-hole
cases.  For the PP.2 seed `m=2`, exact shells of excess `s=0,1,2,3` had
sizes `2,18,74,186` and projective diameter `8`, matching PP.9--PP.10.
PP.3 passed all enumerated linear codes of lengths at most six satisfying
`k>=1,t>=0`.

For augmented cut codes at `n=3,4,5,6`, I enumerated all deep holes and all
shell pairs of excess at most two.  The repaired PP.15--PP.19 passed in every
applicable case (31,224 pairs met the complemented-cut threshold).  The same
enumeration exposed the now-repaired balanced-cut gap in the earlier PP.4.
No theorem in the draft depends on these computations.

## 6. Archive comparison and scope

* PP.1 is the `r=1` specialization of the archived exact flip inequality in
  `nearmin_deterministic_inequalities.md`; its explicit support cover,
  cardinality floor, and projective three-point consequences are the new
  deductions.  It should not be presented as a new independent source of
  exact-minimality information.
* PP.2's direct-sum mechanism is elementary, but no exact scalable generic
  two-cap deep-hole countermodel with this shell statement was found in the
  archive.  It strengthens the earlier synthetic two-cap warning by making
  all coding and flip-cover conditions physical in a genuine code.
* PP.3 is the classical sphere-covering plus Hoeffding argument.  Its new
  contribution here is only to delimit PP.2; it is not a general
  low-dimensional no-localization theorem.
* PP.4 refines the archived scalar positivity inequality AO.20.  AO.20
  already forced projective separation at the `Theta(n^(-1/2))` scale, but
  did not identify the `K_(|S|,n-|S|)` agreement core or its exact bias.
* Cut-DH(3) is exactly a finite-size version of the current
  `L_projective` obligation.  The implication through audited AO.2 has the
  correct factor `d=2s`, projective normalization `d_pr>=gamma N/2`, order
  `n+O(sqrt n)`, cap `O(n^(3/2))`, and response gap `Theta(n^(3/2))` for
  fixed parameters.  It yields contextual state bits for the declared query
  bank, not a cross-order congruence or convergence theorem.

The final draft also repairs two nonblocking scope phrases: its status now
includes PP.4, and the executive conclusion says that a successful argument
must exploit at least the low dimension or the special cut algebra rather
than claiming both have separately been proved necessary.

**Final classification:** theorem-level **PASS after repair**; generic-code
falsifier **PASS**; cut-specific structural reduction **open exactly where
stated**.
