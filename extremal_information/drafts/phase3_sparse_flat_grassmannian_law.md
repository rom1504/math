# Sparse-flat spectra and Hamming Grassmannian response entropy

**Status.** The finite counting identity and its metric-entropy consequences
below are proved.  The identity isolates the information missing between the
separated-host lower certificate and the anticode-quotient upper certificate.
It does not yet determine the binary Hamming Grassmannian exponent.

## 1. Quotient leader geometry

Let `W` be a `D`-dimensional vector space over `F_q`, equipped with a
translation-invariant metric, and let `C_0 <= W` have dimension `k`.  Give
`W/C_0` the quotient norm

```math
\|x+C_0\|_{C_0}=\min_{c\in C_0}\|x+c\|.
```

For `Delta >= 0`, define the quotient leader ball and its **sparse-flat
spectrum** by

```math
L_{C_0}(\Delta)=
 \{\bar x\in W/C_0:\|\bar x\|_{C_0}\le\Delta\},              \tag{SF.1}
```

```math
\Lambda_{C_0}(\ell,\Delta)=
 \#\{U\le W/C_0:\dim U=\ell,\ U\subseteq L_{C_0}(\Delta)\}. \tag{SF.2}
```

Thus `Lambda` does not store labelled coset leaders or the full distance
landscape.  It stores how many linear flats of each dimension survive a
leader threshold.  Its one-dimensional term recovers leader-ball cardinality
(up to the zero vector and scalar orbits), while its higher terms retain
incidence information that cardinality discards.  It is an unlabelled
aggregate of the quotient norm; no universal strict-compression claim is
being assumed.

For two subsets write the directed Hausdorff distance

```math
h^\to(C,C_0)=\max_{c\in C}d(c,C_0).                           \tag{SF.3}
```

The ordinary Hausdorff distance is the maximum of the two directions.

For a subspace `C_0 <= S`, write

```math
\rho_S(C_0)=\max_{s\in S}d(s,C_0)
```

for its covering radius relative to `S`.  Translation invariance gives the
useful exact identification

```math
h^\to(C,C_0)=\rho_{C+C_0}(C_0),\qquad
d_H(C,C_0)=
 \max\{\rho_{C+C_0}(C_0),\rho_{C+C_0}(C)\}.                  \tag{SF.3a}
```

Indeed, writing `s=c+c_0` shows `d(s,C_0)=d(c,C_0)`.  Thus the new
Grassmannian metric is a symmetric pair of **relative covering radii**, not
the usual subspace/injection distance.

## 2. Exact directed Grassmannian balls

### Theorem SF.1 (sparse-flat ball identity)

For every `Delta >= 0`,

```math
\#\{C\in\operatorname{Gr}_k(W):h^\to(C,C_0)\le\Delta\}
=\sum_{\ell=0}^{\min\{k,D-k\}}
 {k\brack\ell}_q q^{\ell^2}
 \Lambda_{C_0}(\ell,\Delta).                                \tag{SF.4}
```

#### Proof

Let `pi:W -> W/C_0`.  The directed-distance condition is exactly
`pi(C) subseteq L_(C_0)(Delta)`.  Put

```math
\ell=\dim\pi(C)=k-\dim(C\cap C_0).
```

Choose `J=C cap C_0`, an `(k-ell)`-subspace of `C_0`, in
`{k bracket ell}_q` ways, and choose
`U=pi(C)`, counted by `Lambda_(C_0)(ell,Delta)`.  After fixing a linear
section of `pi` over `U`, the subspace `C/J` is the graph of an arbitrary
linear map

```math
U\longrightarrow C_0/J.
```

Both spaces have dimension `ell`, so there are `q^(ell^2)` graphs.  The
intersection, image and graph are recovered uniquely from `C`, and summing
over `ell` proves (SF.4). `square`

This is not ordinary Grassmannian injection distance: all metric dependence
is concentrated in which quotient flats lie inside the coset-leader ball.

### Checks at the two extremes

If the covering radius of `C_0` is at most `Delta`, then `L_(C_0)(Delta)` is
the whole quotient and

```math
\Lambda_{C_0}(\ell,\Delta)={D-k\brack\ell}_q.
```

Formula (SF.4) becomes the standard intersection decomposition of
`{D bracket k}_q`; every `k`-subspace is directed-close to `C_0`.

If `L_(C_0)(Delta)` is itself a `t`-dimensional subspace, then
`Lambda_(C_0)(ell,Delta)={t bracket ell}_q`.  This is exactly the algebraic
closure exhibited by a strict metric synchronization: the entire directed
ball depends only on the coarse quotient dimension `t`.

### Proposition SF.2 (the one-sided spectrum is not a response state)

There are two scalable binary Hamming carriers whose quotient normed spaces,
and hence sparse-flat spectra at every threshold, are isometric, but whose
Hausdorff response to the zero carrier differs by a linear amount.

In one four-coordinate block put

```math
C^{(2)}=\operatorname{span}(1100),\qquad
C^{(1)}=\operatorname{span}(1000).
```

Both quotients are isometric to the three-dimensional Hamming cube.  For the
first quotient, the nonzero coordinate records the parity class of the first
two bits and has leader cost one; for the second, quotienting simply deletes
the first bit.  Take `r` direct sums.  The quotient normed spaces are both
the `3r`-dimensional Hamming cube, so all their functions `Lambda(ell,Delta)`
agree.  Nevertheless

```math
d_H((C^{(2)})^{\oplus r},\{0\})=2r,
\qquad
d_H((C^{(1)})^{\oplus r},\{0\})=r.             \tag{SF.4a}
```

Thus even the full unlabelled quotient norm forgets macroscopic **rooted
lift geometry**: the metric carried by the kernel itself.  Sparse-flat data
is an exact state for the directed query *toward its own center*, but cannot
be promoted to a symmetric all-future response state without an additional
rooted datum.

Already at `D=4`, `k=1`, and `Delta=1`, the centers `span(0011)` and
`span(0010)` have the same quotient leader-ball size and identical
sparse-flat spectra, while their symmetric Grassmannian balls contain five
and seven lines respectively.  This finite witness separates the two-sided
failure from the zero-carrier comparison.

## 3. A universal packing consequence

Let

```math
V_W(\Delta)=\#\{x\in W:\|x\|\le\Delta\}.
```

Every quotient leader has a representative in the ambient ball, hence

```math
|L_{C_0}(\Delta)|\le
\min\{q^{D-k},V_W(\Delta)\}.                                  \tag{SF.5}
```

Counting ordered bases inside the leader ball gives, for `ell >= 1`,

```math
\Lambda_{C_0}(\ell,\Delta)
\le {|L_{C_0}(\Delta)|^\ell\over |\operatorname{GL}(\ell,q)|}. \tag{SF.6}
```

Since `q^(ell^2)/|GL(ell,q)| < 4` and
`{k bracket ell}_q <= 4q^(ell(k-ell))`, Theorem SF.1 implies

```math
B^\to_{C_0}(\Delta)
\le 1+16\sum_{\ell=1}^{\min\{k,D-k\}}
 q^{\ell(k-\ell+\lambda)},                                   \tag{SF.7}
```

where

```math
\lambda=\log_q\min\{q^{D-k},V_W(\Delta)\}.                  \tag{SF.8}
```

The symmetric Hausdorff ball is contained in this directed ball.  Greedy
deletion therefore proves

```math
\log_q\operatorname{Pack}(\operatorname{Gr}_k(W),d_H,\Delta)
\ge
\log_q{D\brack k}_q
-\max_{0\le\ell\le\min\{k,D-k\}}
  \ell(k-\ell+\lambda)
-O(\log D).                                                    \tag{SF.9}
```

Here packing means pairwise distance strictly greater than `Delta`.

### Corollary SF.3 (binary Hamming exponent)

Let `W=F_2^D`, `k/D -> kappa`, and `Delta/D -> delta`, with
`0<kappa<=1/2` and `0<delta<1/2`.  Put

```math
\lambda_*=\min\{1-\kappa,H_2(\delta)\}.
```

Then

```math
\liminf_{D\to\infty}{1\over D^2}\log_2
 \operatorname{Pack}(\operatorname{Gr}_k(F_2^D),d_H,\Delta)
\ge
\kappa(1-\kappa)
-\max_{0\le\eta\le\kappa}
 \eta(\kappa-\eta+\lambda_*).                               \tag{SF.10}
```

This follows from (SF.9), the Hamming-ball entropy estimate, and
`log_2 {D bracket k}_2=k(D-k)+O(1)`.

When `kappa <= H_2(delta) < 1-kappa`, the right-hand side is

```math
\kappa(1-\kappa-H_2(\delta)),                                \tag{SF.11}
```

the exponent obtained by placing the carriers in a Gilbert--Varshamov
separated host of rate `1-H_2(delta)`.  When
`H_2(delta)<kappa`, (SF.10) is smaller than that host exponent by
`(kappa-H_2(delta))^2/4`.  Therefore the universal counting argument that
uses only leader-ball cardinality cannot improve the common-host
construction.  A sharper directed-ball estimate must exploit the
linear-flat spectrum (SF.2), not just the number of short syndromes; and
Proposition SF.2 shows that even this does not by itself control symmetric
response geometry.

## 4. The coarse quotient upper mechanism

Puncture any `t` coordinates.  If two carriers have the same projected
subspace, every vector in either carrier has a vector in the other with the
same projection, so their Hausdorff distance is at most `t`.  Every packing
at distance strictly greater than `t` therefore injects into the subspaces of
the `(D-t)`-dimensional quotient.  In particular, when `k<=(D-t)/2`,

```math
\log_2\operatorname{Pack}
 (\operatorname{Gr}_k(F_2^D),d_H,t)
\le k(D-k-t)+O(D).                                             \tag{SF.12}
```

The count includes projected dimensions below `k`; their sum changes only
the lower-order term in the displayed regime.  This is the
anticode/puncturing upper scale from Theorem 13.3.  The interval between
(SF.10) and (SF.12) is the interval between the sparse-flat content of
quotient leader balls and a common-support synchronization quotient.

## 5. Composition and the remaining theorem

For an `ell_1` direct sum `(W_1,C_1) direct-sum (W_2,C_2)`, quotient norms
add.  Product flats give the rigorous supermultiplicative relation

```math
\Lambda_{C_1\oplus C_2}(\ell,\Delta)
\ge
\max_{\substack{\ell_1+\ell_2=\ell\\
                 \Delta_1+\Delta_2\le\Delta}}
 \Lambda_{C_1}(\ell_1,\Delta_1)
 \Lambda_{C_2}(\ell_2,\Delta_2).                             \tag{SF.13}
```

Mixed flats can make the inequality strict.  Thus the sparse-flat spectrum
has a genuine composition operation, but product data need not close it.

The next theorem is now precise: obtain uniform exponential bounds on
`Lambda_C(ell,delta D)` for binary linear codes `C` of rate `kappa`, or show
that two code families with the same leader-ball entropy have macroscopically
different sparse-flat spectra.  Either result would decide whether this
middle invariant compresses the response geometry or merely relocates its
complexity.
