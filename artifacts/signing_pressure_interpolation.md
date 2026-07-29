# Signing-space pressure: interpolation audit

## 0. Definitions and why this pressure would solve convergence

Let

\[
M(A)=
\max_{x\in\{\pm1\}^n}
\left|\sum_{i<j}a_{ij}x_ix_j\right|,
\qquad
F(n)=\min_A M(A),
\]

and let \(E_n=\binom n2\).  It is useful to distinguish the raw
Laplace parameter from the scale-normalized one:

\[
\mathscr Z_n(\lambda)
=\sum_{A\in\{\pm1\}^{E_n}}e^{-\lambda M(A)},
\]

\[
Z_n(\beta)=\mathscr Z_n(\beta\sqrt n),
\qquad
\Phi_n(\beta)=\frac1{n^2}\log Z_n(\beta).
                                                               \tag{0.1}
\]

The factor \(\sqrt n\) is essential: \(M(A)=\Theta(n^{3/2})\), so both
the signing entropy and the energy contribute at speed \(n^2\).

Put \(f_n=F(n)/n^{3/2}\).  The minimum term and the total number of
terms give

\[
e^{-\beta n^2 f_n}
\le Z_n(\beta)
\le 2^{E_n}e^{-\beta n^2 f_n}.
\]

Consequently,

\[
\boxed{
-\frac{\Phi_n(\beta)}{\beta}
\le f_n
\le
\frac{(E_n/n^2)\log2-\Phi_n(\beta)}{\beta}.
}                                                            \tag{0.2}
\]

If \(\Phi_n(\beta)\) converged for an unbounded set of fixed
\(\beta\)'s, then

\[
\limsup f_n-\liminf f_n
\le\frac{\log2}{2\beta}
\]

for every such \(\beta\), proving convergence of \(f_n\).  Thus the
pressure criterion is correct.

---

## 1. Fixed-threshold heredity and exact Shearer inequality

Define

\[
\mathcal A_n(T)=\{A:M(A)\le T\},
\qquad
C_n(T)=|\mathcal A_n(T)|.
\]

Principal restriction cannot increase \(M\).  Hence
\(\{\mathcal A_n(T)\}_n\) is hereditary for every **fixed absolute**
threshold \(T\).

Let \(N>n\), and choose a uniform member of \(\mathcal A_N(T)\).
Applying Shearer's entropy inequality to all induced \(n\)-vertex
subgraphs gives

\[
\boxed{
\frac{\log C_N(T)}{E_N}
\le
\frac{\log C_n(T)}{E_n}.
}                                                            \tag{1.1}
\]

This is an exact theorem.  At the natural threshold
\(T=cN^{3/2}\), however, it becomes

\[
\boxed{
\frac1{E_N}\log C_N(cN^{3/2})
\le
\frac1{E_n}
\log C_n\!\left(
c\left(\frac Nn\right)^{3/2}n^{3/2}
\right).
}                                                            \tag{1.2}
\]

Thus restriction changes the normalized threshold by
\((N/n)^{3/2}\).  This is the same changing-temperature obstruction
that appears in block free-energy arguments.

### Gibbs/relative-entropy form

Normalize by the uniform signing measure:

\[
Y_n(\lambda)=2^{-E_n}\mathscr Z_n(\lambda).
\]

For a distribution \(\mu\) on order-\(N\) signings,

\[
\log Y_N(\lambda)
=
\sup_\mu
\left\{
-D(\mu\|u_N)-\lambda\,\mathbb E_\mu M(A)
\right\}.
\]

Let \(S\) be a uniform \(n\)-vertex subset and
\(\mu_S\) the induced marginal.  Entropy Shearer and
\(M(A_S)\le M(A)\) give, with \(p=E_n/E_N\),

\[
-D(\mu\|u_N)-\lambda\mathbb E M(A)
\le
\frac1p\,
\mathbb E_S
\left[
-D(\mu_S\|u_n)-\lambda p\,\mathbb E M(A_S)
\right].
\]

Taking the variational supremum yields

\[
\boxed{
\frac1{E_N}\log Y_N(\lambda)
\le
\frac1{E_n}\log Y_n(\lambda p).
}                                                            \tag{1.3}
\]

At \(\lambda=\beta\sqrt N\), the smaller-order inverse temperature is

\[
\boxed{
\beta'
=
\frac{\lambda p}{\sqrt n}
=
\beta\,\frac{E_n}{E_N}\sqrt{\frac Nn}
=
\beta\left(\frac nN\right)^{3/2}(1+o(1)).
}                                                            \tag{1.4}
\]

No choice of a fixed proportional restriction closes at the same
\(\beta\).

---

## 2. Layer-cake form: why integrating the threshold does not repair it

For nonnegative energies,

\[
\boxed{
\mathscr Z_n(\lambda)
=
\lambda\int_0^\infty e^{-\lambda T}C_n(T)\,dT.
}                                                            \tag{2.1}
\]

At \(\lambda=\beta\sqrt n\), put \(T=cn^{3/2}\).  Up to subexponential
factors, (2.1) is the variational expression

\[
\frac1{E_n}\log\mathscr Z_n(\beta\sqrt n)
\sim
\sup_c
\left\{
\frac1{E_n}\log C_n(cn^{3/2})-2\beta c
\right\}.                                                    \tag{2.2}
\]

Equation (1.2) controls the entropy profile only after replacing
\(c\) by \(c(N/n)^{3/2}\).  Laplace optimization does not remove that
rescaling; it converts it exactly into (1.3)--(1.4).

There is also no useful equicontinuity of the threshold entropy at
speed \(n^2\) coming from edge Lipschitzness.  Flipping \(r\) edges
changes \(M\) by at most \(2r\), but changing the normalized threshold
by a fixed \(\delta\) permits

\[
r=\Theta(\delta n^{3/2})
\]

edge flips.  A Hamming neighborhood of that radius has logarithmic
volume only \(O(n^{3/2}\log n)=o(n^2)\).  It therefore need not change
the speed-\(n^2\) entropy at all.  The entropy profile may have genuine
jumps or frozen phases at this scale.

---

## 3. Exact gauge-fixed cavity recursion

Let

\[
B^s=\operatorname{diag}(s)B\operatorname{diag}(s).
\]

For a new signed row \(b\), the exact row identity is

\[
M(B,b)
=
\max_x\left(
|H_B(x)|+|b\cdot x|
\right).
                                                               \tag{3.1}
\]

Define

\[
\widehat M(B)=M(B,\mathbf1)
=
\max_x\left(
|H_B(x)|+\left|\sum_i x_i\right|
\right).
\]

Switching the old vertices by \(b\) gives

\[
M(B,b)=\widehat M(B^b).
\]

As \(B\) ranges over all signings, so does \(B^b\).  Therefore the
entire new row can be gauge-fixed:

\[
\boxed{
\mathscr Z_{n+1}(\lambda)
=
2^n\sum_Be^{-\lambda\widehat M(B)}.
}                                                            \tag{3.2}
\]

Equivalently, with

\[
\Delta_B(b)=M(B,b)-M(B),
\]

and the order-\(n\) Gibbs measure

\[
\nu_{n,\lambda}(B)
\propto e^{-\lambda M(B)},
\]

switching invariance gives

\[
\boxed{
\mathscr Z_{n+1}(\lambda)
=
2^n\mathscr Z_n(\lambda)
\,
\mathbb E_{\nu_{n,\lambda}}
e^{-\lambda\Delta_B(\mathbf1)}.
}                                                            \tag{3.3}
\]

This is the exact signing-space analogue of a cavity identity.

Unfortunately,

\[
0\le\widehat M(B)-M(B)\le n,                                \tag{3.4}
\]

and the upper endpoint is real: switch any ground state of \(B\) to
\(\mathbf1\), after which the magnetization term adds exactly \(n\).
Thus the affine correction is not uniformly \(O(\sqrt n)\).

From (3.2)--(3.4),

\[
\boxed{
2^ne^{-\lambda n}\mathscr Z_n(\lambda)
\le
\mathscr Z_{n+1}(\lambda)
\le
2^n\mathscr Z_n(\lambda).
}                                                            \tag{3.5}
\]

After dividing by the uniform signing counts,

\[
\boxed{
e^{-\lambda n}Y_n(\lambda)
\le Y_{n+1}(\lambda)\le Y_n(\lambda).
}                                                            \tag{3.6}
\]

Thus \(Y_n(\lambda)\) is genuinely monotone for every fixed raw
\(\lambda\).  Since it is also decreasing in \(\lambda\),

\[
Y_{n+1}(\beta\sqrt{n+1})
\le Y_n(\beta\sqrt n).
\]

This still does not control
\(-n^{-2}\log Y_n(\beta\sqrt n)\): a positive increasing sequence of
exponents can have an oscillating quotient by \(n^2\), as Section 5
makes explicit.

At \(\lambda=\beta\sqrt{n+1}\), the lower uncertainty in
\(\log\mathscr Z\) is \(O_\beta(n^{3/2})\), or
\(O_\beta(n^{-1/2})\) after division by \(n^2\).  This proves adjacent
continuity of the pressure but is not summable and does not prove a
limit.

More importantly, (3.3) has not closed the state.  The factor

\[
\mathbb E_{\nu_{n,\lambda}}
e^{-\lambda\Delta_B(\mathbf1)}
\]

depends on the complete joint energy-versus-magnetization profile of
\(B\), not on \(M(B)\) or on finitely many spectral moments.  The
already verified order-\(6\) examples with the same \(M\) and different
affine norms are the finite manifestation of this failure.

Gauge symmetry therefore removes the explicit sum over new rows, but
it does not make the cavity increment universal; it is an exact
reindexing rather than an Aizenman--Sims--Starr closure theorem.

---

## 4. Uniform regularity in \(\beta\)

For every \(\beta_0>0\), the family
\(\{\Phi_n(\beta):\beta\ge\beta_0\}\) is uniformly Lipschitz on compact
\(\beta\)-intervals.

Indeed, under the Gibbs measure at inverse temperature \(\beta\),
the variational formula and the single minimizing signing give

\[
\beta\sqrt n\,\mathbb E_\beta M(A)
\le
E_n\log2+\beta\sqrt n\,F(n).
\]

Using \(F(n)=O(n^{3/2})\),

\[
\frac{\mathbb E_\beta M(A)}{n^{3/2}}
\le
O(1)+\frac{\log2}{2\beta}.
                                                               \tag{4.1}
\]

Since

\[
\Phi_n'(\beta)
=
-\frac{\mathbb E_\beta M(A)}{n^{3/2}},
\]

the claim follows.

This regularity lets one replace
\(\beta(1-O(1/n))\) by \(\beta\) in a single adjacent comparison.
It does not repair the harmonic accumulation over all orders, and the
countermodel below shows that no argument using only this regularity
and Shearer can do so.

---

## 5. Exact diagonal countermodel to the pressure axioms

Choose constants

\[
0<c_0-\varepsilon<c_0+\varepsilon<\frac12
\]

and, for sufficiently large \(n\), let

\[
m_n=
\text{the nearest admissible integer to }
n^{3/2}
\left[
c_0+\varepsilon\sin(\log\log n)
\right].
                                                               \tag{5.1}
\]

The harmless rounding can enforce any desired parity.  Since

\[
\frac{d}{dn}
\left(
n^{3/2}[c_0+\varepsilon\sin(\log\log n)]
\right)
=
\sqrt n
\left[
\frac32(c_0+\varepsilon\sin(\log\log n))
+
\frac{\varepsilon\cos(\log\log n)}{\log n}
\right],
\]

\(m_n\) is eventually strictly increasing and

\[
0\le m_{n+1}-m_n=O(\sqrt n)\le n.
\]

Now give every one of the \(2^{E_n}\) abstract signings the same
energy \(m_n\).  This model has:

* vertex-permutation and switching invariance;
* the exact signing entropy \(E_n\log2\);
* hereditary fixed-threshold families;
* monotonicity under restriction;
* the one-vertex increment bound \(m_{n+1}\le m_n+n\);
* the exact threshold Shearer inequality (1.1);
* the Gibbs Shearer inequality (1.3);
* convex, decreasing, locally uniformly Lipschitz pressures.

Its pressure is

\[
\boxed{
\Phi_n^{\rm diag}(\beta)
=
\frac{E_n}{n^2}\log2
-
\beta\frac{m_n}{n^{3/2}},
}                                                            \tag{5.2}
\]

which does not converge for any \(\beta>0\).

For example, (1.3) reduces exactly to \(m_N\ge m_n\):

\[
-\frac{\lambda m_N}{E_N}
\le
-\frac{\lambda p\,m_n}{E_n}
=
-\frac{\lambda m_n}{E_N}.
\]

Thus fixed-threshold heredity, entropy contraction, changing-temperature
Shearer, symmetry, convexity, and local padding are jointly
insufficient to prove pressure convergence.  The countermodel does not
satisfy the special quadratic row profile (3.1); that profile is
precisely where any successful proof must obtain new information.

---

## 6. Verdict and the missing theorem

The signing-pressure strategy remains a valid sufficient criterion for
the original limit, but its proposed generic interpolation mechanisms
do not establish the criterion.

The exact results are:

1. fixed-threshold counts have a monotone entropy density;
2. Gibbs pressures obey the exact rescaled Shearer inequality (1.3);
3. switching symmetry gives the exact gauge-fixed cavity recursion
   (3.2)--(3.3);
4. the pressure is uniformly regular away from \(\beta=0\);
5. all of those properties still permit diagonal oscillation.

An actual convergence proof now needs a specifically quadratic cavity
estimate, beginning with an order-\(n\) asymptotic such as

\[
\log
\mathbb E_{\nu_{n,\beta\sqrt n}}
e^{-\beta\sqrt n\,\Delta_B(\mathbf1)}
=
n\,\Gamma(\beta)+o(n)
                                                               \tag{6.1}
\]

with a limit \(\Gamma(\beta)\) independent of subsequence, together
with compatible control of the \(O(1/n)\) temperature shift in the
core pressure.  Equivalently, one needs a projective variational
principle with a unique value for the full energy--magnetization
profile.  Proving merely
\(\widehat M(B)-M(B)=o(n)\) on average is not enough unless its error is
quantified strongly enough to make the accumulated pressure error
summable.

No pressure limit was proved.  The main progress is that the exact
Shearer and gauge identities isolate the only remaining non-generic
input: convergence of the tilted cavity spectrum.  The diagonal
countermodel proves that this input cannot be replaced by abstract
heredity or convexity.
