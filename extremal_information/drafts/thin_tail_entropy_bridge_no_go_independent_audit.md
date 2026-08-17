# Independent audit: thin-tail entropy bridge no-go

**Verdict: PASS.**

The conference normalization, bridge-output cap event, use of Theorem 36.26,
negative-moment conditioning, and abstract disorder-cube calculation are all
correct.  The theorem really does show that the archived fixed-small-tilt
conference wall persists after restricting to an overwhelming class on
which every complete-sign output has one common bounded-cap thin-tail
certificate.

The source audited and frozen for this verdict is

```text
extremal_information/drafts/thin_tail_entropy_bridge_no_go.md
sha256 4aa3e9cf08c2cf02d200eac4fbe781da0cae66ab1cf93bbcff6f54aac2ff0155
```

No mathematical repair is required.  For maximal self-containment, Section
4 could explicitly write that its abstract pressures are
`L(H)=log E_x cosh(beta H(x)/sqrt(N))`; this is the unique normalization used
in (CT.29) and is already clear from the surrounding text.  The omission is
not a logical or constant error.

## 1. Half-quadratic and conference-pressure mapping: PASS

For the block signing

```math
S_{\epsilon,B}=
\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
```

the repository's unordered-edge convention gives exactly

```math
H_{S_{\epsilon,B}}(x,y)
=H_{A_r}(x)+\epsilon H_{A_r}(y)+x^TBy.
```

There is no factor two in the cross term: it appears twice in `z^TSz` and
is then divided by two.  The normalized partition function is

```math
\overline Z_{2r}(S,t)=2^{-2r}\sum_z\cosh(tH_S(z)),
\qquad t=\frac\beta{\sqrt{2r}},
```

which is exactly the convention in the audited two-temperature theorem.

The source's quantities match Theorem 6.1 of
`artifacts/two_temperature_bridge_audit.md` term for term:

```math
T_r=2\log\overline Z_r(A_r,\beta/\sqrt r),
\qquad
h_\beta=2\psi(\beta/\sqrt2)+\beta^2/4,
```

and

```math
h_\beta-2\psi(\beta)
=\beta^2/4-2\psi(\beta)+2\psi(\beta/\sqrt2)
=\gamma(\beta)>0.
```

Thus the archived conclusion is equivalently

```math
\mathbb E_{\epsilon,B}e^{-\lambda L_{\epsilon,B}}
=e^{-\lambda h_\beta r+o(r)},
```

for every fixed `0<lambda<lambda_0(beta)` in
`0<beta<sqrt(2)/6`, together with `T_r/r->2psi(beta)`.  The source invokes
this only along the same Paley conference sequence; it does not claim the
pressure asymptotic from the matrix identity (CT.3) for every abstract
conference signing.

## 2. Hoeffding event and parent cap: PASS

For fixed Boolean `x,y`, the cross energy `x^TBy` is a sum of `r^2`
independent Rademacher variables.  The standard two-sided Hoeffding bound
at `u=2r^(3/2)` gives

```math
\Pr\{|x^TBy|>u\}
\le2\exp\{-u^2/(2r^2)\}
=2e^{-2r}.
```

There are `2^(2r)` ordered spin pairs, so

```math
\Pr(\mathcal E_r^c)
\le2^{2r}\,2e^{-2r}
=2e^{-(2-2\log2)r}.
```

Hence `c_0=2-2log2` and the prefactor in (CT.13) are correct.

From `A_r^2=(r-1)I`,

```math
|H_{A_r}(x)|
\le\frac12\|A_r\|_{op}\|x\|_2^2
=\frac r2\sqrt{r-1}.
```

The two internal blocks therefore cost at most `r sqrt(r-1)`, while the
cross term costs at most `2r^(3/2)` on `E_r`.  Thus

```math
Q(S_{\epsilon,B})
\le r\sqrt{r-1}+2r^{3/2}
\le3r^{3/2}
<2(2r)^{3/2}.
```

All matrices are complete hollow exact signings of order `N=2r`, for both
orientations.

## 3. Theorem 36.26 constants and quantifiers: PASS

The last display places every retained parent in the class

```math
Q(S)\le C N^{3/2}
```

with the common safe choice `C=2`.  The explicit theorem constant is

```math
d_C=\frac1{200000\max\{C,1\}},
```

so `d_2=1/400000` is correct.  Theorem 36.26 is uniform over all complete
signings in the class and its stated absolute-tail conclusion has already
absorbed the union of the positive and negative endpoint layers.  It
therefore yields (CT.19) simultaneously for every retained `B` and each
`epsilon`; no further probabilistic intersection or orientation loss is
needed.

## 4. Conditioning the negative moment: PASS

The normalized average of `cosh` is at least one, hence
`L_(epsilon,B)>=0` and `e^(-lambda L)<=1`.  Therefore the entire discarded
negative moment is at most

```math
\mathbb E[1_{\mathcal E_r^c}e^{-\lambda L}]
\le2e^{-c_0r}.
```

The full moment has rate `e^(-lambda h_beta r+o(r))`.  If

```math
0<\lambda<\min\{\lambda_0(\beta),c_0/(2h_\beta)\},
```

then `c_0-lambda h_beta>c_0/2`, and the discarded/full ratio is

```math
O\left(e^{-(c_0-\lambda h_\beta)r+o(r)}\right)=o(1).
```

Since `Pr(E_r)=1-o(1)`, conditioning changes neither exponential rate:

```math
\mathbb E_{U_r^E}e^{-\lambda L}
=e^{-\lambda h_\beta r+o(r)}.
```

Taking `-(1/lambda)log` gives `R^E_(lambda,r)/r->h_beta`; subtracting the
child target gives the positive limit `gamma(beta)`.  The strict inequality
`lambda<lambda_*` is important and is correctly stated.  The result rules
out only a nonempty fixed-small-tilt interval, not every fixed or growing
tilt.

## 5. Abstract disorder-cube example: PASS

Read the Section 4 pressures at the SK parameter

```math
L(H)=\log\mathbb E_x
 \cosh\left(\frac{\beta H(x)}{\sqrt N}\right).
```

The good output has two antipodal top states and the bad outputs have

```math
K_N=2\left\lfloor e^{(\log2-\kappa)N}/2\right\rfloor
```

top states, with all other energies zero.  Because `c>d`, precisely these
top states lie in the declared `dN^(3/2)` layer, so every output satisfies
(CT.28).  Since `beta c>log2>kappa`, the top contribution dominates the
unit baseline in both partition functions, and direct calculation gives

```math
L_{good}=(\beta c-\log2)N+o(N),
\qquad
L_{bad}=(\beta c-\kappa)N+o(N),
```

hence (CT.29).

If the unique good output has mass `p_N=2^(-J_N)` and
`Delta_N=L_bad-L_good`, then exactly

```math
\mathcal R_{\lambda,N}
=L_{bad}-{1\over\lambda}
 \log\{1+p_N(e^{\lambda\Delta_N}-1)\}.
```

Using `log(1+v)<=v` and
`(e^u-1)/u<=e^u` gives

```math
0\le L_{bad}-\mathcal R_{\lambda,N}
\le p_N\Delta_Ne^{\lambda\Delta_N}.
```

Here `J_N=Theta(N^2)` and `Delta_N=Theta(N)`.  Thus every positive
`lambda_N=o(N)` makes the right side exponentially smaller than `N`, while
bringing the isolated output into the order-`N` soft minimum requires
`lambda_N Delta_N` to compete with `J_N log2`, hence
`lambda_N=Omega(N)`.  The calculation is exact and the source correctly
labels it an abstract bridge-output model rather than a signing or physical
bridge construction.

## 6. Novelty and non-overclaim: PASS

The archived two-temperature audit already proved the unconditioned
conference wall, its Renyi identity, and the need for a linearly large
low-pressure basin.  CT.1--CT.2 add one precise fact: after discarding every
output outside an explicit `1-exp(-Theta(r))` bounded-cap event, all
remaining actual-signing outputs satisfy the later thin-tail theorem and
the same negative-moment rate and linear defect persist.  That conditional
intersection is not in the archived theorem.

The separate scalar-pressure consequence gives a broader abstract scalar
no-go; this draft appropriately does not claim its pressure sandwich as a
second novelty.  Nor is there a conflict with BCL.0: BCL retains a named
root and constructs a matched bridge, whereas the present entropy tilt sees
only the scalar distribution of bridge-output pressures.

The result does **not** exclude:

- a larger fixed tilt above the proved threshold;
- a growing tilt;
- another bridge distribution or a signing-specific diffuse basin;
- a joint spin/disorder large-deviation state;
- energy-resolved overlap or rooted response information;
- ordinary exact bridge minimization.

It also does not falsify compensation specialized to actual minimizers of
the contracted-temperature pressure: the Paley conference children are a
rigorous structured benchmark, but are not known to be those minimizers.
The final source now states this both in its opening scope and in its route
judgment.

It proves no recurrence, pressure limit, convergence theorem, or change in
the rigorous interval for `M_n/n^(3/2)`.  Its stopping judgment is therefore
properly scoped: a scalar fixed spin-tail deficit cannot by itself supply
the missing bridge-disorder basin rate.
