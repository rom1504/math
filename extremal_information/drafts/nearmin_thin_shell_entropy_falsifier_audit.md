# Independent audit: thin-shell entropy falsifier

Date: 2026-08-17.

Audited file: `drafts/nearmin_thin_shell_entropy_falsifier.md`.

## Verdict

**REPAIR (minor scope and wording repairs; all four mathematical cores
pass).**  TS.1's constants and projectivized count, TS.2's Walsh
normalization and asymptotic count, both inclusions in TS.3, and TS.4's
deterministic inequality are correct.  The repairs below prevent three
readings that are stronger than what is proved: unrestricted PC.3 scope,
sharpness of the edit constants, and a matching upper count for the Walsh
shell.

## TS.1

Let `x` be oriented so that `H_A(x)=Q`.  One-spin optimality gives
`ell_i>=0`, and

```math
\sum_i\ell_i=x^{\mathsf T}Ax=2Q.
```

Consequently at least `floor(n/2)` indices obey `ell_i<=4Q/n`.
For `|S|=k`, the exact identity is

```math
H_A(x^S)=Q-2\sum_{i\in S}\ell_i
             +4\sum_{\{i,j\}\subseteq S}s_{ij}.
```

Since the last sum is at least `-binom(k,2)`, this gives exactly

```math
Q-|H_A(x^S)|
 <=Q-H_A(x^S)
 <=8kQ/n+2k(k-1).
```

For fixed `c` and `k=floor(c sqrt(n))`, the stated deficit is
`(8cC+2c^2+o(1))n`.  The count is also correct:

```math
\log_2 {\lfloor n/2\rfloor\choose k}
=(c/2+o(1))\sqrt n\log_2n.
```

Two members `x^S,x^T` with `|S|=|T|=k<n/2` cannot be global negatives,
so projectivization loses nothing before the orientation pigeonhole.  The
factor `1/2` is therefore valid.  Zero energies can be assigned to either
orientation; alternatively the known `Q=Omega(n^(3/2))` lower bound makes
all these energies nonzero for large `n` at the displayed `O(n)` deficit.

**TS.1 passes without a mathematical patch.**

## TS.2 and the two imported inputs

With Sylvester entries `W_d(u,v)=(-1)^(u dot v)`, one has

```math
W_d^2=NI,\qquad \operatorname{tr}W_d=0
```

for `d>=1`.  Thus for `A_d=W_d-diag(W_d)`,

```math
2H_{A_d}(z)=z^{\mathsf T}A_dz
=z^{\mathsf T}W_dz-\operatorname{tr}W_d
=z^{\mathsf T}W_dz.
```

The operator bound is `|z^T W_d z|<=N^(3/2)`.  Equality with positive
sign holds exactly for Boolean vectors in the `+sqrt(N)` eigenspace, which
are precisely signs of self-dual bent functions under the paper's Walsh
normalization.

The Carlet--Danielsen--Parker--Sole input is correctly used.  In Theorem
4.9 take an arbitrary bent `f`, `f_2=tilde f`,
`g_1(y)=y_1y_2` (self-dual), and
`g_2(y)=y_1y_2+y_1+y_2` (anti-self-dual).  Their indirect sum is self-dual
and, in lexicographic `y` blocks, has sign vector

```math
(F,\widetilde F,\widetilde F,-F).
```

The first block makes the map injective.  Complementary `f` give globally
opposite lifted vectors and no other projective collision occurs, hence
there are at least `b_(d-2)/2` positive projective maximizers.

The source of arXiv:2508.14605v3 was checked directly.  Its main theorem is

```math
\log_2 b_t\ge t2^{t/2}(1+O(1/t))
```

for even `t>=4`.  Substituting `t=d-2` gives

```math
\log_2(b_{d-2}/2)
\ge (1/2+o(1))\sqrt N\log_2N.
```

The PC.3 conjugacy is also exact, but its scope should be stated.  Its seed
is `H_16=D_bW_4D_b`; hence

```math
H_j=H_{16}^{\otimes j}
=(D_b^{\otimes j})W_{4j}(D_b^{\otimes j}).
```

Boolean eigenvectors are carried bijectively by `D_b^(tensor j)`, and the
trace-zero diagonal deletion does not change Boolean quadratic energy.
This proves TS.2 for the **unflipped hollow PC.3 children**
`H_j-diag(H_j)`.  It says nothing about PC.3 sparse-flipped children or
shore-completed parents.

Recommended precise patches:

1. Introduce the asymptotic counting statement for even `d>=6` (or state
   the exact Walsh identity separately for even `d>=2` and the imported
   lower bound for `d>=6`).
2. Replace “the diagonal-conjugate PC.3 tensor tower” by “the unflipped
   hollow PC.3 child tower.”
3. Clarify the order-16 check: there are `20` positive and `20` negative
   **vectors**, hence `10` positive and `10` negative projective classes.

## TS.3

The uniform edge-edit estimate is exactly `|H_A-H_B|<=2r`.  Minimality of
`A` gives

```math
Q(A)=M_n<=Q(B)<=M_n+2r.
```

These yield the stated asymmetric constants:

```math
S_B(Delta)\subseteq S_A(Delta+2r),\qquad
S_A(Delta)\subseteq S_B(Delta+4r).
```

They also justify the Walsh-halo consequence (TS.15).  What is not supplied
is an extremal example proving that `2r` and `4r` are sharp.  Rename the
lemma “black-box shell stability under edge edits,” or add a sharpness
example.  No change to the inclusions is needed.

## TS.4

For `S subseteq J`, nonnegativity of the local fields gives

```math
Q-H_A(x^S)
 =2\sum_{i\in S}\ell_i-4\sum_{\{i,j\}\subseteq S}s_{ij}
 <=2L(J)+4P(J).
```

This is stronger than the displayed absolute-deficit statement and proves
it.  For genuine near-minimizers, the known uniform lower bound
`Q(A_n)=Omega(n^(3/2))` and `L(J_n)+P(J_n)=O(n)` imply that **every** one of
these energies is positive for large `n`.  The `2^|J|` spins therefore have
one orientation, and projectivization leaves at least `2^(|J|-1)` classes
(including the edge case `J=[n]`).  Add this one sentence to the proof; if
the `Omega(n^(3/2))` premise is not intended, only `2^(|J|-2)`
same-orientation projective classes follow in complete generality.

For a quantitatively fixed proposed upper bound
`exp(O(sqrt(n)(log n)^a))`, write the falsifier as

```math
|J_n|/(\sqrt n(\log n)^a)\longrightarrow\infty.
```

The current phrase `gg sqrt(n) polylog(n)` is informal because “polylog”
does not specify `a`.

## Scope and overclaim audit

The final distinction among exact-shell cardinality, thin-shell
cardinality, and contextual response entropy is correct and important.
No response packing follows from TS.1 or TS.2.  In particular:

- the Haugland/Carlet argument proves only a lower bound
  `2^(Omega(sqrt(N) log N))` for the exact Walsh shell;
- it does not prove that this is its true order, and an `exp(cN)` set of
  self-dual Boolean eigenvectors remains compatible with the cited bounds;
- the PC.3 statement is about one unflipped structured Level-4 family, not
  genuine near-minimizers;
- a large set of witnesses may have a concise common generator and hence
  small contextual information.

Accordingly replace “Walsh/PC.3 ... saturates the lower scale” in the
director judgment by “Walsh/PC.3 realizes the universal lower scale already
at exact cap.”  “Saturates” can otherwise be misread as the unproved
matching exact-shell upper bound.
