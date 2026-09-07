# Wave 56: a common near-ground active-face law for all edge blocks

## Outcome

Let `P` be any parent edge universe incident to `v` vertices and let `r` be a
block budget.  Exact signing minimality, Bernstein rounding, and finite
minimax produce **one** parent-state law `mu` satisfying

```math
\mathbb E_\mu\Delta
+4h_r\big((\mathbb E_\mu s_e)_{e\in P}\big)\le\eta_r,
```

where `h_r` is the sum of the `r` largest positive coordinates and

```math
\eta_r=\sqrt{32r(v+1)\log2}+\frac83(v+1)\log2.
```

After conditioning, if `r>>n` and `r=O(n^(3/2))`, the same law may be
supported entirely on an `o(T_n)` parent-deficit face while every edge subset
of size at most `r` has expected signed correlation at most `o(r)`.

Taking `P` to be all parent edges and `r=Theta(n^(3/2))` gives a law
independent of the selector.  Every child ground has an excess-edge set of
size `Theta(n^(3/2))`; at least `1/3-o(1)` of the common law escapes each such
set.  Fubini therefore gives one `o(T_n)`-deficit parent state escaping a
constant fraction of all canonically chosen child fibres.

This is genuine common-law/multiplicity structure.  It still does not prove
bare favorability or a low project row.  A tiny support, including an
opposite-orientation response, can satisfy every marginal balance.  These are
the exact remaining obstructions.

## 1. Fractional perturbations have an approximate common response

Let `A` be an exact order-`n` minimizer, `q=Q(A)`, and on oriented parent
states write

```math
E_A(\omega)=\tau x^{\mathsf T}Ax,qquad
\Delta(\omega)=q-E_A(\omega),qquad
s_e(\omega)=\tau a_ex_ix_j.
```

Fix an edge universe `P`, incident to `v` vertices, and the fractional
uniform-matroid polytope

```math
\mathcal P_r=\{p\in[0,1]^P:\ \sum_ep_e\le r\}.
```

For `p in P_r`, replace `a_e` fractionally by `(1-2p_e)a_e`.  Independently
round it to `-a_e` with probability `p_e` and to `a_e` otherwise.  For any
fixed oriented edge pattern, the rounding error in matrix energy is a sum of
mean-zero terms with absolute value at most four and total variance at most
`16r`.  There are at most `2^v` distinct oriented patterns on the incident
vertices.  Bernstein and a union bound, with `L=(v+1)log2`, give a rounding
whose error is at most

```math
\boxed{
\eta_r=\sqrt{32r(v+1)\log2}+\frac83(v+1)\log2
}
\tag{A56.1}
```

simultaneously on every parent state.  The rounded integral signing has cap
at least `q` by exact minimality, so the fractional signing has cap at least
`q-eta_r`.  Equivalently, for every `p in P_r`,

```math
\boxed{
\min_\omega\left\{
\Delta(\omega)+4\sum_{e\in P}p_es_e(\omega)
\right\}\le\eta_r.
}
\tag{A56.2}
```

No assumption about independence of parent states is used.  The rounding is
only a separation device which imports global integral minimality into the
fractional perturbation game.

## 2. Minimax yields one law balancing every block

For a vector `z in R^P`, define

```math
h_r(z)=\max_{p\in\mathcal P_r}\sum_ep_ez_e.
```

It is the sum of the `floor(r)` largest positive coordinates plus the
appropriate fraction of the next one.  The perturbation polytope is compact
and convex, and the parent state space is finite.  Applying von Neumann/Sion
minimax to (A56.2) proves the **Verified common active-face theorem**:

```math
\boxed{
\exists\mu:\qquad
\mathbb E_\mu\Delta
+4h_r\big((\mathbb E_\mu s_e)_{e\in P}\big)
\le\eta_r.
}
\tag{A56.3}
```

In particular,

```math
\mathbb E_\mu\Delta\le\eta_r,
\qquad
h_r\big((\mathbb E_\mu s_e)_e\big)\le\eta_r/4.
\tag{A56.4}
```

This is strictly stronger than choosing a different migrated state after
each block: the same law controls every set (and every fractional set) of
budget at most `r`.

## 3. Conditioning onto a subproject deficit face

Let `epsilon in (0,1)` and condition `mu` on

```math
G=\{\Delta\le\eta_r/\epsilon\}.
```

Markov gives `mu(G)>=1-epsilon`.  If `m_G=(E[s_e|G])_e`, then, for every
`p in P_r`, the discarded part contributes at worst `-epsilon r`, so

```math
\boxed{
h_r(m_G)\le
\frac{\eta_r/4+\epsilon r}{1-\epsilon}.
}
\tag{A56.5}
```

Suppose `v=O(n)`, `n<<r<=C n^(3/2)`, and `0<c<1/4`.  Then

```math
\eta_r=O(\sqrt{rn}+n)=o(r),
\qquad \eta_r=o(T_n),
\qquad T_n=n^{3/2-c}.
```

Choose `epsilon_n->0` with `eta_r/T_n=o(epsilon_n)`.  The conditioned law,
still denoted `mu`, is supported on

```math
\boxed{
\Delta=o(T_n),qquad
h_r\big((\mathbb E_\mu s_e)_e\big)=o(r).
}
\tag{A56.6}
```

If necessary choose `epsilon_n` also larger than `sqrt(eta_r/r)`; both
requirements are compatible because the two ratios vanish.

## 4. A selector-independent law escapes every child excess block

Take `P` to be the full parent edge set and choose one deterministic budget
`r=C_1n^(3/2)` larger than `q_n` for all large `n`.  On a compact density
window, every child ground has

```math
q_m\le Q(A[S])\le q_n,
\qquad q_m=\Omega(n^{3/2}),quad q_n=O(n^{3/2}).
```

The excess-edge construction (10.1324), now using all child vertices, gives
a child-positive set `E_S` satisfying

```math
\frac12Q(A[S])\le|E_S|\le Q(A[S]).
```

Thus `|E_S|=Theta(r)`, uniformly in `S`, and its indicator belongs to
`P_r`.  Equation (A56.6) gives

```math
\mathbb E_\mu S_{E_S}=o(|E_S|).
\tag{A56.7}
```

Since `S_{E_S}` lies in `[-|E_S|,|E_S|]`,

```math
\boxed{
\mu\{S_{E_S}\le|E_S|/2\}\ge\frac13-o(1)
}
\tag{A56.8}
```

uniformly for every selector and any deterministic choice of one child
ground.  Indeed, if the complementary probability is `theta`, the expectation
is at least `(3theta/2-1)|E_S|`.

Averaging (A56.8) over uniform selectors proves the **Verified common escape
corollary**:

```math
\boxed{
\text{some parent state with }\Delta=o(T_n)
\text{ escapes at least }(1/3-o(1))\text{ of the chosen child fibres.}
}
\tag{A56.9}
```

Here escape means that more than one quarter of the child-positive excess
edges have negative parent sign, because
`S_E=|E|-2|{e in E:s_e=-1}|`.

## 5. Exact scope and obstruction

The theorem balances marginal signed edge responses, not support size,
project row, or child restriction.  A response game with one state satisfying
`Delta=0` and `s_e=-1` on every edge has the left side of (A56.3) equal zero
and satisfies every conclusion with support one.  In the signing problem the
orientation bit makes this obstruction structurally relevant: an opposite-
orientation near-ground state can make an entire child-positive block
negative without sharing the child fibre.

Furthermore, escaping `E_S` controls only one selected internal edge sum.
It neither bounds the complete local deficit

```math
\delta_S(d)=Q(A[S])-\tau x_S^{\mathsf T}A[S]x_S
```

nor the project row `x^T A^2x`.  Consequently the common escape incidence
has `O(1)` selector--state information cost, but is not known to lie in the
bare favorable event and has only the generic row cap `O(n^(5/2))`.

The exact successor target is to add either:

1. a comparison showing that a saved fraction of the common escape incidence
   is bare-favorable and has captured row `O(n^(9/4-c))`; or
2. a two-sided/sector-sensitive perturbation certificate which rules out the
   tiny opposite-orientation response and controls the full local deficit.

Without one of these additions, common-law compatibility alone does not prove
the box witness or selected fixed-cut tail.
