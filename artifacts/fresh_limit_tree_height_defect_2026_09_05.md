# A uniform approximate-fixed-point obstruction from tree height

Date: 2026-09-05. This is a Gaussian function-space theorem for the
hierarchical tree isometry. It does not make a claim about convergence of
the original signing minima.

## 1. Exact height-projection identity

Give each rooted tree its ordinary height, counting edges from the external
root. Thus the single-edge tree has height one. Let `P_L` be the orthogonal
projection in Gaussian first chaos onto coordinates `Z_T` of height at most
`L`; set `P_0=0`. A fixed height contains infinitely many tree types, but
these projections and their Gaussian sigma-fields are well-defined closed
subspaces. Every finite tree has finite height, so `P_L V→V` in `L²`.

Let `V` be any unit-variance first-chaos Gaussian and let `g` be an even
Gaussian `L²` function with squared norm one. Write

\[
 g(z)=\sum_{r\text{ even}} b_r\operatorname{He}_r(z)/\sqrt{r!},
 \qquad K_g(q)=\sum_r b_r^2q^r,\quad 0\le q\le1.
\]

Put `m_L=||P_L V||₂²`, so `m0=0` and `m_L↑1`. The child-polynomial
basis gives the exact identity

\[
 \boxed{\quad
 \|P_L\mathcal U[g(V)]\|_2^2=K_g(m_{L-1}).
 \quad}                                                       \tag{1}
\]

Indeed, parent trees of height at most `L` correspond exactly to all even
Hermite monomials in child coordinates of height at most `L-1`. Parseval
therefore identifies the left side with
`||E[g(V)|𝔽_(L-1)]||₂²`, where `𝔽_(L-1)` is the Gaussian sigma-field
of those child coordinates. The decomposition
`V=P_(L-1)V+(I-P_(L-1))V` has independent Gaussian parts of variances
`m_(L-1)` and `1-m_(L-1)`. The usual Hermite conditional-expectation formula
then gives (1), including `m=0` and `m=1` by continuity.

## 2. Uniform residual barrier

For every such `g`,

\[
 \boxed{\quad
 \inf_{\substack{V\text{ first chaos}\\EV^2=1}}
 \|V-\mathcal U[g(V)]\|_2
 \ge
 \sup_{0<c<1}\bigl(\sqrt c-\sqrt{K_g(c)}\bigr)_+.
 \quad}                                                       \tag{2}
\]

To prove this, write the residual as `ε`. Projection and the reverse
triangle inequality imply, for every `L≥1`,

\[
 \sqrt{m_L}\le\sqrt{K_g(m_{L-1})}+\epsilon.            \tag{3}
\]

If `ε<sqrt(c)-sqrt(K_g(c))` for some `c∈(0,1)`, monotonicity of `K_g`
and induction from `m0=0` show `m_L<c` for every finite `L`. This contradicts
`m_L↑1`. Thus (2) follows. The statement is uniform over all first-chaos
directions and does not assume their coefficients are supported on a fixed
finite set of trees.

In particular, if `Σ_r r b_r²>1`, allowing infinity, then the right side
of (2) is positive. Indeed

\[
 {1-K_g(c)\over1-c}
   =\sum_r b_r^2(1+c+\cdots+c^{r-1})
   \longrightarrow\sum_r r b_r^2>1
\]

as `c↑1`; hence `K_g(c)<c` for some `c<1`. Consequently a supercritical
scalar response cannot have even approximate unit first-chaos fixed points.
This strengthens the decorated Galton--Watson nonexistence criterion for
exact fixed points. No assertion of sharpness of (2) is made.

## 3. A strict gap below the naive one-mask rearrangement ceiling

Recall `J(H)=E|𝒰H|(1-H)` for even `0≤H≤1`. Let

\[
 C_*:=\max_{0<p<1}2\sqrt p\,
           \phi\!\left(\Phi^{-1}((1+p)/2)\right)
        \approx0.44955496156.
\]

Then, in fact,

\[
                    \sup_{0\le H\le1}J(H)<C_*.       \tag{4}
\]

Here the supremum is over the full jointly-even Gaussian mask space, not
only scalar fixed-point masks. This paragraph proves strictness, not an
explicit decimal improvement in the ceiling.

Suppose to the contrary that `J(H_k)→C*`. Set `p_k=EH_k`, `v_k=EH_k²`,
and `V_k=𝒰H_k/sqrt(v_k)`. The rearrangement bounds

\[
 J(H_k)\le2\sqrt{v_k}\phi(\alpha_{p_k})
             \le2\sqrt{p_k}\phi(\alpha_{p_k})\le C_*
\]

force `p_k→p*`, `v_k→p*`, where `p*` is the unique maximizer, and force
asymptotic equality in the first rearrangement. The uniqueness follows
from the strictly increasing stationary equation
`α[2Φ(α)-1]-φ(α)=0`.

For completeness, rearrangement equality is stable uniformly here. With
`K_k=1{|V_k|≤α_(p_k)}`, its nonnegative loss is
`d_k=E[(|V_k|-α_(p_k))(H_k-K_k)]→0`. For any `η>0`,

\[
 E|H_k-K_k|\le d_k/\eta+4\phi(0)\eta.
\]

This follows by separating the `η`-band around the threshold and using
the bounded folded-normal density. Thus `||H_k-K_k||₂→0`. Continuity of
the Gaussian quantile and `v_k→p*` give

\[
 \left\|{H_k\over\sqrt{v_k}}-
 {1_{\{|V_k|\le\alpha_*\}}\over\sqrt{p_*}}\right\|_2\to0.
\]

Apply `𝒰` and set `g_*(z)=1{|z|≤α*}/sqrt(p*)`. This gives
`||V_k-𝒰[g_*(V_k)]||₂→0`, contradicting (2): the nonconstant interval
indicator has infinite Gaussian Dirichlet energy, so its barrier is
strictly positive. Therefore (4) holds.

The conclusion concerns this particular orientation-averaged cross-energy
mechanism. It neither upper-bounds the original normalized minima nor
resolves their convergence.
