# Wave 18: switching groups at the critical scale

This note uses the doubled quadratic normalization from the ledger.  All
temporary computations are in `check_symplectic_switching_r18.py` and
`audit_symplectic_switching_r18.py`.

## 1. Exact Walsh normalization

Let `G <= {+-1}^V` be an elementary abelian switching group.  Write
`chi_i in G^*` for the coordinate character `chi_i(g)=g_i`, and let
`V_alpha={i:chi_i=alpha}` be the nonempty character cells.  If `D` is the
within-cell part of a signing `A` and `C=A-D`, then, for an oriented state
`(sigma,z)`, define

```math
W_eta(z)=
2\sigma\sum_{\substack{\alpha<\beta\\\alpha\beta=\eta}}
z_\alpha^{\mathsf T}A[V_\alpha,V_\beta]z_\beta
\qquad(\eta\in G^*,\ \eta\ne1).
```

Here `alpha<beta` is any fixed ordering; in additive `F_2^r` notation the
condition is `alpha+beta=eta`.  Directly grouping the cross edges gives

```math
Y_g:=\sigma(z\odot g)^{\mathsf T}C(z\odot g)
=\sum_{\eta\ne1}W_\eta(z)\eta(g).
```

Character orthogonality, with normalized counting measure on `G`, gives the
exact normalization

```math
\frac1{|G|}\sum_gY_g=0,
\qquad
\frac1{|G|}\sum_gY_g^2
=\sum_{\eta\ne1}W_\eta(z)^2.
```

Thus the variance in (10.611) is precisely the squared `l_2` norm of the
nontrivial difference-frequency coefficients.  Reynolds averaging alone
says only that this Walsh polynomial has mean zero.  It supplies no bound
making its minimum close to minus its maximum.

There is a sharp distinction between matrix and statewise antipodality.  If
some `g in G` obeyed `C^g=-C` as matrices, write `epsilon_alpha` for its
constant sign on `V_alpha`.  Since every cross-cell entry is nonzero,

```math
epsilon_alpha epsilon_beta=-1
```

for every two distinct nonempty cells.  Three cells are impossible: the
relations for `(alpha,beta)` and `(alpha,gamma)` give
`epsilon_beta=epsilon_gamma`, contradicting the relation for
`(beta,gamma)`.  Hence exact **matrix** antipodality permits at most two
nonempty character cells.  This does not obstruct a state-specific pair
`z,z odot g` with opposite cross energies; the construction below has many
cells and exactly such a pair.

## 2. A quotient-collision sufficient theorem

Let `Omega={+-1}^V/{+-1}` be the projective state group, use the image of
`G` in `Omega` (so a global-sign kernel is removed), and fix one orientation
`sigma`.  For the cross matrix `C`, put

```math
P=\{p\in\Omega:\sigma p^{\mathsf T}Cp=Q(C)\},
```

and

```math
N_t=\{w\in\Omega:\sigma w^{\mathsf T}Cw\le-Q(C)+t\}.
```

For the quotient `Omega/G`, let `p_j=|P cap j|`,
`n_j=|N_t cap j|`, and `L=|Omega/G|`.  Define the relative subgroup
energies

```math
E_G(P)=\sum_jp_j^2,
\qquad
E_G(N_t)=\sum_jn_j^2.
```

Then the number of opposite near-ground pairs whose difference lies in the
actual switching group is

```math
J_G(P,N_t)
=\sum_jp_jn_j
=\sum_{g\in G}|\{(p,w)\in P\times N_t:p\odot w=g\}|.
```

Centering the two quotient count vectors and applying Cauchy--Schwarz gives

```math
J_G(P,N_t)\ge
\frac{|P||N_t|}{L}
-\sqrt{\left(E_G(P)-\frac{|P|^2}{L}\right)
       \left(E_G(N_t)-\frac{|N_t|^2}{L}\right)}.
```

Consequently the strict quotient-mixing condition

```math
\left(E_G(P)-\frac{|P|^2}{L}\right)
\left(E_G(N_t)-\frac{|N_t|^2}{L}\right)
<\frac{|P|^2|N_t|^2}{L^2}
```

implies `J_G(P,N_t)>0`.  Choose a colliding pair `w=p odot g`.  Its
opposite-orbit deficit is at most `t`, because

```math
Q(C)+\sigma w^{\mathsf T}Cw\le t.
```

The Reynolds response argument (10.610) therefore proves

```math
\Gamma_{C,t}(D)\le t/2.
```

This is a fully rigorous sufficient theorem in terms of **relative**
additive energy (equivalently, mixing after quotienting by the proposed
switching group).  It is not a theorem from cardinality alone.  Combined
with maximum cell size `O(sqrt(n))`, `X=Omega(n^(3/2))`, and
`t=o(n^(3/2))`, it proves exactly the spatial target (10.600), but still no
allocation or temporal conclusion.

## 3. What universal low-field near states do and do not give

Here is an explicit scale audit.  Let `M` be symmetric, zero diagonal, with
`|m_ij|<=1`, and suppose `P(M)<=K n^(3/2)`.  Switch a positive ground to
`1`, and write `r_i=sum_j m_ij`.  One-spin optimality gives `r_i>=0` and
`sum_i r_i=P(M)`.  Therefore

```math
L_0=\{i:r_i\le2K\sqrt n\}
```

has at least `n/2` vertices.  If `S subset L_0`, `|S|<=k`, and `x^S`
is obtained by flipping `S`, put

```math
R_S=\sum_{i\in S}r_i,
\qquad
h_S=1_S^{\mathsf T}M[S]1_S.
```

The exact cut expansion and the entrywise bound give

```math
0\le P(M)-(x^S)^{\mathsf T}Mx^S
=4(R_S-h_S)
\le8Kk\sqrt n+4k(k-1).
```

Thus `k=o(n^(3/4))` produces a universally guaranteed cloud of one-sided
`o(n^(3/2))`-near states (and the count is subexponential for such `k`).
This does **not** provide an opposite-oriented cloud or quotient mixing.

There is also a precise large-cell obstruction if one tries to use only
this guaranteed cloud as generators.  Differences of two radius-`k`
members are supported on at most `2k` coordinates.  If `r` such differences
generate a subgroup, every coordinate outside their support union has the
same trivial character.  Hence one character cell has size at least

```math
n-2rk.
```

In particular, for `r=O(log n)` and `k=o(n/log n)`, this is `n-o(n)`, not
`O(sqrt(n))`.  This is only a no-go for the universally guaranteed
low-field flip construction at small rank.  It says nothing against dense
near grounds, a rank `Omega(n/k)` subgroup, or a separately structured
affine near-ground family.

Known universal ground/near-ground information in the ledger does not imply
the quotient-mixing hypothesis.  Exact ground families can be small and
Sidon at finite minimizers; the universal low-field cloud above is
one-sided and `2^{o(n)}`.  When `rank(G)=o(n)`, the quotient has
`L=2^{n-1-rank(G)}` cosets, so cardinality alone is far below even the
random-collision scale `|P||N_t| about L`.  Small aligned sets can still
collide, as the next example shows; what is missing is structural alignment,
not merely more states.

## 4. Exact critical-scale symplectic construction

The three spatial requirements in (10.600) are algebraically compatible at
the exact scale.  This construction is a competitive signing, **not known
to be a global minimizer**.

Let `k=2^m`, `n=k^2`, and index coordinates by
`(u,v) in F_2^m times F_2^m`.  Define

```math
K_{(u,v),(x,y)}=(-1)^{v\mathbin\cdot x+u\mathbin\cdot y},
\qquad A=K-I.
```

The alternating form makes `diag(K)=1`, so `A` is a valid zero-diagonal
signing.  Let

```math
G=\{g_h:g_h(u,v)=(-1)^{h(u)},\ h:F_2^m\to F_2\}.
```

Its `k` character cells are `V_u={u} times F_2^m`, all of order
`k=sqrt(n)`.  Put `r_u(v)=(-1)^{u dot v}` and

```math
P=\bigoplus_u r_ur_u^{\mathsf T}.
```

The Reynolds pieces are exactly

```math
D=P-I,
\qquad C=K-P.
```

Indeed, every within-fibre block is `r_u r_u^T-I`, while averaging over
all functions `h` kills every edge between distinct fibres.  Let
`U=im(P)`.  Direct character orthogonality gives

```math
K^2=nI,
\qquad P^2=kP,
\qquad KP=PK=kP.
```

Thus `K=kI` on `U`, while `U^perp` is invariant, and

```math
C=0\text{ on }U,
\qquad C=K\text{ on }U^\perp,
\qquad \|C\|_{op}=k.
```

For `b ne 0`, define Boolean vectors

```math
z_{a,b}(u,v)=(-1)^{u\mathbin\cdot v+a\mathbin\cdot u+b\mathbin\cdot v}.
```

Their fibre correlation with `r_u` is zero, and summing first over the
`x` variable gives

```math
Pz_{a,b}=0,
\qquad
Kz_{a,b}=k(-1)^{a\mathbin\cdot b}z_{a,b}.
```

Hence `Q(C)=nk=n^(3/2)`.  Choose `a dot b=0` and `c dot b=1`.  Then

```math
w=z_{a+c,b}=z_{a,b}\odot g_c,
\qquad Cz=kz,
\qquad Cw=-kw.
```

This is an exact opposite pair in one actual switching orbit, so its
antipodal gap is zero.  Moreover `Dz=Dw=-z,-w`, respectively.  Under the
two opposite orientations the `D` responses are `-n` and `+n`; hence

```math
\mathcal R_{C,0}(D)\ge n,
\qquad
\Gamma_{C,0}(D)=0.
```

Every fibre block has norm `k(k-1)`.  Therefore, using the ledger's
universal `q_k<=kappa k^(3/2)` bound,

```math
\begin{aligned}
Q(D)=S&=\sum_uQ(A[V_u])=k^2(k-1)=k^3-k^2,\\
X&=S-kq_k
\ge k^3-k^2-\kappa k^{5/2}
=(1-o(1))n^{3/2}.
\end{aligned}
```

Thus `s=sqrt(n)`, `X=Omega(n^(3/2))`, `t=0`, and `Gamma=0` all hold
simultaneously.

For scope, `A` itself has

```math
Q(A)=n(k+1).
```

The upper bound is spectral, and the negative vector `w` attains it.  This
is the correct competitive `n^(3/2)+O(n)` scale, but exact global minimality
already fails at `k=2`, since `Q(A)=12` whereas `q_4=8`, and is not claimed
asymptotically.  (The conference bound `n sqrt(n-1)` is not available at
these orders `n=4^m`; one must not use it here without a conference matrix
of the required order.)  The construction therefore defeats any algebraic
no-go based only on Reynolds projection, signs, or scaling; it does not
solve the minimizer-existence step.

## 5. The same orbit gives no paired raw-shore resource

The construction also sharply separates (10.600) from (10.604).  Every
positive exact `C` ground lies in the `+k` eigenspace and every negative one
in the `-k` eigenspace, both inside `U^perp`; the two eigenspaces are
orthogonal.  Since `D=-I` on `U^perp`, every opposite exact cross-ground
pair `(p,w)` satisfies

```math
p^{\mathsf T}Aw=0.
```

Its full-`A` deficits are `2n` on the positive `C` side and `0` on the
negative side.  Consequently

```math
\mathcal B_{C,0}(A)=0.
```

For the explicit same-orbit pair above, agreement and disagreement each
have `n/2` vertices.  Indeed, `Az=(k-1)z` and
`Aw=-(k+1)w`.  In the `z` gauge, let `a_i` and `b_i` be the signed row
fields from the shore containing `i` and from the opposite shore.  The two
eigenvector equations become

```math
a_i+b_i=k-1,
\qquad
a_i-b_i=-(k+1),
```

so `a_i=-1` and `b_i=k` at every vertex.  Summing within either shore and
taking the `l_1` norm of its opposite-shore field gives

```math
h_X=-\frac n2,
\qquad
L_X=\frac{nk}{2},
\qquad
Q(A)=n(k+1).
```

Thus all four raw decrement-tolled buckets are zero:

```math
[2L_X+h_X-Q(A)]_+
=[2L_X-h_X-Q(A)]_+=0.
```

So even exact antipodality, leading captured excess, macroscopic shores,
and zero active defect need not select paired raw service.  This statement
is about the raw static buckets only; it supplies no allocation or temporal
claim.

## 6. Scoped status

- **Verified algebraically:** the Walsh normalization, the at-most-two-cell
  matrix-antipodality obstruction, the quotient-collision theorem, and the
  low-field sparse-generator bound.
- **Verified construction:** the symplectic family satisfies every spatial
  condition in (10.600) and has zero paired raw functional/buckets.  Matrix
  identities were checked for `k=2,4,8,16`; exhaustive Boolean enumeration
  independently checked `k=2,4`.
- **Open minimizer input:** prove that every exact global minimizer admits a
  critical-cell subgroup satisfying the quotient collision and leading
  excess conditions, or derive an alternative minimizer-specific source of
  a symplectic-like fibre certificate.  No currently known universal
  ground/near-ground count or additive-energy statement supplies this.
