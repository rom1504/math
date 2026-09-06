# Wave 11 memo: exact boundary-state recursion for the augmented cut code

Date: 2026-07-29  
Disposition: **verified structural theorem / exact adaptation; no convergence proof**  
Tracked files changed: none

## 1. Sources searched and selected result

I searched both suggested directions. The selected result is a recent exact
covering-radius/conditional-decoding theorem:

- Karthik Sheshadri, *Trellis State Complexity as an Exact Tropical
  Factorization Rank*, arXiv:2607.23471v1 (26 July 2026),
  <https://arxiv.org/abs/2607.23471> and
  <https://arxiv.org/html/2607.23471>.

This is a very recent v1, and its acknowledgment reports LLM assistance. I
therefore did not treat it as a black box: the complete argument needed here
is reconstructed in Section 2 and independently checked in Section 5.

Contextual primary sources:

- Patrick Solé and Thomas Zaslavsky, *A Coding Approach to Signed Graphs*,
  SIAM J. Discrete Math. 7 (1994), 544--553,
  <https://doi.org/10.1137/S0895480189174374>, author PDF
  <https://people.math.binghamton.edu/zaslav/Tpapers/cas.sidma1994.pdf>.
  This supplies the classical switching-class/cut-code/coset-leader
  dictionary. It has covering-radius bounds, but I found no block recursion
  that removes our boundary state.
- G. David Forney, Jr., *Dimension/Length Profiles and Trellis Complexity of
  Linear Block Codes*, IEEE Trans. Inform. Theory 40 (1994), 1741--1752,
  <https://doi.org/10.1109/18.340452>. This is classical background for the
  trellis state-complexity formula.
- Francesco Guerra and Fabio Toninelli, *The Thermodynamic Limit in Mean Field
  Spin Glass Models*, Comm. Math. Phys. 230 (2002), 71--79,
  <https://arxiv.org/abs/cond-mat/0204280>. Its fixed-temperature Gaussian
  interpolation is compared with our setting in Section 7.

## 2. Selected theorem, reconstructed

Let $C\subseteq\mathbb F_2^e$ be a binary linear code and split its
coordinates as $L\sqcup R$. Define the supported subcodes

$$
C_L=\{c\in C:c_R=0\},\qquad C_R=\{c\in C:c_L=0\},
$$

and the conditional decoding table

$$
W(a_L,a_R)=d((a_L,a_R),C).
$$

Put

$$
s=\dim C-\dim C_L-\dim C_R.
$$

The theorem is that the min-plus factorization rank of $W$ is exactly
$2^s$. In fact its tropical and Kapranov ranks are also $2^s$.

Here min-plus rank is the least $r$ for which

$$
W(a_L,a_R)=\min_{1\leq k\leq r}
\{u_k(a_L)+v_k(a_R)\}.
$$

### 2.1 Upper bound: the door identity

Let $P_R(C)$ be the right projection of $C$ and define the state quotient

$$
\mathcal T=P_R(C)/C_R.
$$

Its dimension is

$$
\dim\mathcal T
=\dim P_R(C)-\dim C_R
=\dim C-\dim C_L-\dim C_R=s.
$$

For $r\in P_R(C)$ let

$$
F(r)=\{c_L:(c_L,r)\in C\}.
$$

If $k\in C_R$, addition of $(0,k)$ gives a bijection
$F(r)\to F(r+k)$, so $F(r)$ depends only on the class
$\tau=r+C_R$. Set

$$
D(a_L,\tau)=d(a_L,F(r)),\qquad
d(a_R,\tau)=\min_{r\in\tau}|a_R+r|.
$$

Grouping the minimization over $c\in C$ by its right projection gives the
exact **door identity**

$$
\boxed{
W(a_L,a_R)=\min_{\tau\in\mathcal T}
\{D(a_L,\tau)+d(a_R,\tau)\}.
}
$$

There are $|\mathcal T|=2^s$ terms, proving rank at most $2^s$.

### 2.2 Lower bound: transversal block and crossing

Choose representatives $r_1,\dots,r_N$ of $\mathcal T$, where $N=2^s$,
and lifts $c^i=(c_L^i,r_i)\in C$. Consider the $N\times N$ submatrix
whose row $i$ is indexed by $c_L^i$ and column $j$ by $r_j$.

Its diagonal is zero. If an off-diagonal entry were zero, then
$(c_L^i,r_j)\in C$. Subtracting $c^i$ would give

$$
(0,r_i+r_j)\in C_R,
$$

contrary to $r_i$ and $r_j$ representing distinct quotient classes.
Thus

$$
W(i,i)=0,\qquad W(i,j)\geq1\quad(i\ne j).
$$

Every term $u_k(i)+v_k(j)$ in a valid min-plus factorization is at least
$W(i,j)$. If one term were tight at diagonal cells $i$ and $j$, then

$$
\begin{aligned}
0
&=[u_k(i)+v_k(i)]+[u_k(j)+v_k(j)]\\
&=[u_k(i)+v_k(j)]+[u_k(j)+v_k(i)]\\
&\geq W(i,j)+W(j,i)\geq2,
\end{aligned}
$$

a contradiction. Each diagonal cell therefore needs a distinct term, so
the rank is at least $N=2^s$.

The same submatrix is tropically nonsingular: the identity permutation has
weight zero, while every other permutation moves at least two indices and
has weight at least two. Hence tropical rank is at least $2^s$; the standard
chain tropical rank $\leq$ Kapranov rank $\leq$ min-plus rank makes all
three equal.

## 3. Exact specialization to the augmented cut code

Let $E_n=\binom n2$ and

$$
\mathcal C_n=
\{c(t,z):c_{ij}=t+z_i+z_j,\ t,z_i\in\mathbb F_2\}.
$$

For $n\geq3$, the parameter map has the one-dimensional kernel obtained by
flipping every $z_i$, so

$$
\dim\mathcal C_n=n.
$$

Split the vertices into shores $S\sqcup R$, of sizes $m,\ell\geq3$, and
split edge coordinates into

$$
I=E(S)\sqcup E(R),\qquad X=E(S,R).
$$

### 3.1 Subcode supported on the internal edges

If $c_X=0$, then

$$
t+z_i+z_j=0\qquad(i\in S,j\in R).
$$

Thus $z$ is constant, say $a$ on $S$ and $b$ on $R$, and $t=a+b$.
Every internal coordinate then equals $t$. The internal-supported subcode is
therefore $\{0,\mathbf1_I\}$ and

$$
\dim (\mathcal C_n)_I=1.
$$

### 3.2 Subcode supported on the cross edges

If $c_I=0$, use any triangle within either shore. Summing its three
equations over $\mathbb F_2$ gives $t=0$, after which all $z_i$ on that
shore are equal. Thus $z$ is constant separately on the two shores and the
cross restriction is either zero or all ones. Consequently

$$
\dim (\mathcal C_n)_X=1.
$$

It follows that the exact boundary-state dimension and count are

$$
\boxed{s=n-2,\qquad |\mathcal T|=2^{n-2}.}
$$

Equivalently, a state is a pair of projective shore-spin patterns

$$
[x]\in\{\pm1\}^S/\{\pm1\},\qquad
[y]\in\{\pm1\}^R/\{\pm1\},
$$

carrying $(m-1)+(\ell-1)=n-2$ bits.

**Consequence.** Any universal exact separable min-plus recursion for the
full conditional value table needs $2^{n-2}$ boundary terms. There is no
constant-state or polynomial-state exact recursion of this form. This is a
representation theorem, not a computational-hardness theorem, and it does
not exclude compression restricted to globally minimizing signings or an
approximation at the $o(n^{3/2})$ scale.

## 4. The door identity is exactly the ledger's block identity

Represent a target word by signs

$$
A_{ij}=(-1)^{a_{ij}}\quad(i,j\in S),\qquad
D_{ij}=(-1)^{a_{ij}}\quad(i,j\in R),\qquad
B_{ij}=(-1)^{a_{ij}}\quad(i\in S,j\in R).
$$

Write

$$
H_A(x)=\sum_{i<j\in S}A_{ij}x_ix_j,
\qquad
H_D(y)=\sum_{i<j\in R}D_{ij}y_iy_j.
$$

For the state $\tau=([x],[y])$, the two door costs are

$$
D_I(a_I,\tau)
=\frac{|I|-|H_A(x)+H_D(y)|}{2},
$$

and

$$
d_X(a_X,\tau)
=\frac{m\ell-|x^\top By|}{2}.
$$

The absolute values occur because each relevant fiber contains a word and
its all-ones complement. Substitution in the door identity gives

$$
W(a_I,a_X)
=\frac{E_n}{2}
-\frac12\max_{x,y}
\left\{|H_A(x)+H_D(y)|+|x^\top By|\right\}.
$$

Since $M(G)=E_n-2W(a)$, this is precisely

$$
\boxed{
M(G)=\max_{x,y}
\left\{|H_A(x)+H_D(y)|+|x^\top By|\right\},
}
$$

the exact identity in ledger Section 1.6. Thus the coding theorem rigorously
identifies the missing boundary state, but its scalar projection returns us
to the already known gluing obstruction.

In particular, choosing optimal internal signings and then a best
rectangular cross signing only yields

$$
M_{m+\ell}\leq M_m+M_\ell+R_{m,\ell},
\qquad
R_{m,\ell}=\min_B\max_{x,y}|x^\top By|.
$$

For proportional shores, $R_{m,\ell}$ is of order
$\sqrt{m\ell(m+\ell)}=\Theta(n^{3/2})$. This error is at the leading scale,
so ordinary scalar Fekete theory still gives no convergence.

## 5. Independent finite audit

The checker is
`/home/math/quadra/tmp/verify_cutcode_boundary_r11.py`. It constructs
$\mathcal C_n$, both supported subcodes, the quotient states and fibers,
then compares the door formula directly with distance to the full code. It
also checks every entry of the transversal block used in the lower bound.

Command:

```bash
.venv/bin/python tmp/verify_cutcode_boundary_r11.py
```

Output:

```text
3+3: dim C=6, dim C_I=1, dim C_X=1, state dimension=4, states=16; door (exhaustive)/block PASS
3+4: dim C=7, dim C_I=1, dim C_X=1, state dimension=5, states=32; door (5000 sampled entries)/block PASS
```

The $3+3$ conditional table was checked exhaustively. For $3+4$, 5,000
deterministic pseudorandom table entries and the entire $32\times32$
transversal block were checked.

## 6. Exact canonical-pressure reformulation: a reciprocal boundary sum

Let $a$ range over all edge words. Since

$$
M(a)=E_n-2W(a),
$$

the canonical signing partition function becomes

$$
\mathfrak Z_n(\beta)
=e^{-\beta\sqrt n E_n}
\sum_a e^{\lambda_n W(a)},
\qquad \lambda_n=2\beta\sqrt n.
$$

For the door costs

$$
c_\tau(a)=D_I(a_I,\tau)+d_X(a_X,\tau),
\qquad W(a)=\min_\tau c_\tau(a),
$$

define the ordinary positive boundary sum

$$
S_a=\sum_{\tau\in\mathcal T}e^{-\lambda_n c_\tau(a)}.
$$

Because one term attains the minimum and every term is at least the minimum,

$$
e^{-\lambda_nW(a)}\leq S_a
\leq |\mathcal T|e^{-\lambda_nW(a)}.
$$

After inversion,

$$
\boxed{
\frac{e^{\lambda_nW(a)}}{|\mathcal T|}
\leq \frac1{S_a}\leq e^{\lambda_nW(a)}.
}
$$

Therefore replacing $e^{\lambda_nW(a)}$ by $S_a^{-1}$ changes
$\log\mathfrak Z_n$ by at most

$$
\log|\mathcal T|=(n-2)\log2=o(n^2).
$$

This gives the pressure-equivalent expression

$$
\boxed{
\Phi_n(\beta)
=-\beta\frac{E_n}{n^{3/2}}
+\frac1{n^2}\log
\sum_{a_I,a_X}
\frac1{
\sum_{\tau}
e^{-\lambda_nD_I(a_I,\tau)}
e^{-\lambda_nd_X(a_X,\tau)}
}
+o(1).
}
$$

The inner boundary sum is a sum of separable products, but the required
**reciprocal** destroys the usual sum-product factorization. Equivalently,
the signing pressure asks for a negative first moment of a boundary
partition function. The exponential state count itself costs only $O(n)$
in the logarithm and is harmless on the $n^2$ pressure scale; the real
obstruction is the inverse inner product, which retains the full geometry
of the boundary profiles. This is the coding analogue of the cap-intersection
state in ledger Section 10.19.

This formula suggests a concrete new target:

$$
\frac1{n^2}\log\sum_{a_I,a_X}
\left[
\sum_\tau e^{-\lambda_n(D_I(a_I,\tau)+d_X(a_X,\tau))}
\right]^{-1}.
$$

One would need either a scale-compatible recursion for this negative moment
or a coarse boundary factorization whose error is $o(n^2)$ in its log.

## 7. Explicit comparison with the two known walls

### 7.1 Scalar Fekete wall

At zero temperature, the door theorem has the correct exact signs and no
temperature issue. But eliminating the $2^{n-2}$ boundary coordinates gives
exactly

$$
|H_A+H_D|+|x^\top By|,
$$

and a scalar bound replaces the last term by $R_{m,\ell}=\Theta(n^{3/2})$.
That is the ledger Section 3.3 obstruction, not a smaller remainder.

### 7.2 Wrong-temperature Shearer wall

Ledger Section 10.48 gives, for restriction from order $N$ to order $m$,

$$
\beta_{\rm Shearer}
=\beta q\sqrt{N/m}
\sim\beta(m/N)^{3/2}.
$$

The boundary formulation has a different but still wrong characteristic.
Its raw distance fugacity is

$$
\lambda_N=2\beta\sqrt N.
$$

If an $m$-vertex child distance is read with this same coefficient, it
corresponds to the child's canonical parameter

$$
2\beta_{\rm child}\sqrt m=2\beta\sqrt N,
\qquad
\boxed{\beta_{\rm child}=\beta\sqrt{N/m}.}
$$

Thus the Shearer restriction sends a smaller child to lower $\beta$, whereas
the raw door cost sends it to higher $\beta$. Neither preserves fixed
$\beta$. Rescaling the child costs to force fixed temperature breaks their
exact additivity and reintroduces a leading-order block correction.

### 7.3 Why Guerra--Toninelli does not directly close this formula

For the SK model, Guerra--Toninelli interpolate a full Gaussian Hamiltonian
against two independently normalized blocks at the **same** $\beta$.
Gaussian integration by parts yields

$$
\frac d{dt}\frac1N\mathbb E\log Z_N(t)
=-\frac{\beta^2}{4}
\left\langle
q_{12}^2-\frac{N_1}{N}(q_{12}^{(1)})^2
-\frac{N_2}{N}(q_{12}^{(2)})^2
\right\rangle\geq0,
$$

where the sign is supplied by convexity of $u\mapsto u^2$. The endpoints
factor at fixed temperature.

Here the disorder is the object being summed over, the energy contains an
optimization over spins, and the exact boundary reduction produces
$S_a^{-1}$ rather than a positive partition sum. There is currently no
Gaussian integration-by-parts identity or overlap convexity that gives a
sign for this negative boundary moment, and the endpoint coefficient has
the child scaling above. Consequently the classical interpolation does not
repair either wall merely by changing notation.

## 8. Final disposition

What is proved:

1. The exact conditional covering-radius recursion for the augmented cut
   code has state dimension $n-2$ and exactly $2^{n-2}$ indispensable
   min-plus states.
2. Its door identity is exactly the known block-gluing identity, with all
   signs and factors of two accounted for.
3. At canonical scale, it yields a pressure-equivalent reciprocal boundary
   partition with only $O(n)$ state-entropy loss.

What is falsified as a route:

- An exact universal constant-state or polynomial-state min-plus boundary
  recursion for the full conditional table.
- A scalar collapse of the exact recursion as a way around the
  $\Theta(n^{3/2})$ Fekete remainder.

What remains genuinely open:

- Approximate/coarse state compression at the much larger allowed error
  scale relevant to $o(n^2)$ log-pressure or $o(n^{3/2})$ energy.
- A negative-moment interpolation or inequality for the reciprocal boundary
  sum.
- Compression restricted to globally minimizing or near-minimizing
  signings. The exact rank lower bound uses only an off-diagonal gap of one,
  so it supplies no robust obstruction at these asymptotic error scales.
