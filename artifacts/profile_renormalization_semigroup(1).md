# Profile-valued finite-fibre renormalization

## Status

This note identifies the smallest natural quadratic profile on which every
finite-fibre lift
\[
\mathcal T_{R,D}(C)=C\otimes R+I_n\otimes D \tag{1}
\]
acts exactly, associatively, and continuously.  It is a strict quotient of
the Boolean spectral profile in
`boolean_spectral_profile_compactification.md`: only overlaps and first
spectral moments are needed for finite-fibre quadratic composition.

The construction is exact, but it does **not** prove convergence of the
minima.  The obstruction is also exact:

- a lift reads higher tuple levels of the input profile;
- iteration shifts to tuple levels \(s^d\);
- compactness does not control this shift;
- and minimizing profiles are not known to be invariant under any lift.

There are concrete period-two and scalar-branching counterexamples below.
Thus an absorbing-minimizer theorem would be a genuinely new theorem, not
a consequence of compact-semigroup formalism.

## 1. Minimal quadratic action profile

Let \(C\) be an order-\(n\) symmetric zero-diagonal signing.  For an
\(r\)-tuple of Boolean vectors, written as
\[
X=(x^1,\ldots,x^r)\in\{\pm1\}^{n\times r},
\]
define its overlap and normalized action matrices by
\[
O_C(X)=\frac1nX^\top X,\qquad
H_C(X)=\frac1{n^{3/2}}X^\top CX. \tag{2}
\]
Set
\[
\mathcal P_r(C)
=\{(O_C(X),H_C(X)):X\in\{\pm1\}^{n\times r}\}. \tag{3}
\]
The full quadratic action profile is
\[
\mathfrak P(C)=(\mathcal P_r(C))_{r\ge1}. \tag{4}
\]

If \(\|C\|_{\mathrm{op}}\le K\sqrt n\), then every entry of \(H\) lies
in \([-K,K]\), while \(O\) is a correlation matrix.  Hence each
\(\mathcal P_r(C)\) lies in a fixed compact finite-dimensional space.
Putting the Hausdorff metric on each level and, for example,
\[
d(\mathfrak P,\mathfrak Q)
=\sum_{r\ge1}2^{-r}
\min\{1,d_H(\mathcal P_r,\mathcal Q_r)\}, \tag{5}
\]
gives a compact product hyperspace after closure.

The normalized Boolean norm is the continuous functional
\[
\boxed{
\Phi(\mathfrak P(C))
=\sup_{(O,H)\in\mathcal P_1(C)}|H_{11}|
=\frac{Q(C)}{n^{3/2}}.
}
\tag{6}
\]

### Relation to the Boolean spectral profile

For the matrix-valued spectral measure \(\mathbf M\) in
`boolean_spectral_profile_compactification.md`,
\[
O=\mathbf M(\mathbb R),\qquad
H=\int t\,d\mathbf M(t). \tag{7}
\]
Thus (4) is a continuous quotient of that richer profile.  The empirical
joint spin law and higher spectral moments are unnecessary for (1):
every term in a finite-fibre quadratic lift is a linear contraction of
an overlap or a bilinear action.

This quotient is minimal in the following operational sense:

- one output spin for an \(s\)-fibre lift requires the jointly attainable
  \(s\)-tuple data \(\mathcal P_s\);
- \(r\) output spins require \(\mathcal P_{sr}\);
- exact finite-\(n\) composition needs both \(O\) and \(H\), because the
  internal filling \(D\) contracts \(O\);
- arbitrary iteration forces all tuple levels \(s^d r\).

A single scalar, or any fixed maximum tuple order, is therefore not closed
under repeated finite-fibre composition.

## 2. Exact renormalization map

Let the fibre size be \(s\), let
\[
R\in\{\pm1\}^{s\times s}
\]
be symmetric, including on its diagonal, and let \(D\) be a symmetric
zero-diagonal signing of order \(s\).  Introduce the scale coordinate
\[
\eta=n^{-1/2}. \tag{8}
\]

For \((O,H)\in\mathcal P_{sr}(C)\), index its rows and columns by
\((p,a)\), where \(1\le p\le r\) and \(1\le a\le s\).  Define
\[
\widetilde O_{pq}
=\frac1s\sum_{a=1}^s O_{(p,a),(q,a)}, \tag{9}
\]
and
\[
\widetilde H_{pq}
=\frac1{s^{3/2}}
\left[
\sum_{a,b=1}^sR_{ab}H_{(p,a),(q,b)}
+\eta\sum_{a,b=1}^sD_{ab}O_{(p,a),(q,b)}
\right]. \tag{10}
\]
Let \(\Lambda_{R,D,\eta}^{(r)}\) denote this linear contraction.

### Exact profile theorem

If \(B=\mathcal T_{R,D}(C)\), of order \(ns\), then
\[
\boxed{
\mathcal P_r(B)
=
\Lambda_{R,D,\eta}^{(r)}(\mathcal P_{sr}(C))
\quad\text{for every }r\ge1.
}
\tag{11}
\]

To prove (11), write the \(p\)-th Boolean vector on \([n]\times[s]\)
as \(y^p_{i,a}=x^{p,a}_i\).  This is a bijection between \(r\) output
vectors and \(sr\) input vectors.  Direct expansion gives
\[
\frac1{ns}(y^p)^\top y^q
=\frac1s\sum_aO_{(p,a),(q,a)}
\]
and
\[
\begin{aligned}
\frac1{(ns)^{3/2}}(y^p)^\top By^q
={}&\frac1{s^{3/2}}\sum_{a,b}R_{ab}H_{(p,a),(q,b)}\\
&+\frac{\eta}{s^{3/2}}\sum_{a,b}D_{ab}O_{(p,a),(q,b)}.
\end{aligned}
\]
These are exactly (9)--(10), and every input tuple arises.

Define the state map
\[
\mathscr R_{R,D}(\eta,\mathfrak P)
=
\left(
\frac{\eta}{\sqrt s},
\bigl(
\Lambda_{R,D,\eta}^{(r)}(\mathcal P_{sr})
\bigr)_{r\ge1}
\right). \tag{12}
\]
Equation (11) says
\[
\mathscr R_{R,D}(n^{-1/2},\mathfrak P(C))
=((ns)^{-1/2},\mathfrak P(\mathcal T_{R,D}(C))). \tag{13}
\]

Each output level in (12) is the continuous linear image of one compact
input level.  Linear images are continuous in Hausdorff distance, so
\(\mathscr R_{R,D}\) is continuous in the product topology (5).

## 3. Exact associativity

For fibre pairs, define
\[
(R_1,D_1)\star(R_2,D_2)
=
\left(
R_1\otimes R_2,\,
D_1\otimes R_2+I_{s_1}\otimes D_2
\right). \tag{14}
\]
Then
\[
\mathcal T_{R_2,D_2}(\mathcal T_{R_1,D_1}(C))
=
\mathcal T_{(R_1,D_1)\star(R_2,D_2)}(C). \tag{15}
\]
The filling in (14) is again zero-diagonal with sign entries off the
diagonal: for two distinct composite coordinates, exactly one of the
two summands is nonzero.

Reshaping the tuple indices in (9)--(10), or simply using (13) and
(15), proves
\[
\boxed{
\mathscr R_{R_2,D_2}\circ\mathscr R_{R_1,D_1}
=
\mathscr R_{(R_1,D_1)\star(R_2,D_2)}.
}
\tag{16}
\]
The scale update \(\eta\mapsto\eta/\sqrt s\) is essential for exact
associativity of the \(D\)-term.

At the boundary \(\eta=0\), the filling disappears at every fixed
renormalization depth, and the macroscopic action depends only on \(R\).
That order of limits does not justify discarding \(D\) when the depth
grows with the target order; (12) is the exact formulation that keeps
the two limits separate.

## 4. Why a minimizing fixed point does not follow

Let
\[
q_n=\min_C\Phi(\mathfrak P(C))
=\frac{2F(n)}{n^{3/2}}. \tag{17}
\]
If \(C_n\) is nearly minimizing, a lift only gives the one-sided
construction bound
\[
q_{ns}
\le
\Phi\!\left(
\mathscr R_{R,D}(n^{-1/2},\mathfrak P(C_n))
\right). \tag{18}
\]
Nothing in compactness or associativity says that the right side is
near \(q_n\), or even that the image profile is minimizing.

For the three-fibre pair
\[
R=J_3-2I_3,\qquad D=J_3-I_3,
\]
the exact calculation in `finite_fibre_renormalization.md` gives
\[
\frac{Q(C_6)}{6^{3/2}}
=0.6804138\ldots,\qquad
\frac{Q(\mathcal T(C_6))}{18^{3/2}}
=1.0213764\ldots. \tag{19}
\]
So this continuous semigroup action moves sharply away from the useful
upper-bound region.

Even an idempotent in a compact closure of the gadget semigroup would
not suffice:

1. an idempotent \(e^2=e\) need not absorb other profiles;
2. it need not lie in the closure of minimizing profiles;
3. an orbit may have idempotent limit maps while failing to converge;
4. a fixed point under one fibre size controls only its multiplicative
   scale tower, not all integer orders;
5. the family of powers is not equicontinuous in (5), because the
   \(d\)-th iterate reads level \(s^d r\).

The last item is the profile version of the changing-scale obstruction.

## 5. Concrete counterexamples

### 5.1 Equal scalar input, different renormalized scalar

Let \(C^-\) and \(C^\triangle\) be the order-three signings with edge
lists
\[
(-1,-1,-1),\qquad(-1,-1,1).
\]
They satisfy
\[
\Phi(\mathfrak P(C^-))
=\Phi(\mathfrak P(C^\triangle))
=\frac6{3^{3/2}}=\frac2{\sqrt3}. \tag{20}
\]
For the same three-fibre map,
\[
Q(\mathcal T(C^-))=24,\qquad
Q(\mathcal T(C^\triangle))=36. \tag{21}
\]
Therefore
\[
\Phi(\mathfrak P(\mathcal T(C^-)))=\frac89,\qquad
\Phi(\mathfrak P(\mathcal T(C^\triangle)))=\frac43. \tag{22}
\]
The scalar input objective neither determines nor contracts the scalar
output objective.  This is an exact realizable counterexample.

### 5.2 A realizable profile two-cycle

Let \(C_+=J_3-I_3\), and use the size-one fibre
\[
R=(-1),\qquad D=(0).
\]
Then \(\mathcal T_{R,D}(C_+)=-C_+\), and applying the same map twice
returns \(C_+\).  Hence
\[
\mathfrak P(C_+)\longleftrightarrow\mathfrak P(-C_+) \tag{23}
\]
is an exact period-two orbit of the continuous associative profile
action.

The one-sided continuous endpoint
\[
\Phi_+(\mathfrak P)=\sup_{(O,H)\in\mathcal P_1}H_{11}
\]
alternates between
\[
\Phi_+(\mathfrak P(C_+))=\frac2{\sqrt3},\qquad
\Phi_+(\mathfrak P(-C_+))=\frac2{3\sqrt3}. \tag{24}
\]
The absolute objective \(\Phi\) is equal on these two particular states,
as it must be under \(C\mapsto-C\).  The example is not offered as a
counterexample to the original limit; it is a counterexample to the
claim that compactness plus an associative profile action forces profile
or endpoint convergence.

### 5.3 The arity-shift countermodel

The exact map (12) has the coordinate pattern
\[
\text{output level }r
\quad\longleftarrow\quad
\text{input level }sr. \tag{25}
\]
The minimal topological model of this pattern is the Cantor space
\[
X=\{0,1\}^{\mathbb N_0},
\]
with the continuous shift \(S(\omega)_j=\omega_{j+1}\) and continuous
observable \(\phi(\omega)=\omega_0\).  For
\[
\omega=(0,1,0,1,\ldots),
\]
\[
\phi(S^d\omega)=0,1,0,1,\ldots . \tag{26}
\]
On this two-point orbit, \(S^2\) is the identity, hence an idempotent
map, but neither the orbit nor its observable converges.

This is an abstract countermodel, not a claimed realizable Boolean
profile.  Its role is precise: any argument that uses only product
compactness, continuity, associativity, and existence of an idempotent
is invalid.  A successful proof must use additional realizability or
extremality structure that excludes the shift behavior.

## 6. What theorem would actually suffice

The profile formalism turns the missing step into a clean target.  One
sufficient statement is:

> **Uniform absorbing recovery.**  There is an extremal limit profile
> \(\mathfrak P_*\) with \(\Phi(\mathfrak P_*)=\liminf q_n\) such that,
> for every sufficiently large target order \(m\), some valid
> finite-fibre completion of a microstate approximating
> \(\mathfrak P_*\) produces an order-\(m\) signing \(B_m\) with
> \[
> \Phi(\mathfrak P(B_m))
> \le\Phi(\mathfrak P_*)+o(1).
> \tag{27}
> \]

A stronger but easier-to-state route would be a contraction toward one
minimizing profile, uniform over a ratio-dense family of fibre gadgets.
The examples above show that neither contraction nor minimizing
invariance is formal.

Equation (27) is the quadratic-action quotient of the all-order recovery
theorem in `boolean_spectral_profile_compactification.md`.  It would
imply
\[
\limsup_m q_m\le\Phi(\mathfrak P_*)=\liminf_m q_m
\]
and hence existence of the desired limit.

## 7. Verdict

The all-level overlap/action profile (4) is enough—and no fixed tuple
truncation is enough—to make finite-fibre composition exact.  Its
renormalization map is associative and continuous.  This clarifies the
algebra but does not supply convergence:

- profile powers act as an arity shift;
- scalar \(Q\) loses decisive information;
- explicit lifts can move away from minimizers;
- compact semigroups can cycle even when idempotents exist.

The convergence problem remains exactly the absorbing-recovery problem
(27).  Any future renormalization proof must establish that theorem
from special structure of extremal Boolean signing profiles, not from
compactness alone.
