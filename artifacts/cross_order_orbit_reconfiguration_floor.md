# Orbit reconfiguration floor for balanced block compilers

Status: **proved zero-temperature method-class obstruction for actual
minimizing children**.

This note concerns the extensive transform

```math
b_n=M_n^{2/3}.
```

It proves that a balanced two-child compiler which works uniformly over
equivalent switching representatives cannot have an `o(n)` defect while its
bridge output changes by only `o(n^{3/2})` signs across that orbit.  The
bridges themselves are otherwise arbitrary; in particular, no rank,
template, randomness, or separate-channel hypothesis is imposed.

## Definitions

For a symmetric hollow signing `A` of order `r`, write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\operatorname{cap}(A)=\max_x|H_A(x)|.
```

Let `A` be an exact minimizer, so
`cap(A)=M_r`.  Its **extended switching orbit** is

```math
\mathcal O(A)=
\{\epsilon\,\operatorname{diag}(s)A\operatorname{diag}(s):
  \epsilon\in\{\pm1\},\ s\in\{\pm1\}^r\}.
```

Every member of `O(A)` is again an exact minimizer.  A deterministic
balanced block compiler on this orbit is an arbitrary map

```math
\mathfrak C:\mathcal O(A)\times\mathcal O(A)
   \longrightarrow \{\pm1\}^{r\times r}.
```

For children `X,Y`, it returns the order-`2r` signing

```math
B(X,Y)=
\begin{pmatrix}
X&\mathfrak C(X,Y)\\
\mathfrak C(X,Y)^{\mathsf T}&Y
\end{pmatrix}.                                      \tag{1}
```

Its worst output defect on the orbit is

```math
E_r^{\mathfrak C}=
\max_{X,Y\in\mathcal O(A)}
\left\{\operatorname{cap}(B(X,Y))^{2/3}-2M_r^{2/3}\right\}. \tag{1a}
```

This is a compiler-output quantity, not
`M_(2r)^(2/3)-2M_r^(2/3)`: minimizing over all order-`2r` signings goes in
the opposite direction.  The theorem rules out a uniform compiler
guarantee, not the existence of one favorable output in its orbit.

For sign matrices `C,C'`, let

```math
d_\pm(C,C')=\min\{d_{\rm Ham}(C,C'),d_{\rm Ham}(C,-C')\}.
```

The global sign is quotiented out because replacing `C` by `-C` is absorbed
by flipping all spins in one child block.  Define the projective Hamming
diameter of the compiler output by

```math
h_r=\max_{X,Y,X',Y'\in\mathcal O(A)}
d_\pm\bigl(\mathfrak C(X,Y),\mathfrak C(X',Y')\bigr). \tag{2}
```

Finally, for an `r` by `r` bridge put

```math
R(C)=\max_{u,v\in\{\pm1\}^r}|u^{\mathsf T}Cv|,
\qquad
\kappa_r={\mathbb E|S_r|\over\sqrt r},
```

where `S_r` is a sum of `r` independent fair signs.

## Theorem: direct defect and reconfiguration bounds

Suppose a compiler satisfies the uniform `b`-scale recurrence certificate

```math
\operatorname{cap}(B(X,Y))^{2/3}
\le 2M_r^{2/3}+e_r
\quad\hbox{for every }X,Y\in\mathcal O(A).          \tag{P}
```

Equivalently, the premise is `E_r^(mathfrak C)<=e_r`.

Put

```math
a_r={M_r\over r^{3/2}},\qquad
\eta_r={h_r\over r^{3/2}}.
```

Then the premise `(P)` has the following two explicit implications:

```math
\boxed{
{e_r\over r}\ \ge
\left[2a_r+(\kappa_r-2\eta_r)_+\right]^{2/3}
-2a_r^{2/3}.}                                      \tag{3}
```

and, equivalently in the direction useful for a target defect,

```math
\boxed{
h_r\ \ge {1\over2}
\left[
r\,\mathbb E|S_r|+2M_r
 -(2M_r^{2/3}+e_r)^{3/2}
\right]_+.}                                        \tag{4}
```

Thus `(P) => e_r` is quantitatively bounded by (3), and
`(P)` together with a proposed value of `e_r` forces the bridge motion (4).

Using the proved upper frontier

```math
a_r\le {1\over2}+o(1)                              \tag{5}
```

and `kappa_r -> sqrt(2/pi)`, equations (3)--(4) give two concrete
consequences.

First, orbit stability `h_r=o(r^{3/2})` forces a positive linear defect:

```math
\liminf_{r\to\infty}{e_r\over r}
\ge
\left(1+\sqrt{2/\pi}\right)^{2/3}-2^{1/3}
=0.218646607120494\ldots .                           \tag{6}
```

Second, the desired `e_r=o(r)` forces macroscopic bridge
reconfiguration:

```math
\boxed{
h_r\ge
\left({\sqrt{2/\pi}-(2^{3/2}-2)/2\over2}-o(1)\right)
r^{3/2}
=\bigl(0.191835499214885\ldots-o(1)\bigr)r^{3/2}.}  \tag{7}
```

For an all-order bridge check not using the central-limit asymptotic, the
Khintchine bound `kappa_r>=1/sqrt(2)` may replace the bridge constant.  If
one also has `a_r<=1/2` at the order in question, the corresponding
constants in (6) and (7) are, respectively,

```math
\left(1+{1\over\sqrt2}\right)^{2/3}-2^{1/3}
=0.168448089030265\ldots
```

and

```math
{1/\sqrt2-(2^{3/2}-2)/2\over2}
=0.146446609406726\ldots .                           \tag{8}
```

The `o(1)` in (7) only comes from the cap upper frontier and the proposed
`e_r=o(r)`; the bridge lower bound itself is exact at every order.

## Proof

### 1. Exact block parity

For arbitrary child signings `X,Y`, bridge `C`, and Boolean spins `x,y`, put

```math
U=H_X(x)+H_Y(y),\qquad V=x^{\mathsf T}Cy.
```

Replacing `y` by `-y` preserves `U` and reverses `V`.  Therefore

```math
\max\{|U+V|,|U-V|\}=|U|+|V|.                       \tag{9}
```

In particular, if both child energies at `x,y` equal `M_r`, the parent cap
is at least `2M_r+|x^T C y|`.  This is an exact joint statement, not a
triangle-inequality payment of separately maximized channels.

### 2. A reference bridge has a large cut

For every sign bridge `C`, first maximize over `v` and then average over a
uniform `u`:

```math
\begin{aligned}
R(C)
&=\max_u\sum_{j=1}^r |(C^{\mathsf T}u)_j|\\
&\ge \mathbb E_u\sum_{j=1}^r
 \left|\sum_{i=1}^r C_{ij}u_i\right|\\
&=r\,\mathbb E|S_r|
=\kappa_r r^{3/2}.                                  \tag{10}
\end{aligned}
```

The last equality holds because multiplying independent fair signs by a
fixed sign column does not change their law.  Explicitly,

```math
\mathbb E|S_{2k}|={2k{2k\choose k}\over2^{2k}},
\qquad
\mathbb E|S_{2k+1}|={(2k+1){2k\choose k}\over2^{2k}},
```

so `kappa_r -> sqrt(2/pi)`.  The all-order estimate
`kappa_r>=1/sqrt(2)` is the `p=1` Khintchine inequality.

### 3. Align the actual child ground states

Choose any reference output `C_0` of the compiler, and choose `u,v` with

```math
|u^{\mathsf T}C_0v|=R(C_0).                         \tag{11}
```

Take a ground state `x_*` of `A`, say
`H_A(x_*)=\sigma M_r` with `\sigma\in\{\pm1\}`.  With
`s=x_*\odot u`, the orbit representative

```math
X=\sigma\,\operatorname{diag}(s)A\operatorname{diag}(s)
```

satisfies `H_X(u)=M_r`.  An independent identical choice gives
`Y in O(A)` with `H_Y(v)=M_r`.  This is why the theorem uses the extended
orbit: switching places a ground state at a prescribed cut vector, and the
global sign orients its energy positively.

Let `C'=\mathfrak C(X,Y)`.  By the diameter definition, there is a
`\tau\in\{\pm1\}` such that `C'` and `\tau C_0` differ in at most `h_r`
entries.  Each changed sign alters a fixed bilinear value by at most `2`, so

```math
|u^{\mathsf T}C'v|
\ge R(C_0)-2h_r.                                    \tag{12}
```

Combining (9), (10), and (12) yields

```math
\operatorname{cap}(B(X,Y))
\ge 2M_r+\bigl(r\mathbb E|S_r|-2h_r\bigr)_+.       \tag{13}
```

Applying `(P)` to these particular, still exactly minimizing,
representatives and taking the `2/3` power proves (3).

Conversely, (13) and `(P)` imply

```math
\bigl(r\mathbb E|S_r|-2h_r\bigr)_+
\le (2M_r^{2/3}+e_r)^{3/2}-2M_r.
```

Whether the expression before the positive part is positive or not, this
rearranges to (4).

### 4. Constants and asymptotics

For fixed `q>0`, the function

```math
f(a,q)=(2a+q)^{2/3}-2a^{2/3}
```

is decreasing in `a` and increasing in `q`.  Taking
`eta_r=o(1)`, `a_r<=1/2+o(1)`, and
`kappa_r->sqrt(2/pi)` in (3) proves (6).

If `e_r=o(r)`, divide (4) by `r^(3/2)` and use

```math
{(2M_r^{2/3}+e_r)^{3/2}\over r^{3/2}}
=2^{3/2}a_r+o(1).
```

Together with (5) and (10), this proves (7).  Direct numerical evaluation
gives the displayed constants; (8) follows in the same way from the
all-order Khintchine bound.

## Relation to earlier obstructions

This statement has different quantifiers from the existing artifacts:

- The switching-template obstruction restricts the structural complexity
  of a parent or its cross blocks.  Here every bridge output can be an
  arbitrary full-complexity sign matrix.
- The universal edit-net obstruction asks for a catalogue approximating
  every signing.  Here only the outputs over one actual minimizer's extended
  switching orbit are compared; no covering of the ambient bridge cube is
  assumed.
- The finite-temperature selector-transport obstruction measures how far a
  bridge law moves from the uniform law.  Here the result is deterministic
  and zero-temperature, and measures how much the compiler must respond to
  equivalent presentations of its actual children.

The new content is the forced response across the child orbit: a uniform
`o(r)` recurrence must change at least
`(0.191835499...-o(1))r^(3/2)` bridge signs (modulo global bridge sign)
between some two equivalent child inputs.

## Scope

This is a method-class obstruction, not a proof or disproof of the desired
cross-order recurrence.

- It applies to balanced block compilers which retain the two supplied
  children as principal blocks and claim a uniform guarantee over their
  extended switching orbit.
- That uniformity is essential.  A non-equivariant argument allowed to pick
  only one favorable representative is not ruled out.
- The compiler is deterministic.  A randomized version would require the
  diameter and recurrence guarantee supportwise; average-law statements are
  not supplied by this proof.
- A switching-equivariant or otherwise global compiler can change
  `Omega(r^(3/2))`--and typically `Theta(r^2)`--of the `r^2` bridge entries.
  Such a compiler escapes the theorem.
- No conclusion is drawn about wholesale parent constructions that do not
  retain the two children, unbalanced splits, or finite-temperature
  functionals.

Thus the theorem rules out orbit-stable and local-repair cancellation, while
leaving precisely the globally responsive bridge architectures open.
