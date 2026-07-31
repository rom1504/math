# Maximum conference excess is not a same-switch quadratic witness

Status: exact normalization/mapping audit. The source is Momihara--Suda,
[Conference matrices with maximum excess and two-intersection sets](https://math.colgate.edu/~integers/r30/r30.pdf),
*Integers* 17 (2017), A30 (also arXiv:1611.01305).

## 1. The two optimization orbits are different

Let `W` be their symmetric Paley conference matrix of order

```math
N=q+1=4m^2+2,qquad WW^{\mathsf T}=qI.
```

For Boolean row and column signs `r,c`, independently switching rows and
columns produces

```math
W'=D_rWD_c,qquad E(W')=\mathbf1^{\mathsf T}W'\mathbf1
                  =r^{\mathsf T}Wc.
```

This is a **bilinear** Boolean optimization. By contrast, the quadratic
energy of the symmetric edge signing is

```math
Q_W(x)=\sum_{i<j}W_{ij}x_ix_j=\frac12x^{\mathsf T}Wx,
```

and its allowed switching orbit is `D_x W D_x`, using the same sign vector
on both sides. Independent row/column equivalence is therefore strictly
broader than the equivalence relevant to `M_N`.

The distinction is present in the paper's construction, not merely a logical
possibility. If `D` is its two-intersection set, the proof of Theorem 2 uses

```text
c_i = -1  iff i is in D,
r_i = -1  iff i is in D^*_alpha,
```

on the `F_q` coordinates (both signs are `+1` at the bordered coordinate).
Here

```math
|D|=2m^2-m+1,qquad |D^*_{\alpha}|=2m^2-m.
```

Thus `r` and `c` cannot be equal. They cannot differ by a global sign either,
as their bordered signs agree and their finite negative-set sizes differ.
Indeed `W'` is not symmetric: symmetry at the bordered entries would already
require `r_i=c_i` for every finite coordinate.

## 2. Exact maximum-excess normalization

For `N=4m^2+2`, Proposition 7 uses the odd integer `2m-1`, since

```math
2m-1\leq \sqrt{N-1}<2m+1.
```

Its equality value is

```math
E(W')
=\frac{N\big((2m-1)^2+2(2m-1)+N-1\big)}{2(2m)}
=2mN.
```

Theorem 2 makes every row sum of `W'` equal to `2m-1` or `2m+1`;
the conference norm forces exactly `N/2` rows of each kind. Hence its theorem
provides an exact independent-sign witness

```math
r^{\mathsf T}Wc=2mN.
```

It does **not** provide an `x` with `x^{\mathsf T}Wx=2mN`. Polarization does
not repair this: `r+c` and `r-c` are ternary vectors with zeros, not Boolean
states. Nor can one use `W'` as an edge signing, because `W'` is nonsymmetric
and its symmetric part need not have off-diagonal entries in `{\pm1}`.
Consequently no bound on our same-switch quadratic cap, and in particular no
landing statement for `M_N`, follows from maximum excess alone.

There is a universal but lossy conversion that makes the limitation
quantitative. For Boolean `x,y`, partition the coordinates into
`I={i:x_i=y_i}` and `J={i:x_i=-y_i}`. Symmetry cancels the two cross terms:

```math
x^{\mathsf T}Wy
=2\{H_{W[I]}(x_I)-H_{W[J]}(x_J)\}.
```

Randomly completing either partial spin proves
`cap(W)>=max(|H_(W[I])|,|H_(W[J])|)`, and hence

```math
\operatorname{cap}(W)
\ge {1\over4}\max_{x,y}|x^{\mathsf T}Wy|.
```

The theorem's independent excess therefore yields only `cap(W)>=mN/2`,
asymptotic constant `1/4`. This is weaker than the universal
top-eigenspace-Gaussian constant `1/pi` and does not approach the required
`1/2` landing scale.

There is also an exact parity witness to the mismatch. A symmetric signing at
this order has an odd number of unordered edges, so every quadratic energy is
odd. Half of the independent excess is `mN=2m(2m^2+1)`, which is even.
Thus no same switch can attain the paper's equality value; the local-field
identity instead gives the symmetric upper bound `cap(W)<=mN-1`.

## 3. Relation to the two-fiber/self-indexed ASDS construction

The common part is only the starting Paley conference matrix at the same
order. The paper's incidence structure has point set the additive group of
`F_q`; its blocks are translates of the nonzero squares, and its set

```math
D\subset F_q,qquad |D|=2m^2-m+1,
```

has intersection sizes `m^2-m` and `m^2` with those blocks. This is not the
self-indexed ASDS datum

```math
P,R\subset Z_s,qquad s=2m^2+1,qquad
|P|=m^2-m,quad |R|=m^2,
```

with

```math
N_P(h)+N_R(h)+1_P(h)=m^2-m\quad(h\ne0).
```

The ambient groups, set sizes, and defining correlations differ.

After a nonsplit-torus gauge, the *unswitched* Paley conference can be put in
the two-fiber cyclic form and therefore gives the already identified Paley
self-indexed ASDS branch. Momihara--Suda then apply different coordinate
signs on its two sides. Those switches generally destroy both symmetry and
the common cyclic action. Their two-intersection theorem therefore neither
constructs a new self-indexed ASDS pair nor supplies a non-Paley two-fiber
class. It is an exact statement about an independently switched matrix in the
same broad row/column equivalence class as Paley.

## 4. Research judgment

The two-intersection set is a genuine, explicit certificate for a nearly
spectral **bilinear** excess. It may be useful only if one proves an additional
rounding theorem converting its particular pair `(r,c)` into a same-sign
Boolean vector with quantitatively controlled loss. No such conversion is in
the cited theorem, and a generic bilinear-to-quadratic reformulation is not a
landing mechanism. It should therefore not be counted as progress on the
structured-family landing obligation without a new same-switch lemma.
