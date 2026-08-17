# Independent audit: flat Gram exposure and universal coordinate pins

**Verdict:** PASS for both drafts.  No mathematical repair is required.
There are two minor presentation recommendations:

1. in `universal_pin_cap_barrier.md`, define
   `||B||_(infinity->1)=max_(x,y in {+-1})|x^TBy|` explicitly, since the
   placement of the domain and codomain can otherwise be read in two ways;
2. keep the existing scope sentence that UP.1 assumes one common planted old
   state for every child.  It is strictly stronger than a contextual metric
   embedding with child-dependent active witnesses.

I ran

```text
experiments/verify_flat_gram_exposed_entropy.py
experiments/verify_universal_pin_cap_barrier.py
```

and both completed successfully.  The checks below are independent of those
finite regressions.

## 1. Hypercontractive sign balance

For a degree-two Boolean polynomial `P`, Bonami--Beckner at `q=4` gives

```math
||P||_4\le(q-1)^{2/2}||P||_2=3||P||_2.
```

The interpolation parameter in

```math
||P||_2\le ||P||_1^theta||P||_4^{1-theta}
```

is exactly `theta=1/3`, because
`1/2=theta+(1-theta)/4`.  Hence, writing
`sigma=||P||_2`,

```math
sigma\le ||P||_1^{1/3}(3sigma)^{2/3}
\quad\Longrightarrow\quad ||P||_1\ge sigma/9.
```

Homogeneity implies `E P=0`, so both positive and negative parts have mean
at least `sigma/18`.  Cauchy--Schwarz therefore gives

```math
Pr(P>0),Pr(P<0)\ge(1/18)^2=1/324.
```

Thus FE.1 and its constants are correct.  Distinct quadratic monomials are
distinct Fourier characters, so a nonzero hollow coefficient matrix really
does give a nonzero polynomial; there is no degenerate exception hidden in
the hypothesis.

## 2. Orientation and selector consequence

For nonzero `D`, one of `H_D` and `-H_D` has minimum exactly `-q(D)`.
Choosing that orientation `s` and a minimizer `p` gives

```math
sH_D(p)=sH_D(-p)=-q(D).
```

FE.1 applied to `sH_D` supplies at least `2^k/324` points with
`sH_D(x)>0`; the draft uses the weak inequality `>=0`, which can only enlarge
the set.  Every such point is at gap at least `q(D)` above both baseline
points.  FE.2 is therefore correct.

To invoke SC.5a under uniform approximation error `eta`, one takes

```math
a=q(D)-2eta>0.
```

This gives exactly

```math
log K\ge -log324+{(q(D)-2eta)^2\over2D_0^2}.
```

There is no omitted entropy term: `log|X|-k log2>=-log324`.

For the short-seed family, Theorem 21.26 supplies simultaneously

```math
q(D_(B,T))\ge c_0k^{3/2},\qquad c_0=\sqrt2/16,
\qquad ||D_(B,T)||_op\le16\sqrt k.
```

Thus FE.9 follows by weakening `q(D)` to `c_0k^(3/2)`.  Because this
reasoning is deterministic after the one seed event already controlling all
ordered pairs, the assertion that no second union bound is needed is also
correct.

If affine slopes have the form `Cy`, then
`||Cy||_2<=||C||_op sqrt(m)`.  For a complete sign `k` by `m` block,

```math
||C||_op\le\sqrt{||C||_1||C||_infinity}=\sqrt{km},
```

so `D_0^2<=km^2`.  Combining FE.9 with `K<=2^m` gives

```math
m log2+log324\ge {c_0^2k^2\over2m^2}
={k^2\over256m^2},
```

and hence the constant in FE.10.  Under
`||C||_op<=L sqrt(k)`, one instead has `D_0^2<=L^2km`, which yields

```math
m\ge {c_0\over\sqrt2L\sqrt{log2}}k-o(k)
={1\over16L\sqrt{log2}}k-o(k).
```

Thus FE.10--FE.11 and their fan-in/operator-norm dependence are correct.

The scope paragraph is essential and accurate.  If a fixed child term
`h(x)` remains outside the affine maximum, SC.5a applies to the residual
`f-h`.  In the displayed compiler, this residual is `-H_(A_T)`, not
`H_(A_B-A_T)`.  FE.9 therefore does not prove an unconditional lower bound
for every one-sided exact-sign compiler.  The draft does not overclaim that
step.

## 3. Absolute-value exposure

Parseval gives

```math
E H_D^2=\sum_(i<j)D_ij^2<2k^2
```

under `|D_ij|<=2`, so a point `p` with `|H_D(p)|<=sqrt(2)k` exists.  For a
maximizer `z` of `|H_D|` and a point `x` at Hamming distance `d`, symmetry of
`D` gives the sharper identity/bound

```math
|H_D(x)-H_D(z)|
={1\over2}|(x-z)^TD(x+z)|
\le2||D||_op\sqrt{d(k-d)}
\le2||D||_op\sqrt{dk}.
```

With `||D||_op<=L sqrt(k)` and
`delta=(c/(4L))^2`, this is at most
`ck^(3/2)/2<=q(D)/2`.  The Hamming-ball entropy estimate and
`delta=1/524288` for `c=sqrt(2)/16,L=16` are correct.  As the draft notes,
this positive-rate set has a linear entropy deficit, so it does not inherit
the full-strength signed selector lower bound.

## 4. Universal coordinate pin

The effective future is even:

```math
g(-x)=\max_y\{-x^TBy+H_C(y)\}
=\max_(y'=-y)\{x^TBy'+H_C(y')\}=g(x).
```

The same holds for every quadratic child.  Hence replacing `x` by `-x`
really does reduce its distance from `u` to the projective distance
`d<=floor(k/2)` without changing either side of UP.3.

Let `S` be the differing coordinates and put
`A_ij=-u_i u_j` on the cut `(S,S^c)`.  On each of its `d(k-d)` edges the
energy changes from `-1` at `u` to `+1` at `x`; all noncrossing monomials are
unchanged.  Therefore

```math
H_A(x)-H_A(u)=2d(k-d),
```

and UP.3 gives UP.4 with exactly the displayed `-eta`.  Choosing
`d=floor(k/2)` proves UP.5.

For two old states,

```math
g(x)-g(x')\le\max_y(x-x')^TBy
\le2||B||_(infinity->1).
```

Applying this to a pair realizing the oscillation gives

```math
||B||_(infinity->1)\ge d_*(k-d_*)-eta/2.
```

Finally, at fixed `(x,y)`, the parent values at `(x,y)` and `(x,-y)` are
`a+b` and `a-b`, with `b=x^TBy`; since
`max(|a+b|,|a-b|)>=|b|`, maximizing over `x,y` proves
`Q(P_A)>=||B||_(infinity->1)`.  Thus UP.6--UP.8 have the correct factors.

No sign assumption on `B,C` is needed for these inequalities.  Exact
hollowness and complete unit signs are needed only for the interpretation as
a complete-sign parent.  The theorem's quantifiers also match its
conclusion: the adversarial child is allowed to depend on the tested old
state because UP.3 is assumed for every child.  Conversely, the theorem says
nothing about a metric compiler whose active old witness depends on that
child.  This is exactly the boundary stated in the draft.

## 5. Final classification

- **FE.1--FE.4:** rigorous, with constants and orientations verified.
- **FE.9--FE.11:** rigorous conditional selector lower bounds; not an
  unconditional one-sided compiler obstruction.
- **UP.1:** rigorous universal-coordinate-pin obstruction with exact
  quadratic leading scale.
- **Needed edits:** no mathematical edits; only the optional explicit norm
  definition above.
