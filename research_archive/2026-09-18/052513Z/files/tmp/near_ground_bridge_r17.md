# Wave 17: two-near-state polarization gives raw shores, not selection

Status: the identities and inequalities below are **exact**.  Every finite
count is exhaustively checked by

```bash
/home/math/quadra/.venv/bin/python tmp/check_near_ground_bridge_r17.py
```

The positive theorem produces raw decrement-tolled buckets at one root cut.
It does not produce an allocation, an endpoint tree, a descendant resource,
or a causal temporal transport.  The finite audits sharply falsify selection
from the old one-state response alone.

## 1. Exact opposite-state polarization

Let `M` be a symmetric zero-diagonal matrix, put `q=Q(M)`, and take two
oriented Boolean states with opposite orientation labels

```math
\omega=(\sigma,z),\qquad \bar\omega=(-\sigma,w).
```

They need not be grounds of `M`.  Write their full-matrix deficits as

```math
\alpha=q-\sigma z^{\mathsf T}Mz,
\qquad
\beta=q+\sigma w^{\mathsf T}Mw.
```

Let `S={i:z_i=w_i}` and `T={i:z_i=-w_i}`.  In the `z` gauge define

```math
h_X=\sigma z_X^{\mathsf T}M[X]z_X\quad(X=S,T),
\qquad
c=\sigma z_S^{\mathsf T}M[S,T]z_T,
\qquad
b=\sigma z^{\mathsf T}Mw.
```

Expansion at the four corners gives

```math
\begin{aligned}
q-\alpha&=h_S+h_T+2c,\\
q-\beta&=-h_S-h_T+2c,\\
b&=h_S-h_T.
\end{aligned}
```

Therefore, with `H=(\beta-\alpha)/2`,

```math
\boxed{
4c=2q-\alpha-\beta,
\qquad h_S+h_T=H,
\qquad h_S-h_T=b.
}
\tag{R17.1}
```

The orientation labels are opposite for every matrix being evaluated.  When
the pair is selected below, the states are specifically opposite-oriented
near grounds of `C`; they need not be near grounds, let alone exact grounds,
of `M=A`.

## 2. Exact four-bucket identities

For a retained shore `X` and its complement `Y`, put

```math
L_X=\lVert M[Y,X]z_X\rVert_1,
\qquad
\ell_X=L_X-c.
```

Since `L_X\ge |c|`, both `\ell_S,\ell_T` are nonnegative.  The raw
decrement-tolled all-successor residual (10.532) at this root split is

```math
r_X^\tau=
\left[2L_X+\tau z_X^{\mathsf T}M[X]z_X-q\right]_+,
\qquad \tau\in\{\sigma,-\sigma\}.
```

Substituting (R17.1), separately for each shore and orientation, gives the
four exact formulas

```math
\boxed{
\begin{aligned}
r_S^\sigma&=[2\ell_S-h_T-\alpha]_+,&
r_S^{-\sigma}&=[2\ell_S+h_T-\beta]_+,\\
r_T^\sigma&=[2\ell_T-h_S-\alpha]_+,&
r_T^{-\sigma}&=[2\ell_T+h_S-\beta]_+.
\end{aligned}
}
\tag{R17.2}
```

For example,

```math
2L_S+h_S-q
=2\ell_S+(2c+h_S-q)
=2\ell_S-h_T-\alpha;
```

the other three lines are identical substitutions, using the second state
for the opposite orientation.  Thus no ground or local-optimality assumption
is hidden in (R17.2).

Let

```math
\mathscr S_M(\omega,\bar\omega)
=r_S^\sigma+r_S^{-\sigma}+r_T^\sigma+r_T^{-\sigma}.
```

Dropping the nonnegative `2\ell_X` terms and using

```math
[-x-\alpha]_++[x-\beta]_+
=\left[\left|x-H\right|-\frac{\alpha+\beta}{2}\right]_+
```

yields

```math
\boxed{
\begin{aligned}
\mathscr S_M(\omega,\bar\omega)
&\ge \frac12\left(
[|b-H|-(\alpha+\beta)]_+
+[|b+H|-(\alpha+\beta)]_+
\right)\\
&\ge [|b|-\alpha-\beta]_+.
\end{aligned}
}
\tag{R17.3}
```

The last step is convexity of `x\mapsto[|x|-(\alpha+\beta)]_+` at the
two points `b-H,b+H`.  If the two states are opposite exact grounds of a
balanced `M`, the endpoint exposure identity gives `\ell_S=\ell_T=0` and
(R17.3) reduces to the known exact endpoint formula
`\mathscr S_M=|z^{\mathsf T}Mw|`.

## 3. Specialization to the deficit-corrected active response

Now let `A=C+D`, put

```math
q=Q(A),\qquad q_C=Q(C),\qquad \Delta=q-q_C,
```

and select

```math
\omega=(\sigma,z),\quad \bar\omega=(-\sigma,w)
\quad\text{with}\quad
u=\delta_C(\omega)\le t,\quad
v=\delta_C(\bar\omega)\le t.
```

Define the two actual deficit-corrected response scores

```math
s(\omega)=e_D(\omega)-u,
\qquad
s(\bar\omega)=e_D(\bar\omega)-v.
```

These are precisely the quantities whose one-state maximum is
`\mathcal R_{C,t}(D)`.  Direct cancellation gives

```math
\boxed{
\alpha=\Delta-s(\omega),
\qquad
\beta=\Delta-s(\bar\omega).
}
\tag{R17.4}
```

Polarization of `C` alone also gives

```math
c_C=\frac{q_C}{2}-\frac{u+v}{4},
\qquad
h_{C,S}+h_{C,T}=\frac{v-u}{2}.
```

More strikingly, the sum of the two response scores becomes the full-matrix
cross exposure of the pair cut:

```math
\boxed{
c_A=\frac{q_C}{2}
+\frac{s(\omega)+s(\bar\omega)}{4}.
}
\tag{R17.5}
```

Equations (R17.3)--(R17.4) prove the genuine paired raw-shore theorem

```math
\boxed{
\mathscr S_A(\omega,\bar\omega)
\ge
\left[
|\sigma z^{\mathsf T}Aw|-2\Delta
+s(\omega)+s(\bar\omega)
\right]_+.
}
\tag{R17.6}
```

Equivalently, either this pair supplies at least `L` raw root residual, or

```math
|z^{\mathsf T}Aw|
<L+2\Delta-s(\omega)-s(\bar\omega).
\tag{R17.7}
```

Thus the precise missing datum is the cross-Gram of two *actual global*
near-`C` states, after their two full `A` deficits are paid.

For clarity, define the orientation-specific responses

```math
\mathcal R_{C,t}^{\sigma}(D)
=\max_{\delta_C(\sigma,z)\le t}
\{e_D(\sigma,z)-\delta_C(\sigma,z)\}.
```

The old response is
`\mathcal R_{C,t}=\max_\sigma\mathcal R_{C,t}^\sigma`, and
`\Gamma=[-\mathcal R_{C,t}]_+`.  A pair of orientationwise maximizers has
minimum possible total full-matrix deficit

```math
\Xi_{C,t}(D)
=2\Delta-
\mathcal R_{C,t}^{+}(D)-\mathcal R_{C,t}^{-}(D).
\tag{R17.8}
```

Even (R17.8), which is stronger than the old `\mathcal R/\Gamma`, does not
control their mutual cross-Gram.  The useful new paired functional is

```math
\mathcal B_{C,t}(A)
=\max_{\substack{\omega=(\sigma,z),\ \bar\omega=(-\sigma,w)\\
\delta_C(\omega),\delta_C(\bar\omega)\le t}}
\left[|z^{\mathsf T}Aw|-\delta_A(\omega)-\delta_A(\bar\omega)\right]_+.
\tag{R17.9}
```

For nondegenerate agreement/disagreement cuts, (R17.3) proves that the best
raw root-shore total is at least `\mathcal B_{C,t}(A)`.  This maximization is
over pairs of actual correlated global states; it never forms Cartesian
products of blockwise states.

There is also an unavoidable orientation-margin wall.  If
`P(C)=\max_z z^TCz`, `N(C)=-\min_z z^TCz`, and
`I(C)=|P(C)-N(C)|`, the minimum deficit in orientations `+` and `-` is,
respectively, `q_C-P(C)` and `q_C-N(C)`.  Hence

```math
\boxed{
\text{both opposite orientation classes meet }\mathcal A_C(t)
\quad\Longleftrightarrow\quad t\ge I(C).
}
\tag{R17.10}
```

The scalar `\Gamma` records neither this two-orientation coverage nor the
cross-Gram in (R17.9).

## 4. `A_8`: a profitable arbitrary cut, but response maximizers are neutral

For the partition

```math
\{0,3,4,5\}\mid\{1,2,6,7\},
```

the cross mosaic has `P(C)=N(C)=16`.  At `t=0` there are 12 projective
`C`-grounds in each orientation.  In either orientation their response-score
histogram is

```text
-4: 4 states,  0: 4 states,  4: 4 states.
```

Thus `\mathcal R^+=\mathcal R^-=4`, `\Gamma=0`, and `\Delta=4`.
The four maximizing states in each orientation are exact `A_8` grounds.
All 16 opposite-orientation maximizing pairs have

```math
\alpha=\beta=0,
\qquad z^{\mathsf T}A_8w=0,
\qquad \mathscr S_{A_8}=0.
\tag{R17.11}
```

These are exactly the old neutral endpoint pairs.  Therefore even knowing
both orientation-specific response maxima does not select raw resource.

The full search has `12\cdot12=144` projective opposite-orientation pairs
(288 if the two choices of `\sigma` are both recorded).  Its paired score
and actual raw total both have maximum eight.  Exactly four projective pairs
attain it.  For every attaining pair,

```math
u=v=0,
\quad s(\omega)=s(\bar\omega)=0,
\quad \alpha=\beta=4,
\quad |z^{\mathsf T}A_8w|=16,
\quad \mathscr S_{A_8}=8.
\tag{R17.12}
```

Each resulting `4+4` cut has `Q(A_8[S])=Q(A_8[T])=q_4=8`; its four raw
bucket residuals are a permutation of `(4,4,0,0)`.  Hence two actual exact
`C`-grounds do uncover a new raw root-shore resource behind the endpoint wall.
But they deliberately sacrifice four response units on each side, are not
`A_8` endpoint pairs, and their cuts do not occur in the existing endpoint
tree or its conserved allocation.

## 5. All forty `A_9` completions

For the `6+1+1+1` cross mosaic,

```math
P(C)=18,
\qquad N(C)=Q(C)=22,
\qquad I(C)=4.
```

Every one of the forty norm-24 completions `H_G=C+G` has the following exact
census.

At `t=0`, there are three `C`-grounds in orientation `-` and none in
orientation `+`.  Thus no opposite pair exists.  Here
`\mathcal R=-6` and `\Gamma=6`.  This realizes (R17.10) sharply.

At `t=4`, there are 15 and 4 near states in the two orientations, hence 60
projective pairs per completion (120 directed records).  In both orientations
`\mathcal R^\sigma=2=\Delta`, so every response maximizer used below is an
exact `H_G` ground.  Twenty-four completions have 12 such projective pairs and
sixteen have 24.  Across all forty, their exact `(|b|,\mathscr S)` histogram is

| `(|b|, raw)` | count |
|---:|---:|
| `(0,0)` | 336 |
| `(4,4)` | 192 |
| `(8,8)` | 144 |

The maximum over *all* 60 near-state pairs is raw residual eight in every
completion.  This is strictly below the compulsory common-mosaic response
rise twelve.  Thus a coefficient-one theorem charging that rise to a paired
raw shore at the first usable margin is false on all forty completions.

At `t=8`, there are 33 and 22 near states, hence 726 projective pairs per
completion (1,452 directed records).  The response-maximizing subset has 108
pairs in twenty-four completions and 156 in sixteen.  Across all forty it is
the full 5,088 endpoint-pair census:

| `(|b|, raw)` | count |
|---:|---:|
| `(0,0)` | 1,944 |
| `(4,4)` | 1,152 |
| `(8,8)` | 1,896 |
| `(16,16)` | 96 |

The best paired raw total is now sixteen.  The larger margin can therefore
find enough raw root resource, but by this point the active symmetric gauge
has already reached `Q(G)=10`, as recorded in (10.597).  This is a finite
possibility, not an asymptotic partition theorem or a temporal transport.

## 6. Scoped conclusion

The exact positive result is (R17.6): two opposite-oriented actual near-`C`
states convert their **deficit-adjusted mutual `A` cross-Gram** into raw
decrement-tolled root-shore residual.  It is the first version of this route
that passes through the `A_8` endpoint-neutral wall, by changing to a
non-endpoint cut.

The exact selection obstruction is independent and decisive:

1. `t` must first cover the cross-mosaic imbalance `I(C)`;
2. the old one-state `\mathcal R/\Gamma` does not provide the other
   orientation;
3. even the two orientation-specific maxima do not control their mutual
   cross-Gram (`A_8`, equation (R17.11)); and
4. at the first usable `A_9` margin, all forty completions have at most eight
   paired raw units against a response rise of twelve.

No claim here concerns descendant raw residuals.  No raw bucket is assigned
the existing conserved endpoint allocation.  No resources from different
pairs are summed.  No internal response excess is reclassified as terminal
credit.  Finally, no temporal prefix, suffix, or antichain inequality follows
from a static pair: a new theorem would still have to choose the non-endpoint
cut in the actual peeling process and route its raw bucket through an honest
allocation and the causal Hall cuts of (10.543).
