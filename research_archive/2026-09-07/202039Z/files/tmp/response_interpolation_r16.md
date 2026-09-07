# A persistent-ground wall for response differentiation

Status: exact algebra and exhaustive finite claims independently checked by
`tmp/audit_response_interpolation_r16.py`.

Let `A=C+D` and `H=C+G` be two symmetric matrices on the same block mosaic,
and suppose

```math
Q(A)=Q(H)=q.
```

For `0\le t\le1`, put `A_t=C+(1-t)D+tG`.  Convexity of `Q` gives
`Q(A_t)\le q`.  If one signed Boolean state `(\sigma,z)` is a common ground,

```math
\sigma z^{\mathsf T}Az
=\sigma z^{\mathsf T}Hz=q,
```

then its energy is identically `q` along the segment.  Hence

```math
\boxed{Q(A_t)=q\quad(0\le t\le1).}
```

Now suppose `D=\bigoplus_iD_i`, `G=\bigoplus_iG_i`, set
`S=\sum_iQ(D_i)`, `B=\sum_iQ(G_i)`, and use the linearly interpolated
deficit baseline.  Its common-mosaic response is

```math
\Phi_t=Q(A_t)-[(1-t)S+tB].
```

On a persistent-ground segment,

```math
\boxed{\Phi_t=\Phi_0+t(S-B).}
```

Thus a response increase can be entirely a baseline effect while the global
support function has zero variation in every direction along the segment.
Any response-to-temporal proposal based only on integrating derivatives or
subgradients of `Q(A_t)` misses this effect.

For the `A_9` capture partition `6+1+1+1` from (10.579), there are exactly 40
order-six minimizing completions with hybrid norm 24.  Every one has at least
four signed projective states that are simultaneous `+24` grounds of the
original and hybrid matrices.  The exact number ranges from 4 to 13.  At
every persistent ground the signed energy profile is one of

```math
(C,D,G)=(14,10,10)\quad\text{or}\quad(18,6,6).
```

Thus the old and new internal energies agree at the persistent state.  The
old/new deficit pairs are respectively `(12,0)` or `(16,4)`: in either case
the deficit falls by exactly `X=12`, without any change in the state energy.

Consequently every one of the 40 interpolation segments has `Q(A_t)=24`,
while

```math
\Phi_t=2+12t.
```

This does not rule out a bridge that explicitly transports the falling local
baseline or the persistent ground's block deficits.  It rules out treating
the compulsory response increase as work done by the global norm along a
replacement interpolation.
