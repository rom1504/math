# Main-agent Wave 40 note: endpoint cost cannot universally control isotonic defect

This is an abstract two-bit exponential-family wall, not an actual signing or
restriction counterexample.

Let `nu` be uniform on `{00,01,10,11}`, set

```math
g_L(00)=0,
\qquad
g_L(01)=g_L(10)=g_L(11)=L.
```

and normalize the endpoint likelihood as
`f_L=e^{g_L}/Z_L`, `Z_L=E_nu e^{g_L}`.

Under `mu_s proportional to nu exp(s g_L)`, put

```math
a_s=\lambda_s\{0\}=\frac1{1+3e^{sL}},
\qquad
v_s=L^2\frac{e^{sL}}{(1+e^{sL})^2}.
```

For either bit, only the context in which the other bit is zero has a
nonzero score gap.  Thus, in state order `00,01,10,11`,

```math
V_s=(2v_s,v_s,v_s,0).
```

Conditioning on the two score levels gives the exact decreasing regression

```math
m_s(0)=2v_s,
\qquad
m_s(L)=\frac23v_s.
```

For the mixture `bar(lambda)_s=(lambda_s+lambda_1)/2`, let
`w_s=(a_s+a_1)/2` be its low-level mass.  Two-point isotonic regression pools
the decreasing pair, so

```math
\overline{R}_s
=w_s(1-w_s)\left(\frac43v_s\right)^2.             (W1)
```

Consequently

```math
K_L
=\frac43\int_0^1
\sqrt{s(1-s)w_s(1-w_s)}\,v_s\,ds.
```

With `u=sL`, dominated convergence after division by `sqrt(L)` gives

```math
\frac{K_L}{\sqrt L}\longrightarrow
\frac43\int_0^\infty
\frac{e^u}{(1+e^u)^2}
\sqrt{\frac{u}{2(1+3e^u)}
\left(1-\frac1{2(1+3e^u)}\right)}\,du
=:c_*>0.                                               (W2)
```

Hence `K_L^2~c_*^2 L`.  By contrast, the endpoint conditional KL receives,
for each bit, mass `(1+e^L)/(1+3e^L)` on the sole active context and binary
KL `D(Ber(e^L/(1+e^L)) || Ber(1/2))`.  Therefore

```math
C_V(1)\longrightarrow\frac23\log2.                    (W3)
```

The parent entropy is also bounded (it tends to `log(4/3)`).  Thus no
universal theorem of the form `K^2=O(C_V(1))`, or more generally one whose
right side stays bounded whenever endpoint cost and parent entropy stay
bounded (such as `K^2=O(C_V(1)+mathscr H)`), holds for arbitrary cube
exponential families.  A positive estimate of (10.1096) needs additional
structure beyond these endpoint scalars; this two-bit example alone does not
identify exact minimality as the unique possible input.

The numerical checker is `tmp/harmonic_abstract_isotonic_wall_r40_main.py`.
