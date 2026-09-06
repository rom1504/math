### Wave 29, Route 3: matched endpoint limits and the parent common mode

**Outcome.** Two endpoint expansions are **Verified**, including all
normalization factors, by tmp/matched_hellinger_r29.py. At exactly matched
temperature \(\gamma=\beta\):

1. the \(\beta\to\infty\) slope of the parent entropy is an explicit mean
   best-extension deficit of maximal child grounds;
2. the first nonzero \(\beta\to0\) common component of the parent likelihood
   is exactly the centered parent row square.

Neither result proves either uniform estimate in (10.800), and the finite
audits below do not falsify them. They identify what a minimizer-specific
comparison must add.

#### 1. Exact matched likelihood

Retain (10.797)--(10.813), and write

~~~math
L_\beta(S,d)=C_\beta e^{-\ell_S(d[S])},
\qquad
f_\beta(d)=\mathbb E_{S\sim U_m}L_\beta(S,d).
~~~

At \(\gamma=\beta\), \(\ell_S=\log K_{\beta,S}\), so

~~~math
\boxed{
L_\beta(S,d)=\frac{C_\beta}{K_{\beta,S}(d[S])},
\qquad
C_\beta=\frac{Z_A(\beta)}{Z_0(\beta)}.
}
\tag{R29.1}
~~~

The normalized parent endpoint law is \(P_\beta=f_\beta\nu_\beta\), and

~~~math
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=D(P_\beta\Vert\nu_\beta).
\tag{R29.2}
~~~

The finite selector-edge measures have density \(L_\beta(S,\cdot)\) with
respect to \(\nu_\beta\):

~~~math
\boxed{
H^2(w_Sq_S,w_Tq_T)
=\frac12\mathbb E_{\nu_\beta}
\left(\sqrt{L_\beta(S,D)}-\sqrt{L_\beta(T,D)}\right)^2.
}
\tag{R29.3}
~~~

#### 2. Verified matched zero-temperature theorem

Fix one finite signing \(A\). Let \(q=Q(A)\), let \(\mathcal G\) be its
oriented parent grounds, and put \(G=|\mathcal G|\). For each selector,
let \(Q_S=Q(A[S])\), and set \(Q_*=\max_SQ_S\). Define

~~~math
\mathcal P_*=
\{(S,y):Q_S=Q_*,\ c_S(y)=Q_*\},
\qquad N_*=|\mathcal P_*|.
~~~

For \((S,y)\in\mathcal P_*\), let

~~~math
e_S(y)=\max_{d:d[S]=y}\langle A,d\rangle,
\qquad
k_S(y)=|\{d:d[S]=y,\ \langle A,d\rangle=e_S(y)\}|.
~~~

Define a probability law on parent cuts by

~~~math
\boxed{
P_\infty(d)=\frac1{N_*}
\sum_{(S,y)\in\mathcal P_*}
\frac{\mathbf1\{d[S]=y,\ \langle A,d\rangle=e_S(y)\}}{k_S(y)}.
}
\tag{R29.4}
~~~

With \(\bar e_*=\sum_dP_\infty(d)\langle A,d\rangle\),

~~~math
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=\beta(q-\bar e_*)+\log G-H(P_\infty)+o_A(1)
\qquad(\gamma=\beta\to\infty).
}
\tag{R29.5}
~~~

**Proof.** Since
\(Z_S(\beta)=g_Se^{\beta Q_S}(1+o_A(1))\), the selector law
\(\pi_0(S)=U_m(S)Z_S/Z_0\) retains exactly the selectors with
\(Q_S=Q_*\), weighted by child-ground count \(g_S\). Conditioned on \(S\),
the child Gibbs law is uniform on those grounds. Thus every pair in
\(\mathcal P_*\) has probability \(1/N_*\). The parent conditional law
then becomes uniform on its \(k_S(y)\) best extensions, proving
\(P_\beta\to P_\infty\). Finally,

~~~math
D(P_\beta\Vert\nu_\beta)
=-H(P_\beta)-\beta\mathbb E_{P_\beta}\langle A,D\rangle
+\log Z_A(\beta),
~~~

and \(\log Z_A(\beta)=\beta q+\log G+o_A(1)\).

The slope vanishes iff every maximal child ground has a parent-ground best
extension. When it vanishes, \(P_\infty\) is supported on \(\mathcal G\)
and the constant is exactly

~~~math
D(P_\infty\Vert U_{\mathcal G})=\log G-H(P_\infty).
\tag{R29.6}
~~~

This is a sharp signing-realizable falsifier for a theorem claimed uniformly
in low temperature. It is not itself a falsifier in the intended regime
\(\beta=\Theta(n^{-1/2+c})\).

The finite Hellinger limit is also explicit. On an active selector with
\(g_S\) child grounds,

~~~math
w_S^\infty=\frac{|\{S'\}|g_S}{N_*},
~~~

and \(q_S^\infty\) chooses a child ground uniformly, then one of its best
extensions uniformly. Inactive selectors have \(w_S^\infty=0\). Substitution
in (R29.3) gives the audited limit.

#### 3. Exact finite audit

All five exact minimizers pass the zero-slope condition.

| signing | \(q\) | \(Q_*\) | \(N_*\) | parent grounds | slope | limiting parent entropy | limiting mean finite \(H^2\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| \(A_4\) | 8 | 6 | 4 | 2 | 0 | 0 | \(2/3\) |
| \(A_5\) | 8 | 8 | 10 | 10 | 0 | 0 | \(3/4\) |
| \(A_6\) | 10 | 8 | 60 | 12 | 0 | 0 | \(1/5\) |
| \(A_8\) | 20 | 18 | 24 | 8 | 0 | 0 | \(5/7\) |
| \(A_9\) | 24 | 24 | 28 | 25 | 0 | \(0.035686818862212\ldots\) | \(185/224\) |

For \(A_9\), the limiting parent masses are eight states of mass \(3/56\),
fifteen of mass \(1/28\), and two of mass \(1/56\). Hence its constant is

~~~math
8\frac3{56}\log\frac{75}{56}
+15\frac1{28}\log\frac{25}{28}
+2\frac1{56}\log\frac{25}{56}
=0.035686818862212\ldots.
\tag{R29.7}
~~~

The checker compares these limits with direct enumeration at
\(\beta=\gamma=8\). Thus the linear divergence in (10.812) is genuinely
caused by \(\gamma=0\): it disappears on all audited examples when the
temperatures match. The nonzero \(A_9\) constant shows that ground
extendibility alone does not make the parent barycenter uniform.

#### 4. Verified high-temperature common-mode identity

Let \(k=n-m\), and keep \(A,n,m\) fixed while \(\gamma=\beta\to0\).
Represent \(d=\sigma xx^{\mathsf T}\), and define

~~~math
B_S(d)=\sum_{u\notin S}
\left(\sum_{i\in S}a_{ui}x_i\right)^2,
\qquad
g_S(d)=mk-B_S(d).
\tag{R29.8}
~~~

Under a uniform outside completion, the external energy has mean zero and
variance \(4[B_S(d)+\binom{k}{2}]\). Also

~~~math
\log C_\beta=k\log2+
2\beta^2\left[mk+\binom{k}{2}\right]+O_A(\beta^3).
~~~

Therefore

~~~math
\boxed{
L_\beta(S,d)=1+2\beta^2g_S(d)+O_A(\beta^3).
}
\tag{R29.9}
~~~

Put

~~~math
\theta_{n,m}
=\frac{(n-m)m(m-1)}{n(n-1)(n-2)}.
~~~

Expanding \(B_S\) and averaging the selector gives

~~~math
\boxed{
\mathbb E_{S\sim U_m}g_S(d)
=\theta_{n,m}\left[n(n-1)-R_2(d)\right].
}
\tag{R29.10}
~~~

Consequently

~~~math
\boxed{
f_\beta(d)
=1+2\beta^2\theta_{n,m}
\left[n(n-1)-R_2(d)\right]+O_A(\beta^3),
}
\tag{R29.11}
~~~

and, since
\(\operatorname{Var}_{U}(R_2)=4\sum_{i<j}(A^2)_{ij}^2\),

~~~math
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=8\beta^4\theta_{n,m}^2
\sum_{i<j}(A^2)_{ij}^2+O_A(\beta^5).
}
\tag{R29.12}
~~~

Every adjacent selector edge similarly has

~~~math
\boxed{
H^2(w_Sq_S,w_Tq_T)
=\frac{\beta^4}{2}
\mathbb E_U(g_S-g_T)^2+O_A(\beta^5).
}
\tag{R29.13}
~~~

Thus selector Hellinger sees variation of \(g_S\) around its selector mean,
whereas parent entropy sees the common selector mean, exactly a centered
row-square statistic. A selector-edge theorem cannot control this common
mode without a minimizer-specific input.

The \(O_A(\cdot)\) remainder is fixed-\(A\), not uniform in \(n\). Since
the target has \(\beta=\Theta(n^{-1/2+c})\), this is a structural tangent,
not an asymptotic proof.

The exact \(\beta^4\) coefficients
\((\text{parent entropy},\text{mean finite }H^2)\) at one deletion are:

~~~text
A4: (4, 40/3)
A5: (16/5, 28)
A6: (0, 48)
A8: (6, 648/7)
A9: (352/27, 356/3).
~~~

#### 5. A matched signing-realizable generic obstruction

Minimality is essential. Take the actual complete signing \(A=J-I\), with
matched temperatures and \(m=n-1\). Then \((A^2)_{ij}=n-2\) for
\(i\ne j\), so (R29.12) has coefficient
\(4(n-1)(n-2)^2/n\). A direct Rademacher calculation in (R29.13) gives
coefficient \(4(n-2)\). Therefore

~~~math
\frac{\operatorname{Ent}_{\nu_\beta}(f_\beta)}
{\mathbb E_{S\sim S'}H^2(w_Sq_S,w_{S'}q_{S'})}
\longrightarrow
\frac{(n-1)(n-2)}n
\qquad(\beta\downarrow0).
\tag{R29.14}
~~~

No \(o(n)\)-coefficient absorption can hold for all actual matched
complete-signing endpoints. This strengthens the abstract two-bit warning
(10.813), but does not touch a minimizer-specific fixed-density theorem.

#### 6. Consequences for (10.800)

- **Verified advance:** the old \(\gamma=0\) linear wall does not survive
  matching on \(A_4,A_5,A_6,A_8,A_9\). Its exact matched replacement is
  (R29.5), with a sharp extendibility criterion and residual ground-
  barycenter KL.
- **Verified obstruction:** even for a genuine complete signing at matched
  temperature, generic parent-to-selector Hellinger absorption can need an
  order-\(n\) coefficient. Exact minimality and the intended density regime
  cannot be discarded.
- **Open target:** prove both bounds in (10.800), uniformly for selected
  exact minimizers at \(\beta=\gamma=\Theta(n^{-1/2+c})\). A parent input
  must control the common reciprocal-partition mode; infinitesimally this is
  the row-square fluctuation in (R29.12).
- **Falsification:** a minimizer family directly violating either scaling in
  (10.800) falsifies the package. For a theorem claimed uniformly as
  \(\beta\to\infty\), positive \(q-\bar e_*\) in (R29.5) is already an
  exact obstruction.

No convergence proof or asymptotic minimizer counterexample is obtained.
