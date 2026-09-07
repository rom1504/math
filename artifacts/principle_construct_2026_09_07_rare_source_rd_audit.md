# Independent rare-source rate-distortion audit

2026-09-07. **PASS** for the analytic lemma sent by the invention agent.
This audit does not yet certify any claimed actual cap application.

Let the fixed distinct coarse atoms be 0,a_1,...,a_d, with a_j nonzero,
and P(Z!=0)<=Kr. Suppose |X-Z|<=Kr almost surely. For fixed tau>0 set
t=tau log(1/r). Then, uniformly over all allowed atom probabilities,

    E_t(X)=-log(1/r) sum_(a!=0) P(Z=a) min(tau a^2,1)
                  +o(r log(1/r)).

## Uniform lower bound for the rate-distortion cost

For sufficiently small r, the disjoint neighborhoods of the fixed atoms
make Z a measurable function of X. Thus every channel L satisfies
I(X;L)>=I(Z;L)>=I(1_(Z!=0);L). Write p=P(Z!=0),
q(L)=P(Z!=0|L), and delta=(log(1/r))^(-1/2). The p=0 case is direct.

On labels with q>delta, binary relative entropy is at least
`q log(1/p)-h(q)`. Since
`h(q)<=q[log(1/delta)+1]` there, and log(1/p)>=log(1/r)-log K, the
information cost is at least rare-high-label mass times log(1/r),
minus O_K(r log log(1/r)). The low-label relative entropy is nonnegative
and can be dropped. This argument remains uniform if p is much smaller
than r; it does not replace log(1/p) by a smaller expression incorrectly.

On low labels the conditional mean of Z has magnitude at most Amax q.
Their total conditional variance is therefore at least

    sum_a a^2 P(Z=a,low)-Amax^2 delta p.

Conditional variance under |X-Z|<=Kr obeys
`Var(X|L)>=(1-delta)Var(Z|L)-delta^(-1)E[(X-Z)^2|L]`.
After averaging and multiplying by t, its additional loss is
O(t r^2/delta), and all preceding variance errors are o(r log(1/r)).
Combining high-label information and low-label variance, each rare atom
is paid at least min(1,tau a^2) times its mass and log(1/r), up to the
uniform stated remainder. Taking the infimum over channels is legitimate
because none of these errors depends on the channel.

## Matching upper bound and passage to the envelope

Reveal only atoms with tau a^2>1 and merge all other atoms into one label.
The label is a deterministic function of X. Its entropy is
`p_exp log(1/r)+O(r)` uniformly: writing each rare probability as ru,
the correction r u log(1/u) is uniformly bounded for u in [0,K], with
only finitely many atoms. The remaining averaged conditional variance is
`sum_cheap p_a a^2+O(r^2)`. This proves the matching upper bound for J_t.

Every conditional variance is at most Var(X)=O(r). The Gaussian reward
satisfies `g_t(V)=-tV+O(t^2 r^2)` uniformly on this interval, and
t^2r^2=o(r log(1/r)). Hence E_t(X)=-J_t(X)+o(r log(1/r)), proving the
displayed equality. Critical atoms tau a^2=1 may be assigned either way.

The proof requires fixed separated coarse atoms (or an explicitly uniform
separation version). It does not cover arbitrary moving rare amplitudes
without additional error control. No actual signing, common-temperature
profile comparison, or diagonal completion is inferred from the lemma
alone.
