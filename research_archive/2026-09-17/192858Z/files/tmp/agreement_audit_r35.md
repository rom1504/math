# Independent audit: Wave 35 singleton-anchored agreement

## Verdict

The elementary decoder and its exponent are **Accepted**.  The singleton
anchor removes the multi-coordinate compatibility gap: every projective label
has a unique representative equal to `+1` at `v_*`.  The finite checker passes.

The multi-block minimax and row-qualified inference are **Accepted after two
explicit quantifier/domain clarifications**:

1. In (10.1004), state
   `xi^S in [-1,1]^{E(S)}` for every selector.  Without this cube domain the
   displayed minimum is formally unbounded below.
2. Require `kappa>0` in (10.1004)--(10.1007).  Division by `kappa` is what
   turns (10.1006) into the conflict bound.  At `kappa=0`, (10.1007) does not
   imply any bound on conflict.

Also say that (10.1006) holds for every **optimizing** dual law `mu`; it is not
an identity for an arbitrary supported law.

## Decoder rederivation

For an outside coordinate `i`, put

```math
u_i=Pr(i in S,y_i^S=+1),\quad
v_i=Pr(i in S,y_i^S=-1),\quad p_i=u_i+v_i.
```

Its plurality-error contribution is `min(u_i,v_i)`.  Its contribution to
the independently sampled, degree-corrected conflict is
`2u_i v_i/p_i`.  The latter dominates the former because the larger of
`u_i,v_i` is at least `p_i/2`.  If `p_i=0`, both contributions vanish.  A tie
gives equality.  The anchor contribution vanishes because every representative
equals `+1` there.  This proves (10.1002).

Integer Markov gives a fraction at least
`(1-C/(d+1))_+` of the selected family at local error at most `d`.  Since the
family has density `beta` within the singleton-anchored slice, whose uniform
slice mass is exactly `m/n`, (10.1003) follows.  Ordinary Hamming error to a
selected favorable representative upper-bounds projective distance.  Under
`C<=(d+1)/2`, the resulting mass is at least `beta(m/n)/2`, hence
`exp{-O(TL_0)}` when `m/n` is bounded below and
`log beta^{-1}=O(TL_0)`.

## Minimax and accounting

After mixing over entire assignments, finite minimax and independent
minimization of each fractional block give

```math
min_{xi^S in [-1,1]^{E(S)}}2 sum_e xi_e^S m_e^S
=-2 sum_e |m_e^S|,
```

which is (10.1005).  For an optimizing law,

```math
E Delta=q_n-E Z-2 E_S sum_e a_e m_e^S.
```

Subtracting the dual value uses
`2(|m|+a m)=4(a m)_+` for `a in {+-1}`, proving (10.1006), including the
factor four and all signs.

For `kappa>0`, (10.1007) and nonnegativity imply
`E_mu C_{v_*}<=(d+1)/2`.  Some assignment in the support has conflict no
larger than this mean.  Every support assignment is row-qualified, and all
coordinatewise plurality tie choices are retained in `R_med`; choose its
low-row median and apply (10.1003).  This proves (10.967).

With `kappa=4(n-1)` and
`d=O(T n^{1/2-2c_0}/(log n)^2)`, the permitted loss is exactly

```math
2(n-1)(d+1)
=O(T n^{3/2-2c_0}/(log n)^2+n),
```

so (10.1008) is correct.

## Literature wording

The revised Dikstein--Dinur wording is accurate: `1/log m` is explicitly the
introductory regime, while the formal theorem has additional suitability and
base-test-soundness assumptions.  Its formal accuracy parameter is not simply
fixed at `0.01`.

For Dinur--Steurer, a maximally precise phrase is that its displayed guarantee
only gives an `O(TL_0)` bound on the expected number of local coordinate
errors under the favorable parameter substitution; this stated guarantee is
too coarse to force `d=o(TL_0)`.  Saying it "has error of order `TL_0`" could
be misread as a lower bound, but this is only a wording issue.  The decisive
list-input and parameter mismatches are correctly stated.

