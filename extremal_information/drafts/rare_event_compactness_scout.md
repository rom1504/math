# Orthogonal scout: rare-event / extremal compactness

Status: literature-grounded scout report.  Portfolio judgment: **keep warm,
do not promote to the main theory yet**.

## Candidate object

For a landscape on `Sigma_N={-1,1}^N` with compact context mark
`q_N:Sigma_N->K`, form the empirical action--energy measure

```math
\mu_N=2^{-N}\sum_\sigma
 \delta_(q_N(\sigma),H_N(\sigma)/N).                           \tag{RE.1}
```

Weak convergence forgets an individual maximizer because its mass is
`2^{-N}`.  The proposed leading-scale object is instead the speed-`N` large
deviation action frontier: a rate `J(q,e)` obtained through the
large-deviation/Gamma topology developed by Mariani.  Under an LDP and the
usual exponential-tightness/no-escape hypotheses, a Laplace principle gives

```math
{1\over N}\max_\sigma\{H_N(\sigma)+Ng(q_N(\sigma))\}
\longrightarrow
\sup_{J(q,e)<\infty}\{e+g(q)\}.                               \tag{RE.2}
```

This is not finite-state response compression.  It is an asymptotic topology
whose finite-action support remembers exponentially rare states.

At order-one resolution, even `J` is insufficient.  One must additionally
retain a marked extremal point process on each exposed face.  The resulting
two-scale candidate is

```math
(\text{action frontier},\ \text{marked boundary point process}).
```

The first component controls extensive maxima; the second controls ties,
clusters, and fluctuations.  This is deliberately orthogonal to a reusable
finite composition quotient.

## Benchmark imported from primary literature

For the random energy model, Bovier--Klimovsky's extremal-process theorem
(their GREM result specialized to zero field and one level) gives, with
`beta_c=sqrt(2 log 2)` and the standard logarithmically corrected centering
`c_N`,

```math
\sum_\sigma\delta_{\beta_c(H_N(\sigma)-c_N)}
\Longrightarrow PPP(e^{-x}\,dx).                              \tag{RE.3}
```

Thus the leading action frontier locates the ground-state density while the
point process supplies the Gumbel fluctuation and multiplicity law.  A paired
REM, obtained by duplicating every random level and its context mark, has the
same leading frontier and normalized maximum but a size-two cluster process.
This is an exact falsifier of “the leading response frontier contains all
extremal information.”

Primary sources:

- [Mariani, *A Gamma-convergence approach to large deviations*](https://arxiv.org/abs/1204.0640)
- [Bovier--Klimovsky, *Fluctuations of the partition function in the GREM
  with external field*](https://arxiv.org/abs/0805.1478)

## Decisive boundary

There is no universal full-sequence compactification for arbitrary landscape
sequences.  Put one height-`N` state at all-plus for even `N` and all-minus
for odd `N`, all other heights zero, and mark by magnetization.  For fixed
`0<h<1/2`, the normalized field response alternates between `1+h` and `1-h`.
Neither a joint frontier nor its marked extremal process has a full-sequence
limit.  At most subsequential compactness or a metastate is automatic.

This boundary is consonant with the fixed-model phenomenon proved by
[Chazottes--Hochman](https://arxiv.org/abs/0907.0081): zero-temperature Gibbs
states need not converge even for a fixed Lipschitz shift potential, and in
dimension at least three the phenomenon can occur for finite-range
interactions.

## Promotion judgment

Keep the direction warm because it solves a benchmark that leading contextual
response deliberately cannot: order-one extremal clustering.  Do not promote
it during this campaign because:

1. it presently gives only subsequential compactness for arbitrary
   deterministic landscapes;
2. it has no all-order realization or composition theorem;
3. the motivating normalized signing objective asks first for a leading-scale
   limit, where the rate frontier rather than the point-process refinement
   would have to acquire additional structure.

Promotion requires one theorem showing that a constrained deterministic
class has a unique all-order action frontier closed under its natural
composition, or an imported recovery theorem that contextual-response states
cannot express.
