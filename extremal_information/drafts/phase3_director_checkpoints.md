# Phase 3 research-director checkpoints

This file records steering judgments, not additional theorem claims.  A
checkpoint changes the next experiment; it does not end the research run.

## Checkpoint 1: unrestricted syndrome responses

**Question.** Did the macroscopic syndrome target reveal a general law, or
only another code-specific dynamic program?

It produced two positive compression mechanisms and a stronger obstruction.

1. A word-length profile `lambda_S` is one-Lipschitz in the Hamming chart of
   any basis contained in `S`.  Storing its values on a radius-`r` covering
   code determines every adversarial appended-support response to error `r`,
   because min-plus convolution and the final maximum are sup-norm
   nonexpansive.  At `r=1` this uses

   ```math
   O(2^w\log w/w+w^2)=o(2^w)
   ```

   bits.  This answers the requested full-context quotient question much
   more strongly than a fixed-relative-error net, although composing two
   independently approximated sketches can add their errors.

2. Thresholding supplies a separate, exactly closed algebra.  Supports of
   radius below `R` form an absorbing ideal under union.  Collapse that ideal
   and retain every radius-at-least-`R` support exactly.  Kneser's theorem
   implies that such a retained support has at most
   `2|F_2^w|/(R+1)` elements.  Hence the quotient uses

   ```math
   O(2^w\log R/R)
   ```

   bits, composes with no accumulated error, and has radius error at most
   `R/2-1`.  Taking `R -> infinity` with `R=o(w)` gives the stronger target:
   an `exp(o(2^w))`-state exact composition algebra with `o(w)` response
   error.  The algebraic operation is the classical Rees quotient; the new
   mathematical content is the high-radius sparsity theorem and its response-
   compression consequence.

3. These upper bounds do not make the response class low-dimensional.  A
   Grassmannian family of dense-carrier supports has pairwise future-response
   distance controlled by subspace injection distance.  Constant-dimension
   packing yields an independently audited lower bound

   ```math
   ((1/2-2\epsilon)^2-o(1))w^2
   ```

   bits at error `epsilon*w` for every fixed `epsilon<1/4` (with the displayed
   coefficient understood as a supremum over a strict packing margin).  This
   is a moving global-carrier obstruction, not the earlier fixed direct-sum
   source.  It does not approach the exponential-in-`w` upper description.

The selected dichotomy is therefore resolved positively, while its optimal
rate remains open.  The useful abstraction is not “store the dynamic-program
table.”  Two different laws appeared:

- **contractive profile sampling:** regularity on a low-covering-entropy
  interface plus nonexpansive continuation;
- **hard-core ideal compression:** collapse an absorbing easy ideal and count
  the surviving difficult objects.

Both mechanisms are candidates for a general future-response law.  Neither
by itself explains when local information becomes macroscopically exposed.

### Immediate next theorem

The next phase will test generality rather than optimize the syndrome rate.
It must do both of the following:

1. formulate the hard-core quotient as a response theorem and validate it on
   a non-code model with a strict quotient; and
2. determine whether selective neutralization admits a non-tautological
   packing theorem that explains the block and Grassmannian lower bounds and
   applies to a second model.

If these reduce only to the definitions of a Rees quotient and a metric
packing, the framework has learned a code theorem but not a general law; the
next checkpoint must pivot to query-mass-sensitive posterior geometry.
