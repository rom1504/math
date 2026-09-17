# Primary-source reconstruction around the sign regularizer

2026-09-17. This note separates classical inputs from the campaign's
actual-sign application. It is not a claim of external novelty.

## 1. Read-k information is classical

The director read the complete eight-page original paper of Gavinsky,
Lovett, Saks and Srinivasan, *A Tail Bound for Read-k Families of
Functions*, [arXiv:1205.1478v1](https://arxiv.org/abs/1205.1478).
Its Corollary 2.8 bounds the sum of marginal KL divergences to a product
uniform reference by k times the full divergence. Apply it to the
conditional input law given U, then average over U, to obtain precisely
the read-k mutual-information inequality used in our sign regularizer.
The proof is entropy chain rule plus the coordinate multiplicity count;
independence of the reference coordinates is essential.

The director also read pp.1--7 of Lau--Nair--Ng,
[*A mutual information inequality and some applications*](https://chandra.ie.cuhk.edu.hk/pub/papers/HC/MI-Ineq.pdf),
including the full proofs of Lemmas1--4 and their fractional-information
theorem. Lemma1(iii), with a constant auxiliary conditioning variable and
fractional weights1/k, is the exact more general statement. Their
independence-conditioned supermodularity proof is valid even when the
raw variables are not uniform. Their layered-function extension requires
the stated Markov sufficiency; it cannot silently be applied to an
optimizer after discarding parts of its energy landscape.

## 2. Information-selected subGaussian bias is classical

Russo--Zou, [*Controlling Bias in Adaptive Data Analysis Using
Information Theory*](https://proceedings.mlr.press/v51/russo16.html),
Proposition1 and its supplementary proof, were read directly. The
variational-entropy argument gives selection bias at most
sqrt(2 sigma^2 I), exactly the ingredient needed here. In our application
the selected label is a spin vector, the input is one random incident
edge row, and the reference query is its absolute signed sum. Every fixed
query has the SAME mean and subGaussian bound. We use the elementary
joint-versus-product proof, reproduced in the canonical theorem, rather
than assert that the selected row stays independent of the optimizer.

## 3. What the combination adds in this problem

The [actual-sign theorem](paper_director_exchangeable_sign_regularization_2026_09_17.md)
uses exchangeability of q freshly rewritten vertices to divide optimizer
information among their rows. Each new internal edge belongs to TWO rows,
not one, and the classical read-two inequality pays that overlap exactly.
The absolute cap's exact vertex-addition identity then turns this bias
estimate into one increment that controls ALL near-level Bernoulli widths.
Sauer counting supplies the simultaneous entropy bound. Finally, exact
physical sign frames and scalar field comparison give the restricted-parent
error O(N^(4/3)(log N)^(1/3)). None of the cited information lemmas alone
supplies the full-sign construction or this application.

## 4. Nearby perturbation-stability literature: scope checked

The director read pp.1--22 of Chatterjee--Ray,
[*On the stability of solutions to random optimization problems under
small perturbations*](https://arxiv.org/abs/2410.21513), including the
statements of Theorems1.4.8,2.1.6,2.1.9 and the Metropolis--Hastings
stationary-perturbation proof used there. Its near-optimum packing
mechanism combines stationarity, Taylor control and Gaussian comparison.
The continuous-input density and rejection-control assumptions of that
particular theorem are not satisfied by exact discrete sign disorder.
It therefore does not directly replace the read-two construction. Our
separate all-energy Gaussian stability theorem concerns weighted
perturbations AFTER preparing one exact-sign instance; it does not pretend
those weighted neighbors are themselves full sign matrices.

These readings are not a reconstruction of every theorem in the long
Chatterjee--Ray paper. Preliminary searches did not establish the external
priority of the star-regularization application; no novelty claim is made.
Downloaded primary PDFs are preserved by the dated research-archive
snapshot, with originals under the campaign's director temporary directory.
