# Discrepancy critique of terminal coset drift

## Verdict

The foreign mechanism is an **arbitrary-root coset graph** with an approximately radial intersection number: the outward degree of every root should asymptotically depend only on its distance layer. Discrepancy normally bounds a supremum or finds one good coloring; it does not analyze coordinate-flip geometry for every coloring modulo a code.

This is genuinely distinct from adversarial-pressure interpolation, but it is not yet a new proof route. The coding tools do not prove the local law, and the law hides a uniform margin theorem for `2^n` Walsh tests. The Paley audit below is more serious: `L_drift` would force the limit to be exactly `1/2`, so it contains the full sharp lower-constant problem. My truth confidence falls from `0.08` to below `0.01`.

## Foreign coding mechanism and its boundary

The relevant import is a rooted Krawtchouk/Terwilliger certificate with coverage localizers at **every ambient root**. Delsarte's annihilator polynomial has the right all-coset quantifier, and Gijswijt--Polak show how all-root covering inequalities enter Terwilliger/Lasserre localizers. This is unlike MacWilliams data based at zero or Schrijver's codeword-root packing SDP.

It is not an available black box:

* Delsarte external distance is all-coset but coarse.
* With dual distance `4`, Tietavainen yields only deficit `Theta(sqrt(N_n))=Theta(n)`, below the target `Theta(n^(3/2))`.
* Bazzi's limited-independence obstruction warns that low or slowly growing degree cannot generally force the needed all-root improvement.
* Existing all-root SDP bounds optimize unrestricted covering-code size, not the prescribed augmented cut code or neighboring-coset degrees.

A successful import must therefore be cut-code-specific, growing-degree, and all-root. That would be new, but it is almost the whole missing proof.

## Exact discrepancy content of `b_n(U)`

Represent a coset by an edge-sign vector `a`, with
\[
 Q(a)=\max_{c\in\mathcal C_n^+}\langle a,c\rangle=N_n-2r_n(U).
\]
For an edge `e`, set
\[
 q_e^\pm(a)=\max\{\langle a,c\rangle:c\in\mathcal C_n^+,\ a_ec_e=\pm1\}.
\]
Flipping `e` gives
\[
 Q(a^{(e)})=\max\{q_e^+(a)-2,q_e^-(a)+2\},
\]
and hence
\[
 \boxed{e\in b_n(U)\iff q_e^+(a)=Q(a)\ \text{and}\ q_e^-(a)\le Q(a)-4.}
\]
Thus `b_n(U)` counts coordinates on which every active maximizer has the favorable sign and every oppositely signed near-maximizer has a four-unit margin. `L_drift` says this active-face statistic, uniformly for **every** root in the terminal window, is asymptotically determined by the scalar `Q(a)/n^(3/2)`. Its information state is small, but its geometric assertion is strong.

## Exact `n<=8` experiment

| `n` | `rho_n` | `M_n` | Departure from radial drift |
|---:|---:|---:|---|
| 3 | 0 | 3 | one coset |
| 4 | 1 | 4 | one `b` per layer; no nondeep dead end |
| 5 | 3 | 4 | one `b` per layer; no nondeep dead end |
| 6 | 5 | 5 | `Q=11`: `b in {10,13}`; `Q=9`: `b in {0,4}`, with 25 nondeep dead ends |
| 7 | 6 | 9 | `Q=13`: `b in {6,9,12,14,17}`; `Q=11`: `b in {0,2,4,5,10}`, with 1,260 dead ends |
| 8 | 9 | 10 | `Q=16`: `b=0..22` (35 dead); `Q=14`: `b=0..21` (10,920); `Q=12`: `b in {0,1,3}` (119,700) |

At `n=8`, the `Q=12` layer has 168,420 cosets, so over 71% are nondeep dead ends. Its normalized value is `12/8^(3/2)=0.530330...`, only `0.02033...` above `[0.33,0.51]`; the deepest value is `10/8^(3/2)=0.441942...`. This does not formally refute an asymptotic lemma on the smaller interval, but it refutes a finite near-equitability heuristic.

## Square-field Paley audit: a decisive hidden consequence

Let `q=s^2`, where `s` is an odd prime power, and let `C_q` be the symmetric Paley conference matrix of order `n=q+1`. It satisfies `C_q^2=qI`. Partition `F_q` into the `s` additive cosets of a square one-dimensional `F_s`-subspace. Within a block the Paley sign is `+1`, while the character sum from a vertex to each other block is `-1`.

Choose block signs `t_1,...,t_s in {+1,-1}` with `T=sum t_i=+1` or `-1`; put the infinity coordinate equal to `T` and make the field coordinates constant `t_i` on block `i`. The block character sums give
\[
 C_qx=sx.
\]
Hence these Boolean vectors attain the spectral bound
\[
 Q(C_q)=\frac{ns}{2}=\frac n2\sqrt{n-1},\qquad
 z(C_q)=\frac12\sqrt{1-1/n}\longrightarrow\frac12.
\]

Moreover every edge can be opposed by one such active eigenvector. For an infinity--field edge, choose its block sign opposite `T`. For a field edge, choose a square-line partition in which the endpoints lie in different blocks and prescribe their block-sign product; the remaining odd number of block signs can still make `T=+/-1`. Flipping that edge raises the energy of this active vector by `2`. Since one flip changes every energy by at most `2`, every flip has new discrepancy exactly `Q(C_q)+2`. Therefore
\[
 \boxed{b_n(C_q)=0\quad\text{for infinitely many }n=q+1,\qquad z_n(C_q)\to1/2.}
\]
The order-10 case (`q=9`) checks exactly: `Q=15`, there are 24 Boolean eigenvectors, and all 45 edge flips give `Q=17`; meanwhile `M_10=13`, so this is a nondeep dead end.

Because every admissible interval contains `1/2`, equation (10), continuity, and these Paley roots force `beta(1/2)=0`. The unique-zero clause then forces `c=1/2`. Applying the same clause to deepest roots yields
\[
 \lim_{n\to\infty}M_n/n^{3/2}=1/2.
\]
Thus `L_drift` implies the full sharp lower bound `liminf >=1/2`, not merely unspecified convergence. It remains strictly weaker than full finite coset histograms as information, but it is **not** a strict reduction in asymptotic difficulty: proving it must prove the sharp constant and exclude every other asymptotic dead-end level.

## Strongest failure mechanism and discrepancy rescue

The quotient distance is not geodesically convex. Different exposed faces of `Q(a)=||W_na||_infinity` have incompatible active-row signs, producing coordinatewise local minima away from the global minimum. Scalar radius does not determine the constrained second maxima `q_e^-`. Paley supplies an infinite nonoptimal zero-drift family at `z->1/2`; any other dead-end family converging inside the interval to `z!=1/2` immediately falsifies `L_drift`.

The closest discrepancy rescue is Lovett--Meka: for one fixed root, treat active/near-active codewords as rows and seek a multi-edge perturbation with nonuniform thresholds. Its entropy budget can charge only dangerous rows rather than union-bound `2^n` tests. But it rescues at most a **block-drift** variant: it does not create improving single-coordinate flips, assumes the missing uniform near-maximal tail bound, is existential for one root rather than all roots, and its constants/recursion do not preserve a leading coefficient. Bednorz--Latala can price Bernoulli-supremum geometry, but remains expected, multi-coordinate, and nonuniform under deterministic root selection.

## Broad cross-domain check

The action spectral-regularization lemma (SR) asks for `o(n^(3/2))`-near-minimizers with `||A||op=O(sqrt n)`. Discrepancy sees a simultaneous matrix-signing problem: operator directions and all Walsh rows must be controlled by the same edge signs. Random/noncommutative bounds control the operator norm and Spencer-type bounds control Walsh discrepancy separately, but no theorem combines them at zero leading cost. This adds no new no-go beyond the recorded fixed-`K` tradeoff, and no rescue of action recovery. The pressure proposal already pays for all tests by log-sum-exp; discrepancy adds only the known adversarial-minimum obstruction. The Pythagorean composition lemma has a clearer vector-balancing floor: its cross block `D` must satisfy `2^(n+m)` offset constraints jointly, while generic discrepancy for `nm` variables is `Theta(sqrt(nm(n+m)))`, exactly the leading `(n+m)^(3/2)` scale. Absolute control of `x^TDy` therefore pays a fixed leading term; only profile-aware cancellation with `H_A(x)+H_B(y)` could prove the claimed sublinear squared error, and no standard theorem does so.

## Fastest settling test

Track
\[
 S_n(Q)=\frac{\max_{U:Q(U)=Q}b_n(U)-\min_{U:Q(U)=Q}b_n(U)}{N_n}
\]
and search for nondeep `b=0` roots. Exact `n=9` should inspect `Q=14` (`z=0.518518...`); targeted SAT/MILP at `n=10` should recover the Paley `Q=15` trap inside the interval. The decisive asymptotic test is an explicit lift of any second dead-end family to a limit `z!=1/2`, or proof that all such families collapse to `1/2` and same-layer scatter vanishes.

**Bottom line:** arbitrary-root localizers are genuinely foreign, but Paley roots show that the unique-zero drift conjecture secretly asserts the sharp constant `1/2`. Preserve the disagreement: this is a low-confidence falsification target, not a weaker substitute for pressure convergence.
