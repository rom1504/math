# Independent audit: same-frame finite-cell response ascent

Date: 2026-09-06. Status: the structural objective and its first variation
pass independent reconstruction. The floating optimizer is a policy
discovery tool, not a numerical certificate. Its resulting rectangle
policy can be checked by the previously audited exact birth certificate.

Audited source:
`computations/resumed_response_rich_core_cell_ascent_2026_09_06.py`.
The underlying actual rich-core frame and conditional projection are
audited in `resumed_bound_audit_rich_core_response_birth_2026_09_06.md`.

## 1. Fixed actual frame and feasible cell policies

Let V,Z be the fixed independent standard Gaussian coordinates from the
rich-core birth construction. Their inverse features are g and
h_Z=(f-A g)/nu, which are orthonormal in L2. No new feedback equation
or Gaussian coordinate is introduced during the optimization.

Positive V bins and Z bins form rectangles C. Their reflections under
(V,Z)->(-V,-Z) receive the odd/even reflected policy. On each paired
cell choose constants f_C,q_C with |f_C|<=q_C<=1. Independent even
Gaussian gates realize such a conditional feasible pair; in particular
they preserve first-chaos coefficients and realize support mass q_C.
One must use this realization, not the squared norm of the pointwise
convex interpolation of response functions.

The fixed tail policy is sign(V) on |V|>2, and sign(Z) on
|V|<=2, |Z|>8. Its exact support and first-chaos contributions are

    p_tail=2 Phi(-2)+(1-2 Phi(-2))*2 Phi(-8),
    aV_tail=2 phi(2),
    aZ_tail=(1-2 Phi(-2))*2 phi(8).

These agree with the source. Write w_C for paired-cell probabilities and
v_C,z_C for coordinate means. Then

    p=p_tail+sum_C w_C q_C,
    aV=aV_tail+sum_C w_C f_C v_C,
    aZ=aZ_tail+sum_C w_C f_C z_C.

The new inverse feature is K1=aV*g+aZ*h_Z and its TRUE squared norm is
aV^2+aZ^2, not a conditional-cell approximation to that norm.

## 2. Exact Gram and Jensen objective

Let K0,p0,H0 denote the previously feasible baseline. Mix baseline and
new response using an independent even gate of probability theta. Then

    K=(1-theta)K0+theta K1,
    t=(1-theta)p0+theta p-||K||_2^2,
    H_C=(1-theta)H0_C+theta(1-q_C).

The baseline mask is constant on every V bin because the alpha boundary
is explicitly included. The source evaluates ||K||^2 from the exact
old norm and true old/new inner products:

    <K0,g>=lambda+gamma*c,
    <K0,h_Z>=gamma*(<u,f>-A*c)/nu.

With exact conditional means K_C=E[K|C], the objective

    J_cell=sum_C w_C H_C Psi(K_C,t)

is a rigorous lower bound on the full response value by convexity of
Psi(k,t)=E|k+sqrt(t)N|. The omitted tails have nonnegative objective
contribution. Their full support and true Gram contributions remain in
t, as required. No independence between the center mask and g is used.

The numerical source substitutes a finite Z-Hermite approximation to
g_C. Consequently its printed value is not itself rigorous. The
previously proved principal-angle tail bound supplies the needed
uniform L2 correction when freezing a policy.

## 3. First variation and feasible ascent oracle

For fixed theta, put s_C=2 Phi(K_C/sqrt(t))-1 and

    B=sum_C w_C H_C phi(K_C/sqrt(t))/sqrt(t),
    rV=sum_C w_C H_C s_C g_C-2B<K,g>,
    rZ=sum_C w_C H_C s_C hZ_C-2B<K,h_Z>.

Direct differentiation, including support dependence of t, gives

    dJ=theta sum_C w_C [
       (rV*v_C+rZ*z_C) df_C+(B-Psi(K_C,t)) dq_C].

Thus the exact cellwise linear oracle is

    q_C*=1{|rV*v_C+rZ*z_C|>Psi(K_C,t)-B},
    f_C*=sign(rV*v_C+rZ*z_C) q_C*.

The source matches every term and factor. Convex interpolation of cell
pairs stays feasible, and scalar line search is legitimate policy
ascent. Neither the small final oracle gap nor the numerical optimizer
proves global optimality.

## 4. Concrete output and inexpensive certificate handoff

The saved floating run reaches .4333545239115277 in a 256-by-890 grid.
Independent inspection of its last pure oracle gives 256 monotone Z
rows: 28 have one transition and 228 have two; none has more than two.
It therefore compresses to one lower/upper center interval per V bin.

The coarser weighted-V Jensen value from this pure policy is
.43332537570. Freezing that policy at theta=1 permits reuse of the
existing exact 256-rectangle certificate, with the same true Gram and
conditional-Hermite tail control.

## 5. Independent exact replay: improved universal bound

The independent replay completed successfully using

    .venv/bin/python computations/resumed_response_rich_core_birth_certificate_2026_09_06.py
      --policy computations/results/resumed_response_rich_core_optimized_rectangle_policy_2026_09_06.json
      --theta 1 --target 4333/10000
      --output /home/math/quadra/tmp/resumed_bound_audit_optimized_rich_birth_replay_2026_09_06.json

The resulting JSON is BYTE-IDENTICAL to the response agent's canonical
`computations/results/resumed_response_rich_core_optimized_policy_certificate_2026_09_06.json`.
The only certificate-source parameterization change from the previously
audited mixed-policy certificate is the optional theta argument.

Exact interval endpoints include

    p_new >= .465968240263521924604959372492350930979274860709375133179917,
    t_new >= .092031354063631380452808124083501002341145243864005647776561,
    pre-tail bin sum >= .433325375702943007494944152418094290701027670940273741230068,
    tail charge <= .000003264038862254079131223520514944142393022906580328123072.

After that charge, the independently reproduced universal lower bound is

    liminf_n M_n/n^(3/2)
      >= .433322111664080753415812928897579346558634648033693413106996.

For a short downward-rounded statement, use .4333221116640807.
This claim depends on the previously reconstructed bounded-operator
full-response theorem, actual rich-core realization, exact frozen
rectangle policy, and principal-angle tail bound; the floating ascent
history is not a proof dependency.
