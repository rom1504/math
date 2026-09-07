# Independent audit: paid Gaussian-sign first-chaos MGF

2026-09-07. PASS for the stated bounded-correlation-operator scope.

For covariance matrices with a positive spectral lower bound, this valid estimate is now superseded by the sharper exact-covariance theorem audited in `flatify_adversary_2026_09_07_exact_covariance_sign_mgf_audit.md`. The Holder inflation below is a limitation of that earlier certificate, not of the actual sign law.

Audited `flatify_independent_2026_09_07_first_chaos_gaussian_sign_mgf.md` independently. Read the primary [Chen--Dafnis--Paouris Theorem 1(i)](https://arxiv.org/html/1306.2410v2): scalar blocks with marginal variance one and exponents all K give precisely the asserted product bound whenever R<=KI. The statement permits positive-semidefinite covariance and nonnegative measurable functions; the functions here have finite Gaussian moments. This is the only imported multivariate inequality in this audit.

The Cameron--Martin identity is exact: separating a<h,G> with a=sqrt(2/pi) shifts G by m=aRh and pays a²h^T R h/2. It does not require R to be invertible. The residual r(z)=sign(z)-az has shifted mean f(m)-am, and its shifted variance is

    1-f(m)²+a²-4a phi(m),  f(m)=2Phi(m)-1.

In particular the baseline is 1-a². Both deviations are bounded by C m² globally: Taylor expansion handles |m|<=1, while linear growth of the mean and boundedness of the variance handle |m|>=1. No small-coordinate bound on Rh has been assumed.

After centering, |r(Z+m)-E r(Z+m)|<=2+a|Z| uniformly in m. This supplies a dimension-free third-order scalar log-MGF remainder on each fixed compact source interval. Applying Gaussian Holder with source K h_i and dividing by K yields total errors bounded by constants times

    ||h||infinity ||m||²,
    K ||h||infinity² ||m||²,
    K² ||h||infinity ||h||².

Since ||m||²<=a²K²||h||², these are O_K(||h||infinity ||h||²) when the source sup norm is bounded. Thus the actual exponential theorem, not merely a covariance estimate, follows.

For the opposite-spectral construction, |A|_ii<=sqrt(n-1) by scalar Cauchy--Schwarz in the spectral measure at coordinate i. Hence K_A dominates both signs of A and has the stated constant diagonal. The sum-of-two-tensor-products formula makes R positive semidefinite with unit diagonal. Its operator norm need not be bounded for arbitrary minimizing children; that missing premise is correctly explicit.

In the flat specialization, v=x tensor y has squared norm n² and

    v^T R_rho v/n²=1-4rho[n/(n-1)]e_A e_D.

Substitution gives exactly the proposed proxy

    V=1+rho(1-2/pi)-(8rho/pi)[n/(n-1)]e_Ae_D.

The error is O(sqrt(n)) because the source has sup norm |t|/sqrt(n) and squared norm t²n. The flat spectral bound |e_A|,|e_D|<=sqrt((n-1)/n)/2 implies V>=2-4/pi>0 uniformly in rho. Therefore optimizing a bounded Chernoff source for bounded b is justified. Replacing h by -h proves the two-sided bound with the same proxy.

At rho=1, any diagonal Holder majorant P>=R has Tr P=Tr(PR)>=Tr(R²)=2n². Thus the stated average exponent obstruction is correct within this diagonal Holder method. It is not a lower bound on the actual nonlinear MGF or on every conceivable concentration argument.

The finite parent shell implication retains the essential missing obligations: all shell cardinalities and the positive spectral responses must be paid, and bounded K is not a consequence of the currently known actual-child cap bound. No original recurrence follows from this audited theorem alone.
