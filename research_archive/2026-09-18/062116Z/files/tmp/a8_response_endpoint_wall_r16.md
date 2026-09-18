# A8 response with zero endpoint service

Status: exact and independently checked by `tmp/audit_a8_response_r16.py`;
the endpoint statement is also checked by `tmp/check_endpoint_crossgram_r14.py`.

For the balanced order-eight minimizer (10.445), use the partition

```math
\{0,3,4,5\}\sqcup\{1,2,6,7\}.
```

Both internal blocks have norm 12, while `q_4=8`, so

```math
S=24,\qquad B=16,\qquad X=8.
```

The cross-only matrix has norm 16.  Exhausting all labelled order-four
minimizers in both blocks gives

```math
\min_{G_1,G_2}Q(C\oplus_{\rm diag}(G_1,G_2))=20=q_8;
```

exactly 16 labelled pairs attain it.  Hence the original and best-hybrid
common-mosaic responses are

```math
\Phi_D=20-24=-4,\qquad
\Phi_G=20-16=4,
```

and the compulsory response increase is exactly `X=8`, with zero hybrid
excess.

Nevertheless every positive/negative endpoint pair of this same `A_8` is
neutral: (10.555)--(10.559) give zero tolled endpoint residual on every shore
and every endpoint tree.  Thus no unconditional inequality can lower-bound
the established endpoint-residual service by a positive multiple of the
common-mosaic response increase.  Any response-to-temporal bridge must use a
new resource, a grouped temporal/terminal term, or an asymptotic structural
exclusion; scalar response alone is insufficient already at order eight.
