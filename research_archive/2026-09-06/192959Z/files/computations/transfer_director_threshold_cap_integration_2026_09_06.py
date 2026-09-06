"""Exact interval integration of the audited discrete-threshold credit.

The annealed certificate a_* is an already verified input, not recomputed here.
All comparisons use rational arithmetic; displayed decimals are diagnostic.
"""
from fractions import Fraction as F
import json
from pathlib import Path

import continued_feedback_ternary_latent_interval_certificate_2026_09_06 as iv
from transfer_seed_threshold_entropy_resampling_verify_2026_09_06 import rational_certificate


def main():
    rational_certificate()
    iv.SCALE = 10 ** 40
    p, delta = F(31, 32), F(1, 2 ** 24)
    q, v = 1 - 2 * delta, 4 * delta * (1 - delta)
    a_star = F(91470529542342299, 20460000000000000000)
    sqrtlo, sqrthi = iv.sqrt_interval(p)
    hlo, hhi = iv.entropy_interval(delta)
    base_lo = (F(1, 2) - (a_star + p * hhi) / (8 * sqrtlo)) / q ** 2
    base_hi = (F(1, 2) - (a_star + p * hlo) / (8 * sqrthi)) / q ** 2
    credit_new = F(18725, 32768) * v ** 2
    credit_clip = F(1, 2) * v ** 2 / p  # beta^2/128, beta=8/sqrt(p).
    assert credit_new > credit_clip
    assert credit_new - credit_clip > F(1301, 32768) * v ** 2
    newlo = base_lo - sqrthi * credit_new / (8 * q ** 2)
    newhi = base_hi - sqrtlo * credit_new / (8 * q ** 2)
    cliplo = base_lo - sqrthi * credit_clip / (8 * q ** 2)
    cliphi = base_hi - sqrtlo * credit_clip / (8 * q ** 2)
    assert newhi < cliplo < cliphi < F(499432211, 10 ** 9)
    result = {
        'status': 'EXACT RATIONAL INTERVAL CERTIFICATE; decimals diagnostic only',
        'annealed_input': str(a_star),
        'new_credit_rational': str(credit_new),
        'clipped_credit_rational': str(credit_clip),
        'new_cap_interval': [str(newlo), str(newhi)],
        'clipped_cap_interval': [str(cliplo), str(cliphi)],
        'strict_cap_improvement_lower': str(cliplo - newhi),
        'new_cap_decimal_diagnostic': float((newlo + newhi) / 2),
        'improvement_decimal_diagnostic': float(cliplo - newhi),
        'displayed_upper_bound_unchanged': '0.499432211',
        'convergence': 'OPEN',
    }
    output = Path(__file__).resolve().parent / 'results/transfer_director_threshold_cap_integration_2026_09_06.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if not k.endswith('_interval')}, indent=2))


if __name__ == '__main__':
    main()
