"""Exact three-inequality obstruction to a uniform additive RD majorant."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    entropy_interval,log_interval,sqrt_interval)


def certificate():
    p=F(31,32);t=F(4);s=F(9,10);cut=F(9,10)
    v=1-p*s*s
    entropy_lower=entropy_interval((1+s)/2)[0]
    sqrtp_upper=sqrt_interval(p)[1]
    radical_low,radical_high=sqrt_interval(1+16*t*t)
    rho_low=(radical_low-1)/(4*t)
    rho_high=(radical_high-1)/(4*t)
    gaussian_lower=-t*(1-rho_low)+log_interval(1-rho_high*rho_high)[0]/4
    base=p*entropy_lower+t*(1-sqrtp_upper)
    low_temperature=base-cut*v
    high_left=base+gaussian_lower+F(1,2)+log_interval(2*cut)[0]/2-cut*v
    high_right=base+gaussian_lower+F(1,2)+log_interval(2*t)[0]/2-t*v
    answer=min(low_temperature,high_left,high_right)
    assert answer>0
    values={'lambda_at_most_9over10':low_temperature,
        'concave_interval_left_endpoint':high_left,
        'concave_interval_right_endpoint':high_right,
        'gaussian_potential_lower':gaussian_lower,
        'uniform_forced_exponent_lower':answer}
    return {'p':str(p),'t':str(t),'channel_correlation':str(s),
        'conditional_variance':str(v),
        'bounds':{key:{'rational':str(value),'display':float(value)}
                  for key,value in values.items()},
        'strictly_positive':True,
        'scope':'uniform additive fixed-temperature terminal majorants only; not a Bellman lower policy or exclusion of the candidate'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=certificate()
    print(json.dumps(result,indent=2))
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
