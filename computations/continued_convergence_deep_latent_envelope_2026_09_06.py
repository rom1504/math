"""Interval checks for the exact deep latent-copy lower envelope."""
import argparse
import json
from pathlib import Path
from mpmath import iv


def gaussian_phi(z):
    if z==0:return iv.mpf(0)
    rho=4*z/(iv.sqrt(1+16*z*z)+1)
    return -z*(1-rho)+iv.log(1-rho*rho)/4


def example(precision):
    iv.dps=precision
    t=iv.mpf(4);eps=iv.mpf(1)/100;ln2=iv.log(2)
    total_variance=1+eps*eps
    gfull=gaussian_phi(t*total_variance)
    gnoise=gaussian_phi(t*eps*eps)
    limit_lower=gnoise-ln2
    records=[]
    for r in [3,4,5,8,12]:
        n=2**r
        # X_r=sqrt(n)S+eps/sqrt(n)sum T_i. The sign S is recoverable,
        # and the sum has at most n+1 values.
        lower=(1-iv.mpf(1)/n)*(gnoise-ln2)-(ln2+iv.log(n+1))/(2*n)
        records.append({'depth':r,'Bellman_lower':str(lower),'gap_over_Gaussian_candidate':str(lower-gfull),
          'strict_gap':bool((lower-gfull).a>0)})
    assert gfull.a > (-2*ln2).b
    assert (limit_lower-gfull).a>0
    return {'source':'X=S+T/100, independent uniform signs','t':4,'precision':precision,
      'Gaussian_candidate':str(gfull),'total_entropy_candidate':str(-2*ln2),
      'deep_latent_lower':str(limit_lower),'gap':str(limit_lower-gfull),
      'finite_depth_bounds':records,
      'scope':'falsifies max{Gaussian potential, -H} as a universal deep upper envelope; source is four-atom, not original ternary'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=60)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    result=example(args.precision);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
