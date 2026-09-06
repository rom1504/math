"""Bounded temperature-segregation test; all numerical values are diagnostics."""
import argparse
import json
import math
from pathlib import Path

from continued_convergence_latent_supersolution_test_2026_09_06 import SourceEnvelope,g
from continued_feedback_ternary_concave_envelope_2026_09_06 import diagnostic,entropy


def ternary_value(retention,temperature):
    if temperature/retention<=.5:return g(temperature)
    result=diagnostic(retention,temperature,80)
    return result['maximum_diagnostic']-retention*math.log(2)-temperature*(1-math.sqrt(retention))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();records=[]
    for weight in (.1,.25,.5,.75,.9):
        for small in (.1,.4,.8):
            for temperature in (.25,1.,4.):
                try:
                    source=SourceEnvelope([small,1.],[weight,1-weight],temperature,121).upper(.003,120)
                    plus=ternary_value(weight,temperature*2*weight*small*small)
                    minus=ternary_value(1-weight,temperature*2*(1-weight))
                    policy=.5*(plus+minus-entropy(weight)-math.log(2))
                    record={'weight_small':weight,'small_magnitude':small,'t':temperature,
                      'source_envelope':source,'segregated_policy':policy,
                      'child_values':[plus,minus],'gap_vs_source_upper':policy-source['upper']}
                except (ValueError,OverflowError,ZeroDivisionError) as error:
                    record={'weight_small':weight,'small_magnitude':small,'t':temperature,'error':str(error)}
                records.append(record);print(json.dumps(record),flush=True)
                args.output.write_text(json.dumps({'records':records,'qualification':'floating diagnostics; no global non-falsification theorem'},indent=2)+'\n')


if __name__=='__main__':main()
