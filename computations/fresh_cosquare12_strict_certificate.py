"""Exact scalable cosquare gap from a full symmetric order-12 seed pair.

No optimizer is used.  Integer exhaustive cube checks and Fraction LDL
verify the certificate suggested by a separate diagnostic SDP.
"""
from fractions import Fraction
import itertools
import json
import numpy as np

from fresh_limit_cosquare_certificate import exact_positive_ldl


C=np.array([
 [1,1,1,1,-1,-1,-1,1,-1,-1,1,1],
 [1,1,-1,-1,-1,1,-1,-1,1,1,1,1],
 [1,-1,1,1,1,-1,1,-1,-1,-1,1,1],
 [1,-1,1,1,-1,1,-1,-1,1,-1,1,1],
 [-1,-1,1,-1,1,-1,-1,1,1,1,1,1],
 [-1,1,-1,1,-1,1,1,1,-1,-1,1,1],
 [-1,-1,1,-1,-1,1,1,1,-1,1,1,1],
 [1,-1,-1,-1,1,1,1,1,-1,-1,1,1],
 [-1,1,-1,1,1,-1,-1,-1,1,1,1,1],
 [-1,1,-1,-1,1,-1,1,-1,1,1,1,1],
 [1,1,1,1,1,1,1,1,1,1,-1,-1],
 [1,1,1,1,1,1,1,1,1,1,-1,-1]],dtype=np.int64)

SCALE=1000
P_SCALED=np.array([
 [3510,180,353,245,-8,457,70,-297,218,-732,-211,-131],
 [180,2827,-620,636,-189,167,88,2,575,562,-264,-214],
 [353,-620,3360,476,196,-174,-160,482,-39,121,-214,-128],
 [245,636,476,3780,87,-9,-35,-118,-591,-474,-66,-276],
 [-8,-189,196,87,2826,-583,629,132,475,662,-215,-263],
 [457,167,-174,-9,-583,3435,344,284,-302,376,-125,-217],
 [70,88,-160,-35,629,344,3682,436,-594,-464,-227,-115],
 [-297,2,482,-118,132,284,436,3552,-170,-306,-183,-159],
 [218,575,-39,-591,475,-302,-594,-170,4398,375,-289,-258],
 [-732,562,121,-474,662,376,-464,-306,375,4223,-257,-290],
 [-211,-264,-214,-66,-215,-125,-227,-183,-289,-257,4637,481],
 [-131,-214,-128,-276,-263,-217,-115,-159,-258,-290,481,4637]],dtype=np.int64)


def verify():
    cp=C.copy();cp[:10,:10]*=-1
    assert np.array_equal(C,C.T) and np.array_equal(cp,cp.T)
    assert np.all(np.abs(C)==1) and np.all(np.abs(cp)==1)
    assert np.all(C[:10,:10].sum(axis=0)==0)
    assert np.array_equal(C@C,cp@cp)
    spins=np.array(list(itertools.product((-1,1),repeat=12)),dtype=np.int64)
    caps=[]
    for a in (C,cp):
        energy=np.einsum('bi,ij,bj->b',spins,a,spins)
        index=int(np.argmax(np.abs(energy)))
        beta=int(np.max(np.sum(np.abs(spins@a),axis=1)))
        caps.append(dict(q_twice=int(abs(energy[index])),
                         signed_q_twice=int(energy[index]),beta=beta,
                         witness=spins[index].tolist()))
    assert caps[0]['q_twice']==52 and caps[0]['beta']==52
    assert caps[1]['q_twice']==60 and caps[1]['beta']==60
    minus=exact_positive_ldl(P_SCALED-SCALE*C)
    plus=exact_positive_ldl(P_SCALED+SCALE*C)
    pvalues=np.einsum('bi,ij,bj->b',spins,P_SCALED,spins)
    upper=Fraction(int(np.max(pvalues)),2*SCALE)
    assert upper<30
    report=dict(order=12,all_sign_inputs=4096,cosquare=True,
                seed_caps=caps,majorant_positive=True,
                minimum_LDL_pivots=[str(min(minus)),str(min(plus))],
                T_upper=str(upper),gap_above_T=str(Fraction(30)-upper),
                regular_hadamard_lift_lower=30,
                scalable_regularized_gap_verified=True)
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    verify()
