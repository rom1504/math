"""Exact scalable cosquare cap separation; no numerical solver used."""
from fractions import Fraction
import itertools
import json
import numpy as np
from fresh_limit_cosquare_certificate import exact_positive_ldl


C = np.array([
 [1,-1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1],
 [-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,1,1],
 [-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,1,1],
 [-1,-1,1,1,-1,1,1,-1,-1,1,-1,1,1,1],
 [-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,1],
 [1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1,1],
 [1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,1,1],
 [-1,-1,-1,-1,1,1,1,1,1,-1,1,-1,1,1],
 [1,1,1,-1,-1,-1,-1,1,1,1,-1,-1,1,1],
 [-1,-1,-1,1,1,-1,1,-1,1,1,1,-1,1,1],
 [1,1,1,-1,-1,-1,-1,1,-1,1,1,-1,1,1],
 [1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,1,1],
 [1,1,1,1,1,1,1,1,1,1,1,1,-1,-1],
 [1,1,1,1,1,1,1,1,1,1,1,1,-1,-1]],dtype=np.int64)

SCALE=1000000
P_SCALED=np.array([
 [3025085,889307,-154342,225446,-921348,-452072,347319,1105907,144018,-23561,144018,553457,-197227,-197227],
 [889307,4097513,1056413,-29912,1257882,-265162,-742676,-4215,-613813,-154425,-613813,305727,-249352,-249352],
 [-154342,1056413,4097519,305719,-265240,1257875,-742676,-4311,-613765,889328,-613765,-29907,-249356,-249356],
 [225446,-29912,305719,2667345,509230,68139,1083873,-354386,-480490,553394,-480490,757,-2,-2],
 [-921348,1257882,-265240,509230,3642988,-592497,343089,-175806,806964,-452177,806964,68081,-228368,-228368],
 [-452072,-265162,1257875,68139,-592497,3642868,343133,-176029,807015,-921365,807015,509128,-228345,-228345],
 [347319,-742676,-742676,1083873,343089,343133,2607799,-609663,40553,347276,40553,1083874,0,0],
 [1105907,-4215,-4311,-354386,-175806,-176029,-609663,4610482,219349,1105880,219349,-354553,-383600,-383600],
 [144018,-613813,-613765,-480490,806964,807015,40553,219349,4055601,144097,1679917,-480409,-410034,-410034],
 [-23561,-154425,889328,553394,-452177,-921365,347276,1105880,144097,3025237,144097,225539,-197247,-197247],
 [144018,-613813,-613765,-480490,806964,807015,40553,219349,1679917,144097,4055601,-480409,-410034,-410034],
 [553457,305727,-29907,757,68081,509128,1083874,-354553,-480409,225539,-480409,2667333,2,2],
 [-197227,-249352,-249356,-2,-228368,-228345,0,-383600,-410034,-197247,-410034,2,3209746,1914269],
 [-197227,-249352,-249356,-2,-228368,-228345,0,-383600,-410034,-197247,-410034,2,1914269,3209746]],dtype=np.int64)
P_SCALED += 100*np.eye(14,dtype=np.int64)


def verify():
    n=len(C)
    cp=C.copy();cp[:12,:12]*=-1
    assert np.array_equal(C,C.T) and np.array_equal(cp,cp.T)
    assert np.all(np.abs(C)==1) and np.all(C[:12,:12].sum(axis=0)==0)
    assert np.array_equal(C@C,cp@cp)
    spins=np.array(list(itertools.product((-1,1),repeat=n)),dtype=np.int64)
    caps=[]
    for a in [C,cp]:
        energies=np.sum((spins@a)*spins,axis=1)
        index=int(np.argmax(np.abs(energies)))
        bilinear_index=int(np.argmax(np.abs(spins@a).sum(axis=1)))
        f=spins[bilinear_index]
        g=np.where(a@f>=0,1,-1)
        caps.append(dict(q_twice=int(abs(energies[index])),
                         beta=int(np.max(np.abs(spins@a).sum(axis=1))),
                         quadratic_witness=spins[index].tolist(),
                         bilinear_left=f.tolist(),bilinear_right=g.tolist(),
                         signed_twice_energy=int(energies[index])))
    exact_positive_ldl(P_SCALED-SCALE*C)
    exact_positive_ldl(P_SCALED+SCALE*C)
    maxp=int(np.max(np.sum((spins@P_SCALED)*spins,axis=1)))
    upper=Fraction(maxp,2*SCALE)
    assert caps[0]['beta']==72 and caps[1]['beta']==80
    assert upper<40
    h4=np.ones((4,4),dtype=np.int64)-2*np.eye(4,dtype=np.int64)
    f=np.array(caps[1]['bilinear_left'],dtype=np.int64)
    g=np.array(caps[1]['bilinear_right'],dtype=np.int64)
    z=np.concatenate([f,f,g,g])
    lifted_twice=int(z @ np.kron(h4,cp) @ z)
    assert lifted_twice==640
    # Exact one-step conversion: Q(H4 tensor Cprime)>=4 beta(Cprime).
    print(json.dumps(dict(order=n,spin_inputs=2**n,cosquare=True,full_caps=caps,
                         P_plus_minus_C_positive=True,T_upper=str(upper),
                         T_upper_decimal=str(float(upper)),
                         outer4_q_witness=lifted_twice//2,
                         outer4_spin_witness=z.tolist(),
                         normalized_outer_lower=40,unnormalized_gap=str(40-upper),
                         scalable_separation_established=True),sort_keys=True))


if __name__=='__main__':
    verify()
