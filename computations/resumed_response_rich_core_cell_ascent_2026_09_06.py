#!/usr/bin/env python3
"""Floating feasible policy ascent in the SAME actual rich-core (V,Z) frame.

The old response remains an affine baseline with its true inverse Gram.
Two-dimensional cell Jensen is stronger than the preceding weighted-V
certificate. Every iterate is a conditional feasible pair realized by
independent even gates; no new Gaussian coordinate is created here.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr

from resumed_response_rich_core_birth_2026_09_06 import run, psi
from resumed_response_scalar_full_optimizer_2026_09_06 import hermites, phi


def interval_moments(left,right,degree):
    values = np.zeros((len(left),degree+1))
    values[:,0] = np.where(left>=0,ndtr(-left)-ndtr(-right),ndtr(right)-ndtr(left))
    hl,hr = hermites(left,degree-1),hermites(right,degree-1)
    for k in range(1,degree+1):
        values[:,k] = (phi(left)*hl[:,k-1]-phi(right)*hr[:,k-1])/np.sqrt(k)
    return values


class CellObjective:
    def __init__(self,context,policy,z_mesh=64,levels=24):
        self.context = context
        self.levels = levels
        self.left = np.array([float(Fraction(row["left"])) for row in policy])
        self.right = np.array([float(Fraction(row["right"])) for row in policy])
        lower = np.array([float(Fraction(row["lower_Z"])) for row in policy])
        upper = np.array([float(Fraction(row["upper_Z"])) for row in policy])
        self.zendpoints = np.unique(np.r_[np.linspace(-8,8,16*z_mesh+1),lower,upper])
        zl,zr = self.zendpoints[:-1],self.zendpoints[1:]
        zmiddle = (zl+zr)/2
        core = context["core"]
        vm = interval_moments(self.left,self.right,core.degree)
        zm = interval_moments(zl,zr,levels-1)
        self.weight = 2*vm[:,0,None]*zm[None,:,0]
        self.vmean = vm[:,1]/vm[:,0]
        self.zmean = zm[:,1]/zm[:,0]
        self.oldH = (self.right<=context["alpha"]+1e-12).astype(float)[:,None]
        binm = vm@core.mcoef/vm[:,0]
        self.g = (vm@context["conditional_coefficients"][:,:levels])@zm.T / (vm[:,0,None]*zm[None,:,0])
        poly = context["polynomial"]
        polyh = np.array([poly[0]+poly[1]+3*poly[2]+15*poly[3],
                          np.sqrt(2)*(poly[1]+6*poly[2]+45*poly[3]),
                          np.sqrt(24)*(poly[2]+15*poly[3]),np.sqrt(720)*poly[3]])
        meanpoly = vm[:,:7:2]@polyh/vm[:,0]
        self.f = -.036*binm+self.oldH[:,0]*meanpoly
        self.u = context["tail_scale"]*binm+self.oldH[:,0]*(context["center_linear"]*binm+context["center_constant"])
        self.hZ = (self.f[:,None]-context["A"]*self.g)/context["nu"]
        self.Kold = context["lambda"]*self.g+context["gamma"]*self.u[:,None]
        self.oldnorm = 1-context["mass"]-context["old_variance"]
        self.old_g_inner = context["lambda"]+context["gamma"]*context["c"]
        self.old_z_inner = context["gamma"]*(context["uf"]-context["A"]*context["c"])/context["nu"]
        self.tail_p = 2*ndtr(-2)+(1-2*ndtr(-2))*2*ndtr(-8)
        self.tail_aV = 2*phi(2)
        self.tail_aZ = (1-2*ndtr(-2))*2*phi(8)
        self.initial_f = np.where(zmiddle[None,:]>upper[:,None],1.,
                                  np.where(zmiddle[None,:]<lower[:,None],-1.,0.))
        self.initial_q = self.initial_f*self.initial_f

    def moments(self,f,q):
        wf = self.weight*f
        return np.array([self.tail_p+np.sum(self.weight*q),
                         self.tail_aV+np.sum(wf*self.vmean[:,None]),
                         self.tail_aZ+np.sum(wf*self.zmean[None,:])])

    def evaluate(self,f,q,theta,details=False):
        p,av,az = self.moments(f,q)
        K = (1-theta)*self.Kold+theta*(av*self.g+az*self.hZ)
        norm = ((1-theta)**2*self.oldnorm
                +2*theta*(1-theta)*(av*self.old_g_inner+az*self.old_z_inner)
                +theta*theta*(av*av+az*az))
        t = (1-theta)*(1-self.context["mass"])+theta*p-norm
        if t <= 0:
            raise ValueError("Nonpositive true residual variance")
        H = (1-theta)*self.oldH+theta*(1-q)
        values = psi(K,t)
        value = float(np.sum(self.weight*H*values))
        if not details:
            return value
        B = float(np.sum(self.weight*H*phi(K/np.sqrt(t))))/np.sqrt(t)
        signs = 2*ndtr(K/np.sqrt(t))-1
        rV = (float(np.sum(self.weight*H*signs*self.g))
              -2*B*((1-theta)*self.old_g_inner+theta*av))
        rZ = (float(np.sum(self.weight*H*signs*self.hZ))
              -2*B*((1-theta)*self.old_z_inner+theta*az))
        score = rV*self.vmean[:,None]+rZ*self.zmean[None,:]
        qstar = (np.abs(score)>values-B).astype(float)
        fstar = np.sign(score)*qstar
        gap = theta*float(np.sum(self.weight*(score*(fstar-f)-(values-B)*(qstar-q))))
        return value,{"variance":t,"B":B,"rV":rV,"rZ":rZ,"gap":gap,
                      "moments":[p,av,az],"qstar":qstar,"fstar":fstar}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection",type=Path,default=Path("computations/results/resumed_response_conditional_v_certificate_2026_09_06.json"))
    parser.add_argument("--policy",type=Path,default=Path("computations/results/resumed_response_rich_core_rectangle_birth_policy_2026_09_06.json"))
    parser.add_argument("--output",type=Path)
    parser.add_argument("--iterations",type=int,default=15)
    parser.add_argument("--z-mesh",type=int,default=32)
    parser.add_argument("--levels",type=int,default=24)
    args = parser.parse_args()
    source = json.loads(args.projection.read_text())
    context = run(source,nodes=384,polynomial_surrogate=True,return_context=True)
    policy = json.loads(args.policy.read_text())["first"]["rectangle_policy"]
    obj = CellObjective(context,policy,args.z_mesh,args.levels)
    f,q = obj.initial_f.copy(),obj.initial_q.copy()
    theta = .437
    records=[]
    policies=[]
    for iteration in range(args.iterations+1):
        opt_theta = minimize_scalar(lambda th:-obj.evaluate(f,q,th),bounds=(0,1),method="bounded",
                                    options={"xatol":2e-10})
        theta = float(opt_theta.x)
        value,info = obj.evaluate(f,q,theta,True)
        record={key:info[key] for key in ["variance","B","rV","rZ","gap","moments"]}
        record.update({"iteration":iteration,"theta":theta,"value":value})
        if iteration==args.iterations or info["gap"]<1e-10:
            records.append(record)
            print(json.dumps(record),flush=True)
            break
        fs,qs = info["fstar"],info["qstar"]
        eps=1e-5
        record["finite_difference_gap"]=(obj.evaluate((1-eps)*f+eps*fs,(1-eps)*q+eps*qs,theta)-value)/eps
        opt = minimize_scalar(lambda step:-obj.evaluate((1-step)*f+step*fs,(1-step)*q+step*qs,theta),
                              bounds=(0,1),method="bounded",options={"xatol":2e-10})
        step = float(opt.x)
        record["step"] = step
        record["line_gain"] = -float(opt.fun)-value
        row_runs=[]
        for row in fs.astype(int):
            ends=np.r_[np.flatnonzero(np.diff(row))+1,len(row)]
            starts=np.r_[0,ends[:-1]]
            row_runs.append([[int(end),int(row[start])] for start,end in zip(starts,ends)])
        policies.append({"weight_step":step,"row_runs":row_runs})
        f,q = (1-step)*f+step*fs,(1-step)*q+step*qs
        records.append(record)
        print(json.dumps(record),flush=True)
    result={"status":"floating_diagnostic_not_certificate","z_mesh":args.z_mesh,
            "Z_Hermite_levels":args.levels,"grid_shape":list(f.shape),
            "initial_old_conditional_2d_value":obj.evaluate(obj.initial_f,obj.initial_q,0),
            "initial_mixed_conditional_2d_value":obj.evaluate(obj.initial_f,obj.initial_q,.437),
            "history":records,"z_endpoints":obj.zendpoints.tolist(),
            "policies":policies}
    if args.output:
        args.output.write_text(json.dumps(result)+"\n")


if __name__=="__main__":
    main()
