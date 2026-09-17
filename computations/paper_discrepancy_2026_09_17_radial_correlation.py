"""Radial correlation-polytope LPs for archived finite minimizers.

Exact rational primal and pointwise dual certificates are reconstructed
and verified. Global optimality labels of the matrices are imported.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]


def integerize(values):
    denominator = math.lcm(*(value.denominator for value in values))
    return denominator,np.asarray([value.numerator*(denominator//value.denominator)
                                   for value in values],dtype=object)


def recover_primal(result, features, edge_signs, direction, module):
    support = np.flatnonzero(result.x[:-1] > 1e-9)
    for bound in (1000,1000000,1000000000):
        values = [F(float(result.x[i])).limit_denominator(bound) for i in support]
        alpha = F(float(result.x[-1])).limit_denominator(bound)
        if sum(values) != 1 or any(value < 0 for value in values):
            continue
        den, nums = integerize(values+[alpha])
        if np.array_equal(features[support].T.astype(object)@nums[:-1],
                          direction*edge_signs.astype(object)*nums[-1]):
            return support,values,alpha
    rows = [[1]*len(support)+[0]]
    rows += [features[support,e].tolist()+[-direction*int(edge_signs[e])]
             for e in range(features.shape[1])]
    values = module.exact_linear_solution(rows,[1]+[0]*features.shape[1],
                                         np.r_[result.x[support],result.x[-1]])
    assert values is not None and all(value >= 0 for value in values)
    return support,values[:-1],values[-1]


def recover_dual(result, features, edge_signs, direction, alpha):
    # Nonnegative polynomial c+b·z has c=alpha and -direction*A·b=1.
    for bound in (1000,1000000,1000000000):
        dual = [F(float(-value)).limit_denominator(bound)
                for value in result.eqlin.marginals]
        if dual[0] != alpha:
            continue
        den, nums = integerize(dual)
        normal = -direction*sum(int(a)*b for a,b in zip(edge_signs,dual[1:]))
        if normal != 1:
            continue
        margins = nums[0]+features.astype(object)@nums[1:]
        if min(margins) >= 0:
            return dual,int(min(margins)),den
    raise AssertionError("Exact rational dual reconstruction failed")


def certify_mixture(result, block, overlaps, words, ground, alpha_plus, alpha_minus,
                    edge_signs, positive_weight, full, module, radial_faces):
    count = len(words)
    eq_count = len(block)
    rhs_exact = [F(1)]+[alpha_plus*int(a) for a in edge_signs]
    rhs_exact += [F(1)]+[-alpha_minus*int(a) for a in edge_signs]
    support = [np.flatnonzero(result.x[offset:offset+count] > 1e-8)
               for offset in (0,count)]
    primal = None
    for bound in (1000,1000000,1000000000):
        weights = [[F(float(result.x[offset+i])).limit_denominator(bound) for i in indices]
                   for offset,indices in zip((0,count),support)]
        valid = True
        for side,(indices,values) in enumerate(zip(support,weights)):
            den,nums = integerize(values)
            actual = block[:,indices].astype(object)@nums
            target = rhs_exact[side*eq_count:(side+1)*eq_count]
            valid = valid and all(F(int(a),den) == b for a,b in zip(actual,target))
        if not valid:
            continue
        weighted = [[scale*w for w in values]
                    for scale,values in zip((positive_weight,1-positive_weight),weights)]
        den,nums = integerize(weighted[0]+weighted[1])
        response = np.column_stack((overlaps[:,support[0]],overlaps[:,support[1]])).astype(object)@nums
        value = F(int(max(response)),den)
        primal = (weights,value)
        break
    if primal is None:
        left_count = len(support[0])
        columns = left_count+len(support[1])+1
        rows = []
        targets = []
        for side in range(2):
            for row_index in range(eq_count):
                row = [F(0)]*columns
                offset = 0 if side == 0 else left_count
                for j,index in enumerate(support[side]):
                    row[offset+j] = F(int(block[row_index,index]))
                target = rhs_exact[side*eq_count+row_index]
                den,nums = integerize(row+[target])
                rows.append([int(v) for v in nums[:-1]])
                targets.append(int(nums[-1]))
        approximate_responses = (float(positive_weight)*overlaps@result.x[:count]
                                 +float(1-positive_weight)*overlaps@result.x[count:2*count])
        tight = np.flatnonzero(np.abs(approximate_responses-result.x[-1]) < 1e-7)
        for index in tight:
            row = [positive_weight*int(v) for v in overlaps[index,support[0]]]
            row += [(1-positive_weight)*int(v) for v in overlaps[index,support[1]]]+[F(-1)]
            den,nums = integerize(row)
            rows.append([int(v) for v in nums])
            targets.append(0)
        approximate = np.r_[result.x[support[0]],result.x[count+support[1]],result.x[-1]]
        values = module.exact_linear_solution(rows,targets,approximate)
        assert values is not None and min(values)>=0
        weights = [values[:left_count],values[left_count:-1]]
        weighted = [[scale*w for w in local]
                    for scale,local in zip((positive_weight,1-positive_weight),weights)]
        den,nums = integerize(weighted[0]+weighted[1])
        response = np.column_stack((overlaps[:,support[0]],overlaps[:,support[1]])).astype(object)@nums
        value = F(int(max(response)),den)
        assert value == values[-1]
        primal = (weights,value)
    assert primal is not None,(words.shape[1],"Secondary primal rational reconstruction failed")
    weights,value = primal
    exact_dual = None
    for bound in (1000,1000000,1000000000):
        dual = [F(float(z)).limit_denominator(bound) if abs(z)>1e-9 else F(0)
                for z in result.eqlin.marginals]
        code_weights = [F(float(-z)).limit_denominator(bound) if abs(z)>1e-9 else F(0)
                        for z in result.ineqlin.marginals]
        if min(code_weights) < 0 or sum(code_weights)>1:
            continue
        objective = sum(a*b for a,b in zip(dual,rhs_exact))
        if objective != value:
            continue
        valid = True
        for side,scale in enumerate((positive_weight,1-positive_weight)):
            local = dual[side*eq_count:(side+1)*eq_count]
            den,nums = integerize(local+[scale*w for w in code_weights])
            margins = overlaps.T.astype(object)@nums[eq_count:]-block.T.astype(object)@nums[:eq_count]
            valid = valid and min(margins)>=0
        if valid:
            exact_dual = (dual,code_weights)
            break
    face_certificate = None
    if exact_dual is None:
        uniform_responses = overlaps.sum(axis=0)
        face_minima = [F(int(np.min(uniform_responses[face])),len(ground)) for face in radial_faces]
        face_value = positive_weight*face_minima[0]+(1-positive_weight)*face_minima[1]
        if face_value == value:
            face_certificate = {"ground_law":"uniform on the complete projective absolute ground code",
                                "positive_radial_zero_face_minimum":str(face_minima[0]),
                                "negative_radial_zero_face_minimum":str(face_minima[1]),
                                "radial_zero_face_sizes":[int(np.count_nonzero(face)) for face in radial_faces]}
    lower_value = value
    rounded_fallback = False
    if exact_dual is None and face_certificate is None:
        # Preserve the failure of exact endpoint reconstruction, but obtain
        # an exact nearby dual by rational rounding and constant correction.
        precision = 10**12
        dual = [F(round(float(z)*precision),precision) for z in result.eqlin.marginals]
        code_weights = [F(max(0,round(float(-z)*precision)),precision)
                        for z in result.ineqlin.marginals]
        scale = max(F(1),sum(code_weights))
        dual = [z/scale for z in dual]
        code_weights = [z/scale for z in code_weights]
        for side,weight in enumerate((positive_weight,1-positive_weight)):
            local = dual[side*eq_count:(side+1)*eq_count]
            den,nums = integerize(local+[weight*w for w in code_weights])
            margins = overlaps.T.astype(object)@nums[eq_count:]-block.T.astype(object)@nums[:eq_count]
            correction = min(F(0),F(int(min(margins)),den))
            dual[side*eq_count] += correction
        for side,weight in enumerate((positive_weight,1-positive_weight)):
            local = dual[side*eq_count:(side+1)*eq_count]
            den,nums = integerize(local+[weight*w for w in code_weights])
            margins = overlaps.T.astype(object)@nums[eq_count:]-block.T.astype(object)@nums[:eq_count]
            assert min(margins)>=0
        lower_value = sum(a*b for a,b in zip(dual,rhs_exact))
        assert lower_value <= value and value-lower_value < F(1,1000000)
        exact_dual = (dual,code_weights)
        rounded_fallback = True
    report = {"optimized_ground_response_upper_exact":str(value),
              "optimized_ground_response_lower_exact":str(lower_value),
              "optimized_ground_response_normalized":float(value)/math.sqrt(words.shape[1]),
              "optimized_response_status":("exact rational interval; endpoint recovery failed, rounded dual corrected pointwise"
                                           if rounded_fallback else "exact two radial primal laws and all-sign-query dual agree")}
    if value == lower_value:
        report["optimized_ground_response_exact"] = str(value)
    if face_certificate is not None:
        report["optimized_face_lower_certificate"] = face_certificate
    if full:
        report["optimized_primal_laws"] = [
            [{"word":words[i].tolist(),"weight":str(weight)}
             for i,weight in zip(indices,values)] for indices,values in zip(support,weights)]
        if exact_dual is not None:
            report["optimized_dual_equality_coefficients"] = [str(v) for v in exact_dual[0]]
            report["optimized_dual_ground_law"] = [
                {"word":word.tolist(),"weight":str(weight)}
                for word,weight in zip(ground,exact_dual[1]) if weight]
    return report


def check_case(name, matrix, provenance, module, full, optimize_mixture):
    a = np.asarray(matrix,dtype=np.int64)
    n = len(a)
    edges = list(itertools.combinations(range(n),2))
    words = np.asarray([(1,)+rest for rest in itertools.product((-1,1),repeat=n-1)],
                       dtype=np.int64)
    features = np.asarray([words[:,i]*words[:,j] for i,j in edges],dtype=np.int64).T
    signs = np.asarray([a[i,j] for i,j in edges],dtype=np.int64)
    energies = features@signs
    cap = int(np.abs(energies).max())
    assert cap == provenance["minimum_cap_imported"]
    ground = words[np.abs(energies) == cap]
    spectral = np.linalg.eigvalsh(a)
    opnorm = float(np.max(np.abs(spectral)))
    reports = []
    laws = []
    radial_faces = []
    for direction in (1,-1):
        equality = np.vstack((np.r_[np.ones(len(words)),0],
                              np.column_stack((features.T,-direction*signs))))
        objective = np.r_[np.zeros(len(words)),-1.0]
        rhs = np.r_[1.0,np.zeros(len(edges))]
        result = linprog(objective,A_eq=equality,b_eq=rhs,bounds=(0,None),method="highs")
        assert result.success,result.message
        support,weights,alpha = recover_primal(result,features,signs,direction,module)
        dual,dual_margin,dual_denominator = recover_dual(result,features,signs,direction,alpha)
        dual_den,dual_nums = integerize(dual)
        radial_faces.append(dual_nums[0]+features.astype(object)@dual_nums[1:] == 0)
        den, nums = integerize(weights)
        response_numerators = np.abs(ground@words[support].T).astype(object)@nums
        responses = [F(int(value),den) for value in response_numerators]
        directional_cap = int(np.max(direction*energies))
        necessary = F(directional_cap,len(edges))
        min_eigenvalue = float(np.min(direction*spectral))
        gaussian_rho = min(1.0,-1/min_eigenvalue)
        gaussian_alpha = 2/math.pi*math.asin(gaussian_rho)
        report = {"direction":direction,"alpha_exact":str(alpha),
                  "t_alpha_sqrt_n":float(alpha)*math.sqrt(n),
                  "energy_upper_alpha_exact":str(necessary),
                  "energy_upper_saturated":alpha == necessary,
                  "directional_psd_upper_alpha":gaussian_rho,
                  "directional_gaussian_arcsine_alpha":gaussian_alpha,
                  "directional_gaussian_arcsine_t":gaussian_alpha*math.sqrt(n),
                  "operator_norm_linearized_gaussian_t":2*math.sqrt(n)/(math.pi*opnorm),
                  "primal_support_size":len(support),
                  "primal_ground_mean_response_min":str(min(responses)),
                  "primal_ground_mean_response_max":str(max(responses)),
                  "primal_ground_mean_response_max_normalized":float(max(responses))/math.sqrt(n),
                  "dual_minimum_integer_margin":dual_margin,
                  "dual_denominator":dual_denominator,
                  "certification":"exact radial law and full-cube dual agree"}
        if full:
            report["primal_law"] = [{"projective_index":int(i),"word":words[i].tolist(),
                                     "weight":str(weight)} for i,weight in zip(support,weights)]
            report["dual_constant"] = str(dual[0])
            report["dual_edge_coefficients"] = [str(value) for value in dual[1:]]
        reports.append(report)
        laws.append((alpha,responses))
    positive_weight = laws[1][0]/(laws[0][0]+laws[1][0])
    mixture_responses = [positive_weight*left+(1-positive_weight)*right
                         for left,right in zip(laws[0][1],laws[1][1])]
    mixture = {"positive_radial_law_weight_exact":str(positive_weight),
               "covariance":"I exactly (symmetrize each projective atom globally)",
               "primal_ground_mean_response_min_exact":str(min(mixture_responses)),
               "primal_ground_mean_response_max_exact":str(max(mixture_responses)),
               "primal_ground_mean_response_max_normalized":float(max(mixture_responses))/math.sqrt(n)}
    if optimize_mixture:
        count = len(words)
        block = np.vstack((np.ones(count,dtype=np.int64),features.T))
        equality = np.zeros((2*len(block),2*count+1))
        equality[:len(block),:count] = block
        equality[len(block):,count:2*count] = block
        rhs = np.r_[1,float(laws[0][0])*signs,1,-float(laws[1][0])*signs]
        overlaps = np.abs(ground@words.T)
        inequalities = np.column_stack((float(positive_weight)*overlaps,
                                        float(1-positive_weight)*overlaps,
                                        -np.ones(len(ground))))
        objective = np.r_[np.zeros(2*count),1.0]
        result = linprog(objective,A_eq=equality,b_eq=rhs,A_ub=inequalities,
                         b_ub=np.zeros(len(ground)),bounds=(0,None),method="highs")
        assert result.success,result.message
        mixture["optimized_at_fixed_radial_endpoints_numeric"] = float(result.x[-1])
        mixture["optimized_at_fixed_radial_endpoints_normalized_numeric"] = float(result.x[-1])/math.sqrt(n)
        mixture.update(certify_mixture(result,block,overlaps,words,ground,laws[0][0],laws[1][0],
                                       signs,positive_weight,full,module,radial_faces))
    return {"case":name,"n":n,"cap":cap,"spectral_min":float(spectral[0]),
            "spectral_max":float(spectral[-1]),"directions":reports,"isotropic_mixture":mixture}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full",action="store_true")
    parser.add_argument("--min-n",type=int,default=3)
    parser.add_argument("--max-n",type=int,default=14)
    parser.add_argument("--even-only",action="store_true")
    parser.add_argument("--optimize-mixture",action="store_true")
    args = parser.parse_args()
    path = ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py"
    spec = importlib.util.spec_from_file_location("stored_minimizers",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reports = [check_case(name,matrix,provenance,module,args.full,args.optimize_mixture)
               for name,matrix,provenance in module.cases()
               if args.min_n <= len(matrix) <= args.max_n and (not args.even_only or len(matrix)%2 == 0)]
    print(json.dumps({"status":"PASS exact radial primal and dual certificates",
                      "global_minimizer_provenance_imported":True,"cases":reports},indent=2))


if __name__ == "__main__":
    main()
