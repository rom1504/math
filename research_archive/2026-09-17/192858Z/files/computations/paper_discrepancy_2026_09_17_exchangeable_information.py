"""Exhaustive finite replay of exchangeable-row information regularization.

Graph energies and pointwise comparisons are exact integers. Entropies
and uniform-tie probabilities are evaluated in floating point; this is
an independent replay, not the proof of the information theorem.
"""

from __future__ import annotations

import json
import math

import numpy as np


def spin_words(order):
    indices = np.arange(1 << order, dtype=np.int64)
    return (2*((indices[:, None] >> np.arange(order)) & 1)-1).astype(np.int64)


def energy_data(old, new):
    order = old+new
    words = spin_words(order)
    old_edges = [(i,j) for i in range(old) for j in range(i+1,old)]
    random_edges = [(i,j) for i in range(order) for j in range(i+1,order)
                    if j >= old]
    base = sum((words[:,i]*words[:,j] for i,j in old_edges),
               np.zeros(len(words),dtype=np.int64))
    products = np.asarray([words[:,i]*words[:,j] for i,j in random_edges],
                          dtype=np.int64).T
    return words, random_edges, base, products


def entropy(probabilities):
    positive = probabilities[probabilities > 0]
    return float(-np.dot(positive,np.log(positive)))


def replay(old, new):
    words, edges, base, products = energy_data(old,new)
    order = old+new
    total_graphs = 1 << len(edges)
    incident = [[edge for edge,(i,j) in enumerate(edges) if vertex in (i,j)]
                for vertex in range(old,order)]
    assert all(len(row) == order-1 for row in incident)
    assert max(sum(edge in row for row in incident) for edge in range(len(edges))) <= 2
    joint = np.zeros((new,1 << order,1 << (order-1)),dtype=float)
    cap_increments = np.zeros(new)
    selected_overlaps = np.zeros(new)
    cap_sum = 0
    for graph in range(total_graphs):
        signs = np.asarray([2*((graph >> i)&1)-1 for i in range(len(edges))])
        energies = base+products@signs
        cap = int(np.max(np.abs(energies)))
        cap_sum += cap
        grounds = np.flatnonzero(np.abs(energies) == cap)
        weight = 1.0/len(grounds)
        for j,row in enumerate(incident):
            row_pattern = sum(((graph >> edge)&1) << position
                              for position,edge in enumerate(row))
            joint[j,grounds,row_pattern] += weight
            contribution = products[:,row]@signs[row]
            deleted_cap = int(np.max(np.abs(energies-contribution)))
            difference = cap-deleted_cap
            assert difference >= 0
            assert np.all(difference <= np.abs(contribution[grounds]))
            cap_increments[j] += difference
            selected_overlaps[j] += float(np.mean(np.abs(contribution[grounds])))
    joint /= total_graphs
    cap_increments /= total_graphs
    selected_overlaps /= total_graphs
    zlaw = joint[0].sum(axis=1)
    hz = entropy(zlaw)
    informations = []
    for j in range(new):
        assert np.max(np.abs(joint[j].sum(axis=1)-zlaw)) < 1e-12
        rowlaw = joint[j].sum(axis=0)
        assert np.max(np.abs(rowlaw-1/(1 << (order-1)))) < 1e-12
        informations.append(hz+entropy(rowlaw)-entropy(joint[j].ravel()))
    assert max(informations)-min(informations) < 1e-11
    assert sum(informations) <= 2*hz+1e-11
    length = order-1
    fixed_mean = sum(math.comb(length,k)*abs(2*k-length)
                     for k in range(length+1))/(1 << length)
    selection_bounds = fixed_mean+np.sqrt(2*length*np.asarray(informations))
    assert np.all(selected_overlaps <= selection_bounds+1e-11)
    final_bound = math.sqrt(length)+2*math.sqrt(length*order*math.log(2)/new)
    assert np.all(cap_increments <= final_bound)

    # Independently average a fresh row over all (new-1)-vertex parents.
    parent_words, parent_edges, parent_base, parent_products = energy_data(old,new-1)
    row_overlaps = parent_words@parent_words.T
    delta_sum = 0.0
    nearlevel_checks = 0
    parent_graphs = 1 << len(parent_edges)
    for graph in range(parent_graphs):
        signs = np.asarray([2*((graph >> i)&1)-1 for i in range(len(parent_edges))])
        energies = parent_base+parent_products@signs
        cap = int(np.max(np.abs(energies)))
        fresh_caps = np.max(np.abs(energies)[None,:]+np.abs(row_overlaps),axis=1)
        delta = float(np.mean(fresh_caps))-cap
        delta_sum += delta
        for window in (0,1,2,4):
            code = np.flatnonzero(cap-np.abs(energies) <= window)
            bernoulli_width = float(np.mean(np.max(row_overlaps[:,code],axis=1)))
            assert bernoulli_width <= window+delta+1e-12
            nearlevel_checks += 1
    conditional_mean = delta_sum/parent_graphs
    assert np.max(np.abs(cap_increments-conditional_mean)) < 1e-12
    return {"old_order":old,"new_vertices":new,"full_graphs":total_graphs,
            "optimizer_entropy":hz,"incident_row_information":informations,
            "read_two_sum_information":sum(informations),"read_two_bound":2*hz,
            "mean_cap_increment":float(cap_increments[0]),
            "independent_fresh_row_increment":conditional_mean,
            "mean_selected_row_overlap":float(selected_overlaps[0]),
            "information_selection_bound":float(selection_bounds[0]),
            "final_universal_increment_bound":final_bound,
            "all_nearlevel_checks":nearlevel_checks}


def main():
    cases = [replay(old,new) for old,new in ((2,1),(2,2),(2,3),(3,2),(3,3))]
    print(json.dumps({"status":"PASS exhaustive finite replay",
                      "integer_energies_exact":True,"entropies_floating_point":True,
                      "total_full_signings":sum(case["full_graphs"] for case in cases),
                      "total_nearlevel_checks":sum(case["all_nearlevel_checks"] for case in cases),
                      "cases":cases},indent=2))


if __name__ == "__main__":
    main()
