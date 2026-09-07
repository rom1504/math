import itertools
import numpy as np
from functools import lru_cache


def matrices(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(edges)):
        a = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(edges):
            a[i, j] = a[j, i] = 1 if (bits >> k) & 1 else -1
        yield a


def states(n):
    # One representative modulo global negation.
    for bits in range(1 << max(0, n - 1)):
        x = np.ones(n, dtype=int)
        for i in range(n - 1):
            if (bits >> i) & 1:
                x[i] = -1
        yield x


def endpoint_data(a):
    vals = [(int(x @ a @ x), x.copy()) for x in states(len(a))]
    q = max(abs(v) for v, _ in vals) if vals else 0
    return q, [(1 if v == q else -1, x) for v, x in vals if abs(v) == q]


def normalize_at_endpoint(a, tau, w):
    return tau * (w[:, None] * a * w[None, :])


def node_layers(c, pred, sigma):
    """Current C has positive absolute ground 1; pred is predecessor spin.

    Returns carried deficit and the two shore layer payoffs. Empty shores have
    no directed layer and are represented by None.
    """
    n = len(c)
    q, eps = endpoint_data(c)
    assert int(np.ones(n, dtype=int) @ c @ np.ones(n, dtype=int)) == q
    h = int(pred @ c @ pred)
    gap = q - sigma * h
    s = np.flatnonzero(pred != pred[-1])
    # pred[-1] is the chosen global-sign anchor. Difference shore can be
    # complemented without changing the quadratic state.
    t = np.array([i for i in range(n) if i not in set(s)], dtype=int)
    out = []
    for core, peel in ((s, t), (t, s)):
        if len(core) == 0:
            out.append(None)
            continue
        x = c[np.ix_(core, core)]
        qx, endpoints = endpoint_data(x)
        inherited = pred[core]
        hx = int(inherited @ x @ inherited)
        delta = qx - sigma * hx
        if len(peel):
            cross = c[np.ix_(peel, core)] @ inherited
            visibility = 2 * int(np.abs(cross).sum())
        else:
            visibility = 0
        out.append({
            "core": core,
            "peel": peel,
            "a": x,
            "q": qx,
            "delta": delta,
            "ell": visibility - delta,
            "endpoints": endpoints,
            "inherited": inherited,
        })
    return q, gap, out


def next_node(branch):
    ans = []
    for tau, w in branch["endpoints"]:
        c1 = normalize_at_endpoint(branch["a"], tau, w)
        pred1 = w * branch["inherited"]
        # Preserving the inherited augmented orientation multiplies it by tau.
        sigma1 = tau
        q1, gap1, layers1 = node_layers(c1, pred1, sigma1)
        assert q1 == branch["q"] and gap1 == branch["delta"]
        ans.append((tau, w, gap1, layers1, c1, pred1, sigma1))
    return ans


def search(n):
    best = None
    count = 0
    for a in matrices(n):
        q, endpoints = endpoint_data(a)
        for tau0, w0 in endpoints:
            c = normalize_at_endpoint(a, tau0, w0)
            # Deduplicate gauges by requiring this normalization only once.
            if not np.array_equal(c, a):
                continue
            if tau0 != 1 or not np.all(w0 == 1):
                continue
            for pred in states(n):
                q0, g, layers = node_layers(c, pred, 1)
                if g <= 0 or any(z is None for z in layers):
                    continue
                immediate = max(0, *(z["ell"] for z in layers))
                level2 = 0
                details = []
                for z in layers:
                    nz = next_node(z)
                    branch_min = None
                    for row in nz:
                        ls = row[3]
                        m = max(0, *(u["ell"] for u in ls if u is not None))
                        branch_min = m if branch_min is None else min(branch_min, m)
                    level2 = max(level2, branch_min or 0)
                    details.append((z, nz, branch_min))
                ratio = max(immediate, level2) / g
                if best is None or ratio < best[0]:
                    best = (ratio, c.copy(), pred.copy(), q0, g, layers, details)
                count += 1
    print("n", n, "nodes", count, "best ratio", None if best is None else best[0])
    if best:
        ratio, c, pred, q, g, layers, details = best
        print("q,g,pred", q, g, pred.tolist())
        print(c.tolist())
        print("level1", [(z["q"], z["delta"], z["ell"], z["core"].tolist()) for z in layers])
        for z, nz, bm in details:
            print(" branch", z["core"].tolist(), "minmax", bm)
            for tau, w, gap1, ls, c1, pred1, sigma1 in nz:
                print("  endpoint", tau, w.tolist(), "sigma", sigma1, "gap", gap1,
                      "layers", [None if u is None else (u["q"], u["delta"], u["ell"], u["core"].tolist()) for u in ls])


def search_four_payoffs(n):
    best = None
    best_sum = None
    count = 0
    for a in matrices(n):
        q, endpoints = endpoint_data(a)
        if int(np.ones(n, dtype=int) @ a @ np.ones(n, dtype=int)) != q:
            continue
        for pred in states(n):
            q0, g, layers = node_layers(a, pred, 1)
            if g <= 0 or any(z is None for z in layers):
                continue
            vals = []
            for z in layers:
                h = int(z["inherited"] @ z["a"] @ z["inherited"])
                vis = z["ell"] + (z["q"] - h)
                vals.extend([vis - (z["q"] - h), vis - (z["q"] + h)])
            ratio = max(0, *vals) / g
            sum_ratio = (max(vals[0], vals[1]) + max(vals[2], vals[3])) / g
            if best is None or ratio < best[0]:
                best = (ratio, a.copy(), pred.copy(), q0, g, vals, layers)
            if best_sum is None or sum_ratio < best_sum[0]:
                best_sum = (sum_ratio, a.copy(), pred.copy(), q0, g, vals)
            count += 1
    print("four n", n, "nodes", count, "best", None if best is None else best[0])
    print("best sum ratio", None if best_sum is None else best_sum[0])
    if best:
        ratio, a, pred, q, g, vals, layers = best
        print("q,g,pred", q, g, pred.tolist(), "vals", vals)
        print(a.tolist())


def search_pure_orientation_child(n):
    for a in matrices(n):
        q, endpoints = endpoint_data(a)
        if int(np.ones(n, dtype=int) @ a @ np.ones(n, dtype=int)) != q:
            continue
        for pred in states(n):
            q0, g, layers = node_layers(a, pred, 1)
            if g <= 0:
                continue
            for z in layers:
                if z is None or z["q"] <= 0 or z["delta"] != 2 * z["q"]:
                    continue
                for tau, w in z["endpoints"]:
                    if tau == -1 and (np.array_equal(w, z["inherited"]) or
                                      np.array_equal(w, -z["inherited"])):
                        print("pure child n", n, "q,g,pred", q0, g, pred.tolist())
                        print(a.tolist())
                        print("core", z["core"].tolist(), "qchild", z["q"],
                              "delta", z["delta"], "ell", z["ell"])
                        return True
    print("no pure child", n)
    return False


def key_matrix(a):
    return tuple(int(a[i, j]) for i in range(len(a)) for j in range(i + 1, len(a)))


def from_key(n, key):
    a = np.zeros((n, n), dtype=int)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = key[k]
            k += 1
    return a


@lru_cache(None)
def min_reset_tree_cost_cached(n, akey, pkey):
    a = from_key(n, akey)
    pred = np.array(pkey, dtype=int)
    q, endpoints = endpoint_data(a)
    assert int(np.ones(n, dtype=int) @ a @ np.ones(n, dtype=int)) == q
    g = q - int(pred @ a @ pred)
    if g == 0 or n <= 1:
        return 0, ()
    _, _, layers = node_layers(a, pred, 1)
    total = 0
    certificate = []
    for z in layers:
        if z is None:
            continue
        best = None
        for tau, w in z["endpoints"]:
            child = normalize_at_endpoint(z["a"], tau, w)
            inherited = w * z["inherited"]
            delta = z["q"] - tau * int(z["inherited"] @ z["a"] @ z["inherited"])
            vis = 2 * int(np.abs(a[np.ix_(z["peel"], z["core"])] @ z["inherited"]).sum())
            ell = vis - delta
            subcost, subcert = min_reset_tree_cost_cached(
                len(child), key_matrix(child), tuple(int(x) for x in inherited))
            row = (max(0, ell) + subcost, tau, w.tolist(), delta, ell, subcert)
            if best is None or row[0] < best[0]:
                best = row
        total += best[0]
        certificate.append(best)
    return total, tuple(certificate)


def search_tree_cost(n):
    worst = None
    for a in matrices(n):
        q, endpoints = endpoint_data(a)
        if int(np.ones(n, dtype=int) @ a @ np.ones(n, dtype=int)) != q:
            continue
        for pred in states(n):
            g = q - int(pred @ a @ pred)
            if g <= 0:
                continue
            cost, cert = min_reset_tree_cost_cached(n, key_matrix(a), tuple(int(x) for x in pred))
            assert cost >= g
            ratio = cost / q
            if worst is None or ratio > worst[0]:
                worst = (ratio, cost, q, g, a.copy(), pred.copy(), cert)
    print("tree n", n, "worst cost/q", None if worst is None else worst[:4])
    if worst:
        print(worst[4].tolist(), "pred", worst[5].tolist())


if __name__ == "__main__":
    for n in range(3, 7):
        search(n)
