"""Bounded independent exact-energy search for the cosquare signing pair."""
import argparse
import contextlib
import io
import json
import time
import numpy as np

from fresh_limit_cosquare_block_probe import run


def pair():
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        run(10, 2, 10, 20260905)
    return json.loads(stream.getvalue())["witnesses"][0]


def search(power, seconds, seed, batch, target):
    witness = pair()
    core = np.array(witness["flipped"], dtype=np.int16)
    h4 = np.ones((4, 4), dtype=np.int16)-2*np.eye(4, dtype=np.int16)
    outer = np.ones((1, 1), dtype=np.int16)
    for _ in range(power):
        outer = np.kron(outer, h4)
    a = np.kron(outer, core)
    n = len(a)
    diag = np.diag(a)
    rng = np.random.default_rng(seed)
    base_spins = ((np.arange(1 << 12)[:, None] >> np.arange(12)) & 1)*2-1
    base_q = np.einsum('bi,ij,bj->b', base_spins, core, base_spins)
    base = base_spins[np.argmax(np.abs(base_q))].astype(np.int16)
    lifted = np.tile(base, len(outer))
    best = -1
    best_spin = None
    best_signed = None
    steps = 0
    start = time.monotonic()

    def check(x):
        nonlocal best, best_spin, best_signed
        values = np.einsum('bi,ij,bj->b', x.astype(np.int64), a, x.astype(np.int64))
        idx = np.argmax(np.abs(values))
        val = int(values[idx])
        if abs(val) > best:
            best = abs(val)
            best_signed = val
            best_spin = x[idx].copy()
            print(json.dumps(dict(elapsed=round(time.monotonic()-start, 3),
                                  n=n, best_q=best/2, signed_q=val/2,
                                  normalized_seed_q=best/(2*len(outer)**1.5))), flush=True)
        return best >= 2*target

    while time.monotonic()-start < seconds:
        x = rng.choice(np.array([-1, 1], dtype=np.int16), size=(batch, n))
        # Plant some starts around the exact seed witness and current best.
        x[:batch//4] = lifted
        flips = rng.random((batch//4, n)) < rng.choice([.03, .08, .15, .3])
        x[:batch//4] *= 1-2*flips.astype(np.int16)
        if best_spin is not None:
            x[batch//4:batch//2] = best_spin
            flips = rng.random((batch//4, n)) < rng.choice([.01, .03, .06, .12])
            x[batch//4:batch//2] *= 1-2*flips.astype(np.int16)
        orientation = np.ones(batch, dtype=np.int16)
        orientation[:batch//2] = -1
        field = x @ a
        temperatures = np.concatenate((np.geomspace(np.sqrt(n)*1.4, .2, 35), np.zeros(8)))
        for temp in temperatures:
            for k in rng.permutation(n):
                local = field[:, k]-diag[k]*x[:, k]
                gain = -2*orientation*x[:, k]*local
                do = gain > 0
                if temp:
                    do |= rng.random(batch) < np.exp(np.minimum(gain/temp, 0))
                delta = np.where(do, -2*x[:, k], 0).astype(np.int16)
                x[:, k] += delta
                field += delta[:, None]*a[k]
            if check(x):
                break
        # Exact 12-spin block responses escape single-spin local maxima.
        for _ in range(12):
            changed = False
            for block in rng.permutation(len(outer)):
                sl = slice(12*block, 12*(block+1))
                intra = int(outer[block, block])*core
                outside = field[:, sl]-x[:, sl] @ intra
                score = (base_spins @ outside.T).astype(np.int64)
                score += int(outer[block, block])*base_q[:, None]//2
                score *= orientation[None, :]
                selected = base_spins[np.argmax(score, axis=0)].astype(np.int16)
                delta = selected-x[:, sl]
                changed |= bool(np.any(delta))
                x[:, sl] = selected
                field += delta @ a[sl]
            if check(x) or not changed:
                break
        steps += batch
        if best >= 2*target:
            break
    assert best_spin is not None
    exact = int(best_spin.astype(np.int64) @ a @ best_spin.astype(np.int64))
    assert exact == best_signed and abs(exact) == best
    print(json.dumps(dict(method="annealed_batched_coordinate_search_exact_verification",
                          power=power, order=n, trials=steps, elapsed=time.monotonic()-start,
                          signed_quadratic_twice=exact, q=abs(exact)/2,
                          spin=best_spin.tolist(), core=core.tolist(),
                          target=target, target_reached=best >= 2*target)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--power', type=int, default=2)
    parser.add_argument('--seconds', type=float, default=90)
    parser.add_argument('--seed', type=int, default=2026090507)
    parser.add_argument('--batch', type=int, default=128)
    parser.add_argument('--target', type=int, default=1946)
    args = parser.parse_args()
    search(args.power, args.seconds, args.seed, args.batch, args.target)
