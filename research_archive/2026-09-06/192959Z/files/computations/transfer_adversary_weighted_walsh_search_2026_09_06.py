"""Bounded lower-witness search; never an upper bound on stabilized R."""

import argparse
import json

import numpy as np


def walsh(values):
    result = values.copy()
    n = result.shape[-1]
    block = 1
    while block < n:
        reshaped = result.reshape(*result.shape[:-1], -1, 2, block)
        left = reshaped[..., 0, :].copy()
        right = reshaped[..., 1, :].copy()
        reshaped[..., 0, :] = left + right
        reshaped[..., 1, :] = left - right
        block *= 2
    return result


def score(inputs):
    spectra = walsh(inputs)
    return np.sum(np.abs(-spectra[:, 0] + 2 * spectra[:, 1])
                  + np.abs(2 * spectra[:, 0] + 4 * spectra[:, 1]), axis=1)


def improve(inputs, rounds=30):
    for _ in range(rounds):
        spectra = walsh(inputs)
        outputs = np.stack((-spectra[:, 0] + 2 * spectra[:, 1],
                            2 * spectra[:, 0] + 4 * spectra[:, 1]), axis=1)
        outputs = np.where(outputs >= 0, 1, -1)
        transformed = walsh(outputs)
        fields = np.stack((-transformed[:, 0] + 2 * transformed[:, 1],
                           2 * transformed[:, 0] + 4 * transformed[:, 1]), axis=1)
        following = np.where(fields >= 0, 1, -1)
        following = np.where(fields == 0, inputs, following)
        if np.array_equal(following, inputs):
            break
        inputs = following
    return inputs


def search(n, batch, cycles, seed):
    assert n > 0 and n & (n - 1) == 0
    rng = np.random.default_rng(seed)
    candidates = rng.choice([-1, 1], size=(batch, 2, n))
    if n.bit_length() % 2 == 1:
        coords = np.arange(n)
        quadratic = np.zeros(n, dtype=np.int64)
        for shift in range(0, n.bit_length() - 1, 2):
            quadratic ^= ((coords >> shift) & 1) & ((coords >> (shift + 1)) & 1)
        candidates[:batch // 2] = (1 - 2 * quadratic)[None, None, :]
    retained = improve(candidates)
    retained_scores = score(retained)
    best = int(np.max(retained_scores))
    best_input = retained[int(np.argmax(retained_scores))].copy()
    for cycle in range(cycles):
        rates = rng.choice([1 / n, 2 / n, 4 / n, 0.125, 0.25, 0.5], size=batch)
        trial = retained * (1 - 2 * (rng.random(retained.shape) < rates[:, None, None]))
        trial = improve(trial)
        trial_scores = score(trial)
        accept = (trial_scores >= retained_scores) | (rng.random(batch) < 0.03)
        retained[accept] = trial[accept]
        retained_scores[accept] = trial_scores[accept]
        if int(np.max(trial_scores)) > best:
            best = int(np.max(trial_scores))
            best_input = trial[int(np.argmax(trial_scores))].copy()
            print(json.dumps({"order": n, "cycle": cycle, "new_integer_score": best,
                              "normalized_lower_witness": best / (2 * n ** 1.5)}), flush=True)
    # Recalculate the winning integer witness independently.
    assert int(score(best_input[None])[0]) == best
    spectra = walsh(best_input)
    return {"order": n, "batch": batch, "cycles": cycles, "seed": seed,
            "integer_score": best, "normalized_lower_witness": best / (2 * n ** 1.5),
            "majorant_value": 5 / np.sqrt(2),
            "input_correlation": int(best_input[0] @ best_input[1]) / n,
            "reciprocal_product_l1_defect": float(np.mean(np.abs(spectra[0] * spectra[1] / n - 0.75))),
            "input": best_input.tolist()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=64)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--cycles", type=int, default=100)
    parser.add_argument("--seed", type=int, default=260906)
    args = parser.parse_args()
    print(json.dumps(search(args.order, args.batch, args.cycles, args.seed)))
