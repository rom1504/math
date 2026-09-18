from itertools import product

import numpy as np


A = [
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
]


def energy(matrix, spin):
    return sum(
        matrix[i][j] * spin[i] * spin[j]
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
    )


def spins(n):
    for tail in product((-1, 1), repeat=n - 1):
        yield (1, *tail)


EDGES = [(i, j) for i in range(9) for j in range(i + 1, 9)]


def word_mask(t, spin_by_vertex):
    mask = 0
    for bit, (i, j) in enumerate(EDGES):
        if i in spin_by_vertex and j in spin_by_vertex:
            value = t * spin_by_vertex[i] * spin_by_vertex[j]
            if value == -1:
                mask |= 1 << bit
    return mask


parent_values = [energy(A, x) for x in spins(9)]
parent_norm = max(map(abs, parent_values))
print("parent norm", parent_norm)
endpoint_row_profiles = []
endpoint_details = []
for x in spins(9):
    value = energy(A, x)
    if abs(value) != parent_norm:
        continue
    t = 1 if value > 0 else -1
    row_sums = tuple(
        sum(t * A[i][j] * x[i] * x[j] for j in range(9) if j != i)
        for i in range(9)
    )
    endpoint_row_profiles.append((t, row_sums))
    endpoint_details.append((max(row_sums), t, x, row_sums))
print("parent absolute endpoints modulo antipode", len(endpoint_row_profiles))
print("endpoint switched row profiles", sorted(set(endpoint_row_profiles)))
print("maximum row endpoint detail", max(endpoint_details))

local_words = []
for deleted in range(9):
    keep = [i for i in range(9) if i != deleted]
    child = [[A[i][j] for j in keep] for i in keep]
    records = []
    words = set()
    for x in spins(8):
        child_energy = energy(child, x)
        if abs(child_energy) != parent_norm:
            continue
        field = sum(A[deleted][j] * x[k] for k, j in enumerate(keep))
        records.append((child_energy, field))
        t = 1 if child_energy > 0 else -1
        words.add(word_mask(t, dict(zip(keep, x))))
    child_norm = max(abs(energy(child, x)) for x in spins(8))
    field_values = sorted({field for _, field in records})
    print(
        deleted + 1,
        "child norm",
        child_norm,
        "absolute grounds modulo antipode",
        len(records),
        "extension fields",
        field_values,
    )
    local_words.append(sorted(words))

global_words = []
for x in spins(9):
    spin_by_vertex = dict(enumerate(x))
    for t in (-1, 1):
        global_words.append(word_mask(t, spin_by_vertex))
assert len(set(global_words)) == 512

local_masks = []
for deleted in range(9):
    mask = 0
    for bit, (i, j) in enumerate(EDGES):
        if deleted not in (i, j):
            mask |= 1 << bit
    local_masks.append(mask)

choice_ranges = [range(len(words)) for words in local_words]
combinations = np.array(list(product(*choice_ranges)), dtype=np.uint8)
print("families", len(combinations), "choice counts", [len(w) for w in local_words])

def popcount(value):
    return bin(value).count("1")


repair_costs = []
for i, words in enumerate(local_words):
    repair_costs.append(np.array(
        [
            [popcount((global_word ^ local_word) & local_masks[i]) for global_word in global_words]
            for local_word in words
        ],
        dtype=np.uint8,
    ))

total_repairs = np.zeros((len(combinations), len(global_words)), dtype=np.uint16)
for i, costs in enumerate(repair_costs):
    total_repairs += costs[combinations[:, i]]
family_repairs = total_repairs.min(axis=1)
print("minimum repair R", int(family_repairs.min()), "families attaining", int(np.count_nonzero(family_repairs == family_repairs.min())))

gammas = np.zeros(len(combinations), dtype=np.uint16)
for i in range(9):
    for j in range(i + 1, 9):
        overlap_mask = local_masks[i] & local_masks[j]
        table = np.array(
            [
                [popcount((wi ^ wj) & overlap_mask) for wj in local_words[j]]
                for wi in local_words[i]
            ],
            dtype=np.uint8,
        )
        gammas += table[combinations[:, i], combinations[:, j]]
print("minimum Gamma", int(gammas.min()), "families attaining", int(np.count_nonzero(gammas == gammas.min())))
print("minimum (Gamma,R) pairs", sorted(set(zip(map(int, gammas), map(int, family_repairs))))[:20])
