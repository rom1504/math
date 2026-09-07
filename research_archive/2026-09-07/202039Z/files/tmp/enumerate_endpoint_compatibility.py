from itertools import product


A = (
    (0, 1, -1, 1, 1, 1, 1, -1, -1),
    (1, 0, 1, -1, -1, -1, 1, -1, -1),
    (-1, 1, 0, 1, 1, 1, 1, 1, -1),
    (1, -1, 1, 0, 1, 1, 1, -1, 1),
    (1, -1, 1, 1, 0, 1, -1, 1, -1),
    (1, -1, 1, 1, 1, 0, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, 0, 1, 1),
    (-1, -1, 1, -1, 1, -1, 1, 0, -1),
    (-1, -1, -1, 1, -1, -1, 1, -1, 0),
)
n = len(A)
edges = tuple((u, v) for u in range(n) for v in range(u + 1, n))


def spin_vectors(vertices):
    """Represent spins modulo global negation by fixing the first to +1."""
    for tail in product((-1, 1), repeat=len(vertices) - 1):
        yield dict(zip(vertices, (1, *tail)))


def signed_energy(vertices, spin):
    return sum(
        A[u][v] * spin[u] * spin[v]
        for pos, u in enumerate(vertices)
        for v in vertices[pos + 1 :]
    )


def child_codewords(deleted):
    vertices = tuple(v for v in range(n) if v != deleted)
    samples = tuple((signed_energy(vertices, x), x) for x in spin_vectors(vertices))
    norm = max(abs(value) for value, _ in samples)
    words = []
    for value, x in samples:
        if abs(value) != norm:
            continue
        orientation = 1 if value > 0 else -1
        word = {
            (u, v): orientation * x[u] * x[v]
            for pos, u in enumerate(vertices)
            for v in vertices[pos + 1 :]
        }
        words.append((orientation, word))
    assert len({tuple(word[e] for e in word) for _, word in words}) == len(words)
    return norm, tuple(words)


families = []
for i in range(n):
    norm, words = child_codewords(i)
    assert norm == 12
    families.append(words)
print("local family sizes", tuple(map(len, families)))


pair_cost = {}
for i in range(n):
    for j in range(i + 1, n):
        common = tuple(e for e in edges if i not in e and j not in e)
        pair_cost[i, j] = tuple(
            tuple(sum(wi[e] != wj[e] for e in common) for _, wj in families[j])
            for _, wi in families[i]
        )


global_words = []
for t in (-1, 1):
    for x in spin_vectors(tuple(range(n))):
        global_words.append(tuple(t * x[u] * x[v] for u, v in edges))
assert len(global_words) == 2**n
edge_index = {edge: k for k, edge in enumerate(edges)}


local_to_global_cost = []
for i, words in enumerate(families):
    local_edges = tuple(e for e in edges if i not in e)
    local_to_global_cost.append(
        tuple(
            tuple(
                sum(word[e] != global_word[edge_index[e]] for e in local_edges)
                for global_word in global_words
            )
            for _, word in words
        )
    )


min_gamma = None
min_gamma_indices = None
min_repair = None
min_repair_indices = None
count = 0
for choice in product(*(range(len(words)) for words in families)):
    count += 1
    gamma = sum(
        pair_cost[i, j][choice[i]][choice[j]]
        for i in range(n)
        for j in range(i + 1, n)
    )
    if min_gamma is None or gamma < min_gamma:
        min_gamma = gamma
        min_gamma_indices = choice

    repair = min(
        sum(local_to_global_cost[i][choice[i]][g] for i in range(n))
        for g in range(len(global_words))
    )
    if min_repair is None or repair < min_repair:
        min_repair = repair
        min_repair_indices = choice

print("families checked", count)
print("minimum Gamma", min_gamma, "choice", min_gamma_indices)
print("minimum repair", min_repair, "choice", min_repair_indices)
