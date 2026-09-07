"""Exact simple-column linear-code counterexample to generic deep-hole balance."""
from collections import deque
import json


def main():
    columns = [29, 63, 57, 12, 41, 49, 19, 39, 13, 40, 42, 27]
    assert len(set(columns)) == len(columns) and 0 not in columns
    anti = 0
    for value in columns:
        anti ^= value
    distances = [-1] * 64
    leaders = [None] * 64
    distances[0], leaders[0] = 0, 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for index, column in enumerate(columns):
            y = x ^ column
            if distances[y] < 0:
                distances[y] = distances[x] + 1
                leaders[y] = leaders[x] ^ (1 << index)
                queue.append(y)
    assert min(distances) >= 0
    # Independent exhaustive check, not only graph traversal.
    exhaustive = [13] * 64
    for word in range(1 << 12):
        syndrome = 0
        for i, column in enumerate(columns):
            if (word >> i) & 1:
                syndrome ^= column
        exhaustive[syndrome] = min(exhaustive[syndrome], bin(word).count("1"))
    assert exhaustive == distances
    radius = max(min(distances[x], distances[x ^ anti]) for x in range(64))
    assert anti == 47 and radius == 2
    assert distances[2] == 2 and distances[2 ^ anti] == 4
    print(json.dumps({"columns": columns, "antipodal_syndrome": anti,
                      "augmented_covering_radius": radius,
                      "deep_syndrome": 2, "two_half_distances": [2, 4],
                      "deep_word": leaders[2], "all_distances": distances,
                      "exhaustive_words_checked": 4096}, indent=2))


if __name__ == "__main__":
    main()
