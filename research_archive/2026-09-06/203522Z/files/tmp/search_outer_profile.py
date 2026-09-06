from collections import defaultdict
from itertools import combinations, permutations


def hamming(a, b):
    return bin(a ^ b).count("1")


def permute_word(word, permutation, n):
    result = 0
    for old in range(n):
        if (word >> old) & 1:
            result |= 1 << permutation[old]
    return result


def transform_mask(mask, translation, permutation, n):
    result = 0
    for word in range(1 << n):
        if (mask >> word) & 1:
            image = permute_word(word ^ translation, permutation, n)
            result |= 1 << image
    return result


def canonical(mask, n, permutations_list):
    return min(
        transform_mask(mask, translation, permutation, n)
        for translation in range(1 << n)
        for permutation in permutations_list
    )


def outer_profile(mask, n):
    code = [word for word in range(1 << n) if (mask >> word) & 1]
    counts = [0] * (n + 1)
    for root in range(1 << n):
        counts[min(hamming(root, word) for word in code)] += 1
    return tuple(counts)


def inner_profile(mask, n):
    code = [word for word in range(1 << n) if (mask >> word) & 1]
    counts = [0] * (n + 1)
    for left in code:
        for right in code:
            counts[hamming(left, right)] += 1
    return tuple(counts)


def main():
    n = 4
    buckets = defaultdict(list)
    for code in combinations(range(1 << n), 4):
        mask = sum(1 << word for word in code)
        key = outer_profile(mask, n)
        buckets[key].append(mask)
    for profile, masks in buckets.items():
        if len(masks) > 1:
            by_inner = defaultdict(list)
            for mask in masks:
                by_inner[inner_profile(mask, n)].append(mask)
            if len(by_inner) > 1:
                print("outer", profile)
                for mask in masks[:6]:
                    code = [x for x in range(1 << n) if (mask >> x) & 1]
                    print(mask, code, "inner", inner_profile(mask, n))
                return
    print("no desired collision", len(buckets))


if __name__ == "__main__":
    main()
