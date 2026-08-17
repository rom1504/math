#!/usr/bin/env python3
"""Exact checks for Theorem 17.3b's static-cover/dynamic-memory gap."""

def internal_cover_size(n, k):
    """Greedy optimum for equally spaced points and radius k."""

    uncovered = 0
    centers = []
    while uncovered <= n:
        center = min(uncovered + k, n)
        centers.append(center)
        uncovered = center + k + 1
    return len(centers)


def exact_predictor(n, k):
    """Transient states followed by the optimal midpoint-decoded sink."""

    cutoff = n - 2 * k
    sink = cutoff
    transition = {
        state: state + 1 for state in range(cutoff)
    }
    transition[sink] = sink
    output = {
        state: state for state in range(cutoff)
    }
    output[sink] = n - k
    return cutoff, transition, output


def verify_instance(n, k):
    assert 1 <= 2 * k <= n
    expected_cover = (n + 1 + 2 * k) // (2 * k + 1)
    assert internal_cover_size(n, k) == expected_cover

    cutoff, transition, output = exact_predictor(n, k)
    expected_memory = n - 2 * k + 1
    assert len(output) == expected_memory

    checks = 0
    for physical_start in range(n + 1):
        predictor = min(physical_start, cutoff)
        physical = physical_start
        for _ in range(2 * n + 3):
            assert abs(output[predictor] - physical) <= k
            physical = min(physical + 1, n)
            predictor = transition[predictor]
            checks += 1

    # The lower-bound inequality forbids a repeated predictor state whose
    # first occurrence lies before n-2k.
    for first_repeat in range(cutoff):
        assert n - first_repeat > 2 * k
        checks += 1

    return expected_cover, expected_memory, checks


def main():
    total = 0
    largest_gap = None
    for n in range(2, 81):
        for k in range(1, n // 2 + 1):
            cover, memory, checks = verify_instance(n, k)
            total += checks
            gap = (memory / cover, n, k, cover, memory)
            if largest_gap is None or gap > largest_gap:
                largest_gap = gap

    _ratio, n, k, cover, memory = largest_gap
    print(f"exact trajectory/error checks: {total}")
    print(
        "largest tested memory/cover gap: "
        f"n={n}, k={k}, static_cover={cover}, dynamic_memory={memory}"
    )


if __name__ == "__main__":
    main()
