#!/usr/bin/env python3
"""Exact intrinsic restricted-contraction check for a width-two Ising strip.

The boundary has four spin assignments, so its projective max-plus state has
dimension three.  The certificate checked here has two parts:

1. A polyhedral gap box ``Y`` is forward invariant under two Ising columns.
2. One column has a common maximizing predecessor throughout ``Y`` and is
   therefore a projective reset.  Removing that reset edge from the legal
   control graph leaves an acyclic graph of longest path one.

Consequently every legal two-column product is 0-contractive on ``Y``.  The
check uses per-column inequalities and a graph longest-path computation; it
does not enumerate legal products.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product


SPINS = tuple(product((-1, 1), repeat=2))
ANCHOR = (1, 1)
ANCHOR_INDEX = SPINS.index(ANCHOR)


@dataclass(frozen=True)
class Column:
    horizontal: tuple[Q, Q]
    fields: tuple[Q, Q]
    vertical: Q = Q(0)


@dataclass(frozen=True)
class GapBox:
    """Bounds on v(++ )-v(x) for the three non-anchor states."""

    bounds: dict[tuple[int, int], tuple[Q, Q]]

    def contains(self, v):
        anchor = v[ANCHOR_INDEX]
        return all(lo <= anchor - v[SPINS.index(x)] <= hi
                   for x, (lo, hi) in self.bounds.items())

    def vertices(self):
        competitors = tuple(x for x in SPINS if x != ANCHOR)
        for choices in product((0, 1), repeat=len(competitors)):
            gaps = {
                x: self.bounds[x][side]
                for x, side in zip(competitors, choices)
            }
            yield tuple(Q(0) if x == ANCHOR else -gaps[x] for x in SPINS)


def local_energy(column, y):
    return (column.fields[0] * y[0]
            + column.fields[1] * y[1]
            + column.vertical * y[0] * y[1])


def kernel_entry(column, y, x):
    return (local_energy(column, y)
            + column.horizontal[0] * x[0] * y[0]
            + column.horizontal[1] * x[1] * y[1])


def transfer(v, column):
    return tuple(max(v[k] + kernel_entry(column, y, x)
                     for k, x in enumerate(SPINS))
                 for y in SPINS)


def unique_selectors(v, column):
    selectors = []
    for y in SPINS:
        values = tuple(v[k] + kernel_entry(column, y, x)
                       for k, x in enumerate(SPINS))
        winner = max(values)
        assert values.count(winner) == 1
        selectors.append(values.index(winner))
    return tuple(selectors)


def normalize(v):
    anchor = v[ANCHOR_INDEX]
    return tuple(z - anchor for z in v)


def hilbert(v, w):
    difference = tuple(x - y for x, y in zip(v, w))
    return Q(max(difference) - min(difference), 2)


def common_predecessor_threshold(column, x):
    """Worst kernel advantage of x over the anchor predecessor."""

    return max(kernel_entry(column, y, x)
               - kernel_entry(column, y, ANCHOR) for y in SPINS)


def certify_reset(box, column):
    """Prove the anchor predecessor maximizes every output row on the box."""

    for x, (lo, _hi) in box.bounds.items():
        assert lo >= common_predecessor_threshold(column, x)

    reset_image = normalize(tuple(kernel_entry(column, y, ANCHOR)
                                  for y in SPINS))
    assert box.contains(reset_image)

    # This finite check is diagnostic.  The inequalities above prove the
    # assertion on the entire continuum box.
    for v in box.vertices():
        assert normalize(transfer(v, column)) == reset_image
    return reset_image


def max_difference_bounds(column, y0, y1):
    """Bounds for max_x(v_x+K(y0,x))-max_x(v_x+K(y1,x))."""

    differences = tuple(kernel_entry(column, y0, x)
                        - kernel_entry(column, y1, x) for x in SPINS)
    return min(differences), max(differences)


def certify_global_image_in_box(box, column):
    """Prove that this column maps the entire projective space into box."""

    for y, (wanted_lo, wanted_hi) in box.bounds.items():
        lo, hi = max_difference_bounds(column, ANCHOR, y)
        assert wanted_lo <= lo <= hi <= wanted_hi

    # Exact vertex diagnostics for the restricted invariant set.
    for v in box.vertices():
        assert box.contains(transfer(v, column))


def longest_path_after_deleting_resets(vertices, edges):
    """Return the longest path in the nonreset graph, rejecting a cycle."""

    adjacency = {v: [] for v in vertices}
    indegree = {v: 0 for v in vertices}
    for source, target, is_reset in edges:
        if is_reset:
            continue
        adjacency[source].append(target)
        indegree[target] += 1

    queue = [v for v in vertices if indegree[v] == 0]
    distance = {v: 0 for v in vertices}
    seen = 0
    while queue:
        source = queue.pop()
        seen += 1
        for target in adjacency[source]:
            distance[target] = max(distance[target], distance[source] + 1)
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
    assert seen == len(vertices), "nonreset legal subgraph contains a cycle"
    return max(distance.values())


def main():
    # Coordinates are normalized at (++).  The three gap intervals are for
    # (-,+), (+,-), and (-,-), respectively.
    box = GapBox({
        (-1, 1): (Q(4), Q(20)),
        (1, -1): (Q(4), Q(20)),
        (-1, -1): (Q(8), Q(40)),
    })

    # A has global image in Y but retains a slope-one selector direction.
    active = Column(horizontal=(Q(4), Q(4)), fields=(Q(6), Q(6)))
    certify_global_image_in_box(box, active)

    # W is weak relative to every incoming gap in Y.  The common predecessor
    # (++ ) wins all four output maximizations, so W is constant projectively.
    weak = Column(horizontal=(Q(1), Q(1)), fields=(Q(3), Q(3)))
    reset_image = certify_reset(box, weak)

    # An exact witness shows that A's restricted Hilbert coefficient is one.
    v0 = (Q(-8), Q(-4), Q(-4), Q(0))
    v1 = (Q(-8), Q(-5), Q(-4), Q(0))
    assert box.contains(v0) and box.contains(v1)
    assert unique_selectors(v0, active) == tuple(range(len(SPINS)))
    assert unique_selectors(v0, weak) == (ANCHOR_INDEX,) * len(SPINS)
    assert hilbert(transfer(v0, active), transfer(v1, active)) == hilbert(v0, v1)

    # The legal language alternates A and W.  Deleting W leaves one edge, so
    # graph structure alone proves that every length-two path contains a reset.
    vertices = (0, 1)
    edges = ((0, 1, False), (1, 0, True))
    reset_free_horizon = longest_path_after_deleting_resets(vertices, edges)
    block_length = reset_free_horizon + 1
    assert reset_free_horizon == 1
    assert block_length == 2

    # Y lies in the Hilbert ball of radius 20 because its largest possible
    # anchored gap is 40.  Thus Theorem 17.2 applies with L=2 and rho=0.
    radius = Q(20)
    assert max(hi for _lo, hi in box.bounds.values()) / 2 == radius

    print("width-two Ising boundary states: 4")
    print(f"restricted Hilbert-ball radius: {radius}")
    print("active-column restricted coefficient: 1 (exact witness)")
    print(f"weak-column reset image: {reset_image}")
    print(f"nonreset longest legal path: {reset_free_horizon}")
    print(f"Theorem 17.2 contraction certificate: L={block_length}, rho=0")


if __name__ == "__main__":
    main()
