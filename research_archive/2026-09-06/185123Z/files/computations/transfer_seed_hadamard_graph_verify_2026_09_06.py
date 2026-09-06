"""Exact graph-contraction and injectivity checks for Hadamard compression.

No input/output files. Graph edges include loops; parallel pairs cancel.
"""
from fractions import Fraction as F
import itertools
import json
import math


def components(vertices, edges):
    seen = set()
    count = 0
    neighbors = [set() for _ in range(vertices)]
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    for root in range(vertices):
        if root in seen:
            continue
        count += 1
        stack = [root]
        seen.add(root)
        while stack:
            for nxt in neighbors[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
    return count


def b_value(vertices, edges):
    return F(vertices - components(vertices, edges)) - F(len(edges), 2)


def quotient(vertices, edges, labels):
    new_edges = set()
    for u, v in edges:
        edge = tuple(sorted((labels[u], labels[v])))
        if edge in new_edges:
            new_edges.remove(edge)
        else:
            new_edges.add(edge)
    return max(labels) + 1, new_edges


def reduce_graph(vertices, edges):
    initial_b = b_value(vertices, edges)
    initial_nonempty = bool(edges)
    steps = 0
    while True:
        candidate = None
        for w in range(vertices):
            adjacent = [v if u == w else u for u, v in edges if w in (u, v)]
            if (w, w) not in edges and len(adjacent) == 2:
                candidate = (w, *adjacent)
                break
        if candidate is None:
            break
        w, u, v = candidate
        assert u != v
        old_b = b_value(vertices, edges)
        old_count = components(vertices, edges)
        keep_edges = {edge for edge in edges if w not in edge}
        groups = {}
        labels = []
        for vertex in range(vertices):
            if vertex == w:
                labels.append(-1)
            else:
                representative = u if vertex == v else vertex
                if representative not in groups:
                    groups[representative] = len(groups)
                labels.append(groups[representative])
        new_vertices, new_edges = quotient(vertices, keep_edges, labels)
        assert new_vertices == vertices - 2
        cancellations = (len(edges) - 2 - len(new_edges)) // 2
        assert len(edges) - 2 - len(new_edges) == 2 * cancellations
        assert components(new_vertices, new_edges) - old_count <= cancellations
        assert old_b - b_value(new_vertices, new_edges) <= 1
        vertices, edges = new_vertices, new_edges
        steps += 1
    assert b_value(vertices, edges) <= 0
    assert steps >= initial_b
    if initial_nonempty and initial_b >= 0:
        assert steps >= 1
    return vertices, edges, steps


def hom(h, vertices, edges, injective=False):
    size = len(h)
    assignments = itertools.permutations(range(size), vertices) if injective else itertools.product(range(size), repeat=vertices)
    total = 0
    count = 0
    for colors in assignments:
        total += math.prod(h[colors[u]][colors[v]] for u, v in edges)
        count += 1
    return F(total, count)


def partitions(n):
    def visit(labels):
        if len(labels) == n:
            yield tuple(labels)
        else:
            for value in range(max(labels) + 2):
                yield from visit(labels + [value])
    yield from visit([0])


def main():
    structural = 0
    examples = []
    for vertices in range(1, 6):
        possible = list(itertools.combinations(range(vertices), 2))
        for code in range(1 << len(possible)):
            nonloops = {edge for index, edge in enumerate(possible) if code & (1 << index)}
            if any(sum(vertex in edge for edge in nonloops) % 2 for vertex in range(vertices)):
                continue
            for loop_code in range(1 << vertices):
                edges = nonloops | {(v, v) for v in range(vertices) if loop_code & (1 << v)}
                reduced_vertices, reduced_edges, steps = reduce_graph(vertices, edges)
                structural += 1
                if vertices <= 4:
                    examples.append((vertices, edges, reduced_vertices, reduced_edges, steps))
    h = [[1 - 2 * (bin(i & j).count('1') % 2) for j in range(4)] for i in range(4)]
    exact_hom = 0
    exact_mobius = 0
    for vertices, edges, reduced_vertices, reduced_edges, steps in examples:
        original = hom(h, vertices, edges)
        assert original == F(1, 4 ** steps) * hom(h, reduced_vertices, reduced_edges)
        exact_hom += 1
        mobius_sum = F(0)
        for labels in partitions(vertices):
            block_count, new_edges = quotient(vertices, edges, labels)
            sizes = [labels.count(index) for index in range(block_count)]
            coefficient = math.prod((-1) ** (size - 1) * math.factorial(size - 1) for size in sizes)
            mobius_sum += coefficient * F(1, 4 ** (vertices - block_count)) * hom(h, block_count, new_edges)
        prefactor = F(4 ** vertices, math.factorial(4) // math.factorial(4 - vertices))
        assert hom(h, vertices, edges, injective=True) == prefactor * mobius_sum
        exact_mobius += 1
    print(json.dumps({'eulerian_loop_graphs_structurally_checked': structural,
                      'exact_degree_two_homomorphism_checks': exact_hom,
                      'exact_partition_mobius_checks': exact_mobius}, indent=2))


if __name__ == '__main__':
    main()
