from random import Random

from search_path_ratio_r7 import evaluate, random_matrix


def flip(A, i, j):
    B = [list(row) for row in A]
    B[i][j] *= -1
    B[j][i] *= -1
    return tuple(tuple(row) for row in B)


def score(A):
    P, N, V = evaluate(A)
    return V / (P + N), P, N, V


def main(n=8, restarts=100, seed=24680):
    rng = Random(seed)
    global_best = (10,)
    for restart in range(restarts):
        A = random_matrix(n, rng)
        cur = score(A)
        changed = True
        while changed:
            changed = False
            candidates = []
            for i in range(n):
                for j in range(i + 1, n):
                    B = flip(A, i, j)
                    row = score(B)
                    candidates.append((row, rng.random(), B))
            row, _, B = min(candidates, key=lambda x: (x[0][0], x[0][3], x[1]))
            if (row[0], row[3]) < (cur[0], cur[3]):
                A, cur, changed = B, row, True
        if cur[0] < global_best[0]:
            global_best = (*cur, A)
            print("new", restart, cur, flush=True)
    print("FINAL", global_best[:4])
    for row in global_best[4]:
        print(row)


if __name__ == "__main__":
    main()
