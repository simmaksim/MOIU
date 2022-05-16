from time import sleep
import numpy as np


def N_W_corner(a, b):
    m = len(a)
    n = len(b)
    a_t, b_t = a[:], b[:]
    X = np.zeros((m, n))
    i, j = 0, 0
    B = []

    while i < m and j < n:
        s = min(a_t[i], b_t[j])
        a_t[i] -= s
        b_t[j] -= s
        X[i][j] = s
        # adding point (i, j) to B
        B.append((i, j))
        if (a_t[i] == 0):
            i += 1
        elif (b_t[j] == 0):
            j += 1
    # maybe go further idk
    if i == m:
        i -= 1
        j += 1
        while j < n:
            B.append((i, j))
            j += 1
        j -= 1
    if j == n:
        i -= 1
        while i < m:
            B.append((i, j))
            j += 1
    # print(f"X:\n{X}\nB={B}")
    return X, B

 # выделить угловые ячейки
def Corners(B: set, m, n) -> set:
    has_removable = True
    while has_removable:
        has_removable = False
        for i in range(m):
            count = 0
            j_rem = -1
            for j in range(n):
                if (i, j) in B:
                    j_rem = j
                    count += 1
            if count == 1:
                B.remove((i, j_rem))
                has_removable = True
        for j in range(n):
            count = 0
            i_rem = -1
            for i in range(m):
                if (i, j) in B:
                    i_rem = i
                    count += 1
            if count == 1:
                B.remove((i_rem, j))
                has_removable = True
    return B


def ProsCos(B: set, m, n, i_fixed, j_fixed):
    B_new = [(i_fixed, j_fixed, "+")]
    B.remove((i_fixed, j_fixed))
    k = 0
    while True:
        if k == len(B_new):
            break
        b = B_new[k]
        k += 1
        # up
        for i in reversed(range(0, b[0])):
            if (i, b[1]) in B:
                B_new.append((i, b[1], "-" if b[2] == "+" else "+"))
                B.remove((i, b[1]))
                break
        # right
        for j in range(b[1] + 1, n):
            if (b[0], j) in B:
                B_new.append((b[0], j, "-" if b[2] == "+" else "+"))
                B.remove((b[0], j))
                break
        # down
        for i in range(b[0] + 1, m):
            if (i, b[1]) in B:
                B_new.append((i, b[1], "-" if b[2] == "+" else "+"))
                B.remove((i, b[1]))
                break
        # left
        for j in reversed(range(0, b[1])):
            if (b[0], j) in B:
                B_new.append((b[0], j, "-" if b[2] == "+" else "+"))
                B.remove((b[0], j))
                break
    print(B_new)
    return B_new



def SecP(a, b, C, X, B: list):
    m = len(a)
    n = len(b)

    flag = False
    count_It = 1
    while not flag:
        print(f"\nIteration: {count_It}")
        b_free = np.zeros((n + m,))
        A = np.zeros((n + m, n + m))
        k = 0
        for b in B:
            i, j = b[0], b[1]
            A[k][i], A[k][m + j] = 1, 1
            b_free[k] = C[i][j]
            k += 1
        A[n + m - 1][0] = 1
        sol = np.linalg.solve(A, b_free)
        u, v = sol[0: m], sol[m: ]
        flag = True
        i_fixed, j_fixed = 0, 0
        for i in range(m):
            if not flag:
                break
            for j in range(n):
                if u[i] + v[j] > C[i][j]:
                    i_fixed, j_fixed = i, j
                    flag = False
                    break
        if flag:
            break
        B.append((i_fixed, j_fixed))
        print(B)

        B_subset = Corners(set(B), m, n)
        B_subset = ProsCos(B_subset, m, n, i_fixed, j_fixed)
        XB_minus_subset = []
        for b in B_subset:
            if b[2] == "-":
                XB_minus_subset.append(X[b[0]][b[1]])
        Q = min(XB_minus_subset)
        
        for b in B_subset:
            i, j = b[0], b[1]
            X[i][j] += Q if b[2] == "+" else -Q
        
        for b in B_subset:
            i, j = b[0], b[1]
            if X[i][j] == 0 and b[2] == "-":
                B.remove((i, j))
                break
        print("X:\n", X)
        count_It += 1
        sleep(1)
    return X


def main():

    a = [0, 0, 0]
    b = [0, 0, 0]
    c = [[8,4,1],[8,4,3], [9,7,5]]

    X, B = N_W_corner(a, b)
    X = SecP(a, b, c, X, B)
    print(f"Optimal X:\n{X}")


if __name__ == "__main__":
    main()