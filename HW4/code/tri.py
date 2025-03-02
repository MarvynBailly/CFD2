import numpy as np

def solve_tridiagonal(a, b, c, f):
    """
    Thomas algorithm to solve tridiagonal system.
    Input arrays a, b, c, and f are 1D NumPy arrays.
    This function modifies f in place and returns the solution.
    (Based on tri.m)
    """
    M = len(a)
    x = np.zeros(M)
    # forward sweep
    x[0] = c[0] / b[0]
    f[0] = f[0] / b[0]
    for j in range(1, M):
        z = 1.0 / (b[j] - a[j] * x[j-1])
        x[j] = c[j] * z
        f[j] = (f[j] - a[j] * f[j-1]) * z
    # backward sweep
    for j in range(M-2, -1, -1):
        f[j] = f[j] - x[j] * f[j+1]
    return f