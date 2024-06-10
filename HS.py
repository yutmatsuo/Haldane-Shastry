import numpy as np
from scipy.linalg import eigh

def spin(a, n):
    s = np.zeros(n, dtype=np.int8)
    for i in range(n):
        s[i] = (a % (2 ** (i+1)) - a % (2 ** i )) // (2 ** i)
    return s

def hs(n):
    h = np.zeros((2 ** n, 2 ** n))
    for a in range(2 ** n):
        for b in range(2 ** n):
            spa = spin(a, n)
            spb = spin(b, n)

            for i in range(1, n):
                for j in range(i):
                    if spa[i] == spb[i] and spa[j] == spb[j]:
                        h[a, b] += 1 / (4 * np.sin(np.pi * (i - j) / n) ** 2)
                    elif spa[i] == spb[j] and spa[j] == spb[i]:
                        h[a, b] -= 1 / (4 * np.sin(np.pi * (i - j) / n) ** 2)
    return h

def hse(n):
    return eigh(hs(n))[0]

